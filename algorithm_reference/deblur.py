import numpy as np
from collections import defaultdict
from operator import attrgetter
import re
import io
import skbio

sniff_fasta = skbio.io.io_registry.get_sniffer('fasta')
sniff_fastq = skbio.io.io_registry.get_sniffer('fastq')

def _get_fastq_variant(input_fp):
    # https://bit.ly/3GEDIxF
    variant = None
    variants = ['illumina1.8', 'illumina1.3', 'solexa', 'sanger']
    for v in variants:
        try:
            next(skbio.read(input_fp, format='fastq', variant=v))
        except Exception:
            continue
        else:
            variant = v
            break

    if variant is None:
        raise ValueError("Unknown variant, unable to interpret PHRED")

    return variant

def sequence_generator(input_fp):
    """Yield (id, sequence) from an input file

    Parameters
    ----------
    input_fp : filepath
        A filepath, which can be any valid fasta or fastq file within the
        limitations of scikit-bio's IO registry.

    Notes
    -----
    The use of this method is a stopgap to replicate the existing `parse_fasta`
    functionality while at the same time allowing for fastq support.

    Raises
    ------
    skbio.io.FormatIdentificationWarning
        If the format of the input file cannot be determined.

    Returns
    -------
    (str, str)
        The ID and sequence.

    """
    kw = {}
    if sniff_fasta(input_fp)[0]:
        format = 'fasta'
    elif sniff_fastq(input_fp)[0]:
        format = 'fastq'

        kw['variant'] = _get_fastq_variant(input_fp)
    else:
        # usually happens when the fasta file is empty
        # so need to return no sequences (and warn)
        msg = "input file %s does not appear to be FASTA or FASTQ" % input_fp
        print(msg)
        return

    # some of the test code is using file paths, some is using StringIO.
    if isinstance(input_fp, io.TextIOBase):
        input_fp.seek(0)

    for record in skbio.read(input_fp, format=format, **kw):
        yield (record.metadata['id'], str(record))

class Sequence(object):
    """Sequence object to represent the aligned reads

    Attributes
    ----------
    label : str
        The sequence label
    sequence : str
        The sequence string
    length : int
        The sequence length (aligned)
    unaligned_length : int
        The unaligned sequence length
    frequency : float
        The number of times the sequence have been seen in the dataset
    np_sequence : numpy array of int8
        An int8 numpy representation of the sequence string

    Methods
    -------
    to_fasta
    """

    def __init__(self, label, sequence):
        self.label = label
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
        self.unaligned_length = self.length - self.sequence.count('-')
        self.frequency = float(
            re.search(r'(?<=size=)\w+', self.label).group(0))
        
        trans_dict = defaultdict(lambda: 5, A=0, C=1, G=2, T=3)
        trans_dict['-'] = 4
        self.np_sequence = np.array(
            [trans_dict[b] for b in self.sequence], dtype=np.int8)

    def __eq__(self, other):
        return (type(self) == type(other) and
                self.sequence == other.sequence and
                self.frequency == other.frequency)

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_fasta(self):
        """Returns a string with the sequence in fasta format

        Returns
        -------
        str
            The FASTA representation of the sequence
        """
        prefix, suffix = re.split(r'(?<=size=)\w+', self.label, maxsplit=1)
        new_count = int(round(self.frequency))
        new_label = "%s%d%s" % (prefix, new_count, suffix)
        return ">%s\n%s\n" % (new_label, self.sequence)

def get_default_error_profile():
    """Return the default error profile for deblurring
    based on illumina run data
    """
    error_dist = [1, 0.06, 0.02, 0.02, 0.01,
                  0.005, 0.005, 0.005, 0.001, 0.001,
                  0.001, 0.0005]
    return error_dist

def get_sequences(input_seqs):
    """Returns a list of Sequences

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format

    Returns
    -------
    list of Sequence

    Raises
    ------
    ValueError
        If no sequences where found in `input_seqs`
        If all the sequences do not have the same length either aligned or
        unaligned.
    """
    try:
        seqs = [Sequence(id, seq) for id, seq in input_seqs]
    except Exception:
        seqs = []

    if len(seqs) == 0:
        print('No sequences found in fasta file!')
        return None

    # Check that all the sequence lengths (aligned and unaligned are the same)
    aligned_lengths = set(s.length for s in seqs)
    unaligned_lengths = set(s.unaligned_length for s in seqs)

    if len(aligned_lengths) != 1 or len(unaligned_lengths) != 1:
        raise ValueError(
            "Not all sequence have the same length. Aligned lengths: %s, "
            "sequence lengths: %s"
            % (", ".join(map(str, aligned_lengths)),
               ", ".join(map(str, unaligned_lengths))))

    seqs = sorted(seqs, key=attrgetter('frequency'), reverse=True)
    return seqs


# Copy from: /opt/conda/envs/qiime2-amplicon-2026.1/lib/python3.10/site-packages/deblur/deblurring.py:71
def deblur(input_seqs, mean_error=0.005,
           error_dist=None,
           indel_prob=0.01, indel_max=3):
    """Deblur the reads

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format. The label
        should include the sequence count in the 'size=X' format.
    mean_error : float, optional
        The mean illumina error, used for original sequence estimate.
        Default: 0.005
    error_dist : list of float, optional
        A list of error probabilities. The length of the list determines the
        amount of hamming distances taken into account. Default: None, use
        the default error profile (from get_default_error_profile() )
    indel_prob : float, optional
        Indel probability (same for N indels). Default: 0.01
    indel_max : int, optional
        The maximal number of indels expected by errors. Default: 3

    Results
    -------
    list of Sequence
        The deblurred sequences

    Notes
    -----
    mean_error is used only for normalizing the peak height before deblurring.
    The array 'error_dist' represents the error distribution, where
    Xi = max frequency of error hamming. The length of this array - 1 limits
    the hamming distance taken into account, i.e. if the length if `error_dist`
    is 10, sequences up to 10 - 1 = 9 hamming distance will be taken into
    account
    """

    if error_dist is None:
        error_dist = get_default_error_profile()
    print('Using error profile %s' % error_dist)

    # Get the sequences
    seqs = get_sequences(input_seqs)
    if seqs is None:
        print('no sequences deblurred')
        return None
    print('deblurring %d sequences' % len(seqs))

    # fix the original frequencies of each read error using the
    # mean error profile
    mod_factor = pow((1 - mean_error), seqs[0].unaligned_length)
    error_dist = np.array(error_dist) / mod_factor

    max_h_dist = len(error_dist) - 1

    for seq_i in seqs:
        # no need to remove neighbors if freq. is <=0
        if seq_i.frequency <= 0:
            continue

        # Correct for the fact that many reads are expected to be mutated
        num_err = error_dist * seq_i.frequency

        # if it's low level, just continue
        if num_err[1] < 0.1:
            continue

        # Compare to all other sequences and calculate hamming dist
        seq_i_len = len(seq_i.sequence.rstrip('-'))
        for seq_j in seqs:
            # Ignore current sequence
            if seq_i == seq_j:
                continue

            # Calculate the hamming distance
            h_dist = np.count_nonzero(np.not_equal(seq_i.np_sequence,
                                                   seq_j.np_sequence))

            # If far away, don't need to correct
            if h_dist > max_h_dist:
                continue

            # Close, so lets calculate exact distance

            # We stop checking in the shortest sequence after removing trailing
            # indels. We need to do this in order to avoid double counting
            # the insertions/deletions
            length = min(seq_i_len, len(seq_j.sequence.rstrip('-')))
            sub_seq_i = seq_i.np_sequence[:length]
            sub_seq_j = seq_j.np_sequence[:length]

            mask = (sub_seq_i != sub_seq_j)
            # find all indels
            mut_is_indel = np.logical_or(sub_seq_i[mask] == 4,
                                         sub_seq_j[mask] == 4)
            num_indels = mut_is_indel.sum()
            if num_indels > 0:
                # need to account for indel in one sequence not solved in the
                # other (so we have '-' at the end. Need to ignore it in the
                # total count)
                h_dist = np.count_nonzero(
                    np.not_equal(seq_i.np_sequence[:length],
                                 seq_j.np_sequence[:length]))

            num_substitutions = h_dist - num_indels

            correction_value = num_err[num_substitutions]

            if num_indels > indel_max:
                correction_value = 0
            elif num_indels > 0:
                # remove errors due to (PCR?) indels (saw in 22 mock mixture)
                correction_value = correction_value * indel_prob

            # met all the criteria - so correct the frequency of the neighbor
            seq_j.frequency -= correction_value

    result = [s for s in seqs if round(s.frequency) > 0]
    print('%d unique sequences left following deblurring' % len(result))
    return result

if __name__ == "__main__":
    ## Outputs from MSA
    alignment_output = "/workspaces/qiime2_blog/support_data/L6S68_30_L001_R1_001.fastq.gz.trim.derep.no_artifacts.msa"
    msa = sequence_generator(alignment_output)

    ## Parameters
    mean_error = 0.005
    error_dist = [1.0, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005]
    indel_prob = 0.01
    indel_max = 3


    ## Dbluring
    seqs = deblur(msa, mean_error, error_dist, indel_prob, indel_max)