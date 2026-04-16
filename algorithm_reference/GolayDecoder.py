# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

"""
module provides golay encoding/decoding of DNA barcodes.
may not be the same golay code as previously used.
Provides mainly decode_nt()

If you wish to assign a read DNA seq barcode to a list of known originals,
correcting for errors if necessary, I recommend you use the generic version
in barcode.py .  Golay decoding assumes that the sequence can be any valid
golay sequence (4096 options), not just those you used in the study.

this module implements *a* golay (24,12,8) code.  That's 24 bit codewords,
12 bits of information, min 8 bits (hamming distance) between any 2 codewords
there are 2**12 = 4096 codewords, all <= 3 bit errors are corrected, 4 bit
errors are detected, and worse errors can cause erroneously corrected
codewords.

Since DNA has 4 nucleotides, each base represents 2 bits (see
DEFAULT_GOLAY_NT_TO_BITS).  The default values represent A <--> C and  G <--> T
as 2 bit errors, all other nucleotide swaps are 1 bit errors.
e.g.:
cw1 = AAACCCGGGTTT (24 bits, 12 nucleotides)
cw2 = AAACCCGGGTTA (bit distance of 2 from cw1)
cw3 = AAACCCGGGTTG (bit distance of 1 from cw1)

A specific golay code was chosen, as there are multiple.

bitvectors as referenced here are often just listlike objects with ints

refs:
http://www-math.mit.edu/phase2/UJM/vol1/MKANEM~1.PDF
http://cwww.ee.nctu.edu.tw/~ftchien/class/ecc08f/topic2.pdf
error correcting coding theory (Rhee)
the art of error correcting coding (Morelos-Zaragoza)
"""


import numpy as np
import math

class GolayDecoder:

    def __init__(self):
        # BEGIN module level constants
        self.NT_TO_BITS = {"A": "11", "C": "00", "T": "10", "G": "01"}

        # We use this matrix as the parity submatrix P
        self.DEFAULT_P = np.array([
            [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ],
            [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, ],
            [1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, ],
            [1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, ],
            [1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, ],
            [1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, ],
            [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, ],
            [1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, ],
            [1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, ],
            [1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, ],
            [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, ],
            [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, ]], dtype='int')
        
        # generator mtx G, where transmitted codewords (24bits) are
        # G.T dot msg or (msg dot G) (msg is 12 bit message)
        # 2**12 = 4096 total codewords, one for each msg
        # (all mod 2 arithmetic)
        # other G matrices give golay (24,12,8) codes, but this one
        # matches existing codes from pre 2010 used in knight lab
        self.DEFAULT_G = np.concatenate((self.DEFAULT_P,
                                         np.eye(12, dtype="int")), axis=1)

        # pairity check matrix H satisfies G dot H.T = zeros (mod 2 arithmetic)
        # also satisfies syn = H dot rec = H dot err (rec is recieved 24 bits,
        # err is 24 bit error string added to transmitted 24 bit vec)
        # (all mod 2 arithmetic)
        self.DEFAULT_H = np.concatenate((np.eye(12, dtype="int"),
                                         self.DEFAULT_P.T), axis=1)

        # len = 2325.  (1 (all zeros) + 24 (one 1) + 276 (two 1s) + 2024)
        self._ALL_3BIT_ERRORS = self._make_3bit_errors()

        # build syndrome lookup table
        self.DEFAULT_SYNDROME_LUT = {}
        for errvec in self._ALL_3BIT_ERRORS:
            syn = tuple(np.dot(self.DEFAULT_H, errvec) % 2)
            self.DEFAULT_SYNDROME_LUT[syn] = errvec

        self.BITS_TO_NT = {v: k for k, v in self.NT_TO_BITS.items()}
    
    def decode(self, seq):
        """Decode a nucleotide string of 12 bases, using bitwise error checking

        Parameters
        ----------
        seq : str
            The nucleotide string to decode

        Returns
        -------
        str, int
            The corrected barcode and the number of observed errors.

            If a 4 bit error is detected, the corrected barcode returned will
            be None.
        """

        if len(seq) != 12:
            raise ValueError("Golay decoding requires 12nt barcodes. The "
                             "barcode attempting to be decoded (%s) is of "
                             "length %dnt." % (seq, len(seq)))

        if not set(seq).issubset({'A', 'T', 'G', 'C'}):
            return None, 4

        received_bits = self._seq_to_bits(seq)
        corrected_bits, num_errors = self.decode_bits(received_bits)

        if corrected_bits is None:
            return None, num_errors
        else:
            return self._bits_to_seq(corrected_bits), num_errors
        
    def _seq_to_bits(self, seq):
        """Convert a nucleotide sequence to a bitvector

        Parameters
        ----------
        seq : str
            The nucleotide sequence

        Returns
        -------
        np.array
            The bit pattern of the nucleotide sequence
        """
        bitstring = list(''.join([self.NT_TO_BITS[nt] for nt in seq]))
        return np.asarray(bitstring, dtype=int)
    
    def decode_bits(self, received_bitvec):
        """Decode a received 24 bit vector to a corrected 24 bit vector

        Parameters
        ----------
        received_bitvec : np.array
            The bit vector to error correct.

        Returns
        -------
        np.array, int
            An array representing the corrected bit vector and the number of
            oobserved errors.

            If the number of errors is >= 4, a None is returned instead of the
            vector
        """
        syn = np.dot(self.DEFAULT_H, received_bitvec) % 2
        try:
            err = self.DEFAULT_SYNDROME_LUT[tuple(syn)]
        except KeyError:
            return None, 4
        corrected = (received_bitvec + err) % 2

        return corrected, np.sum(err)
    
    def _bits_to_seq(self, bits):
        """Convert a bit pattern to a sequence

        Parameters
        ----------
        bits : np.array
            The bit pattern

        Returns
        -------
        str
            The corresponding nucleotide sequence
        """
        seq = ""
        for i in range(0, len(bits), 2):  # take bits in twos
            bit1 = str(bits[i])
            bit2 = str(bits[i + 1])
            seq += self.BITS_TO_NT[bit1 + bit2]
        return seq
    
    def _make_3bit_errors(self, veclen=24):
        """Construct all bitvectors with <= 3 bits as 1's, rest 0's

        Parameters
        ----------
        veclen : int
            The bit length to construct

        Returns
        -------
        np.array
            The array of vectors with errors
        """
        def _comb(n, k):
            fac = math.factorial
            return fac(n) / fac(k) / fac(n - k)

        nvecs = int(veclen + _comb(veclen, 2) + _comb(veclen, 3))

        # +1 for an all zero vector
        errorvecs = np.zeros((nvecs + 1, veclen), dtype=int)

        # one 1
        offset = 1
        for i in range(veclen):
            vec = errorvecs[offset]
            vec[i] = 1
            offset += 1

        # two 1s
        for i in range(veclen):
            for j in range(i + 1, veclen):
                vec = errorvecs[offset]
                vec[i] = 1
                vec[j] = 1
                offset += 1

        # three 1s
        for i in range(veclen):
            for j in range(i + 1, veclen):
                for k in range(j + 1, veclen):
                    vec = errorvecs[offset]
                    vec[i] = 1
                    vec[j] = 1
                    vec[k] = 1
                    offset += 1

        return errorvecs
    
if __name__ == "__main__":
    golay = GolayDecoder()

    print("Decode from ATGCAGCTCAGT: ", golay.decode("ATGCAGCTCAGT"))
    print("Decode from CCCCTCAGCGGC: ", golay.decode("CCCCTCAGCGGC"))