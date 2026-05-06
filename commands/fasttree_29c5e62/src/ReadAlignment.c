#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"

alignment_t *ReadAlignment(/*IN*/FILE *fp, bool bQuote) {
	/* bQuote supports the -quote option */
	int nSeq = 0;
	int nPos = 0;
	char **names = NULL;
	char **seqs = NULL;
	char buf[BUFFER_SIZE] = "";
  	if (fgets(buf,sizeof(buf),fp) == NULL) {
    	fprintf(stderr, "Error reading header line\n");
    	exit(1);
  	}
  	int nSaved = 100;
  	if (buf[0] == '>') {
    	/* FASTA, truncate names at any of these */
    	char *nameStop = bQuote ? "'\t\r\n" : "(),: \t\r\n";
    	char *seqSkip = " \t\r\n";	/* skip these characters in the sequence */
    	seqs = (char**)mymalloc(sizeof(char*) * nSaved);
    	names = (char**)mymalloc(sizeof(char*) * nSaved);

		do {
			/* loop over lines */
			if (buf[0] == '>') {
				/* truncate the name */
				char *p, *q;
				for (p = buf+1; *p != '\0'; p++) {
				for (q = nameStop; *q != '\0'; q++) {
					if (*p == *q) {
					*p = '\0';
					break;
					}
				}
				if (*p == '\0') break;
				}

				/* allocate space for another sequence */
				nSeq++;
				if (nSeq > nSaved) {
					int nNewSaved = nSaved*2;
					seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved, /*copy*/false);
					names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved, /*copy*/false);
					nSaved = nNewSaved;
				}
				names[nSeq-1] = (char*)mymemdup(buf+1,strlen(buf));
				seqs[nSeq-1] = NULL;
			} else {
				/* count non-space characters and append to sequence */
				int nKeep = 0;
				char *p, *q;
				for (p=buf; *p != '\0'; p++) {
				for (q=seqSkip; *q != '\0'; q++) {
					if (*p == *q)
					break;
				}
				if (*p != *q)
					nKeep++;
				}
				int nOld = (seqs[nSeq-1] == NULL) ? 0 : strlen(seqs[nSeq-1]);
				seqs[nSeq-1] = (char*)myrealloc(seqs[nSeq-1], nOld, nOld+nKeep+1, /*copy*/false);
				if (nOld+nKeep > nPos)
				nPos = nOld + nKeep;
				char *out = seqs[nSeq-1] + nOld;
				for (p=buf; *p != '\0'; p++) {
					for (q=seqSkip; *q != '\0'; q++) {
						if (*p == *q)
						break;
					}
					if (*p != *q) {
						*out = *p;
						out++;
					}
				}
				assert(out-seqs[nSeq-1] == nKeep + nOld);
				*out = '\0';
			}
		} while(fgets(buf,sizeof(buf),fp) != NULL);

		if (seqs[nSeq-1] == NULL) {
			fprintf(stderr, "No sequence data for last entry %s\n",names[nSeq-1]);
			exit(1);
		}
		names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nSeq, /*copy*/false);
		seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nSeq, /*copy*/false);
	} else {
		/* PHYLIP interleaved-like format
		Allow arbitrary length names, require spaces between names and sequences
		Allow multiple alignments, either separated by a single empty line (e.g. seqboot output)
		or not.
		*/
    	if (buf[0] == '\n' || buf[0] == '\r') {
      		if (fgets(buf,sizeof(buf),fp) == NULL) {
				fprintf(stderr, "Empty header line followed by EOF\n");
				exit(1);
      		}
    	}
    	if (sscanf(buf, "%d%d", &nSeq, &nPos) != 2 || nSeq < 1 || nPos < 1) {
      		fprintf(stderr, "Error parsing header line:%s\n", buf);
      		exit(1);
   		}
    	names = (char **)mymalloc(sizeof(char*) * nSeq);
    	seqs = (char **)mymalloc(sizeof(char*) * nSeq);
    	nSaved = nSeq;

    	int i;
		for (i = 0; i < nSeq; i++) {
			names[i] = NULL;
			seqs[i] = (char *)mymalloc(nPos+1);	/* null-terminate */
			seqs[i][0] = '\0';
		}
		int iSeq = 0;
    
    	while(fgets(buf,sizeof(buf),fp)) {
      		if ((buf[0] == '\n' || buf[0] == '\r') && (iSeq == nSeq || iSeq == 0)) {
				iSeq = 0;
			} else {
				int j = 0; /* character just past end of name */
				if (buf[0] == ' ') {
					if (names[iSeq] == NULL) {
						fprintf(stderr, "No name in phylip line %s", buf);
						exit(1);
					}
				} else {
					while (buf[j] != '\n' && buf[j] != '\0' && buf[j] != ' ')
						j++;
					if (buf[j] != ' ' || j == 0) {
						fprintf(stderr, "No sequence in phylip line %s", buf);
						exit(1);
					}
					if (iSeq >= nSeq) {
						fprintf(stderr, "No empty line between sequence blocks (is the sequence count wrong?)\n");
						exit(1);
					}
					if (names[iSeq] == NULL) {
						/* save the name */
						names[iSeq] = (char *)mymalloc(j+1);
						int k;
						for (k = 0; k < j; k++) names[iSeq][k] = buf[k];
						names[iSeq][j] = '\0';
					} else {
						/* check the name */
						int k;
						int match = 1;
						for (k = 0; k < j; k++) {
							if (names[iSeq][k] != buf[k]) {
								match = 0;
								break;
							}
						}
						if (!match || names[iSeq][j] != '\0') {
							fprintf(stderr, "Wrong name in phylip line %s\nExpected %s\n", buf, names[iSeq]);
							exit(1);
						}
					}
				}
				int seqlen = strlen(seqs[iSeq]);
				for (; buf[j] != '\n' && buf[j] != '\0'; j++) {
					if (buf[j] != ' ') {
						if (seqlen >= nPos) {
							fprintf(stderr, "Too many characters (expected %d) for sequence named %s\nSo far have:\n%s\n",
							nPos, names[iSeq], seqs[iSeq]);
							exit(1);
						}
						seqs[iSeq][seqlen++] = toupper(buf[j]);
					}
				}
				seqs[iSeq][seqlen] = '\0'; /* null-terminate */
				if(verbose>10) fprintf(stderr,"Read iSeq %d name %s seqsofar %s\n", iSeq, names[iSeq], seqs[iSeq]);
				iSeq++;
				if (iSeq == nSeq && strlen(seqs[0]) == nPos)
					break; /* finished alignment */
			} /* end else non-empty phylip line */
    	}
		if (iSeq != nSeq && iSeq != 0) {
		fprintf(stderr, "Wrong number of sequences: expected %d\n", nSeq);
		exit(1);
		}
  	}
  	/* Check lengths of sequences */
  	int i;
	for (i = 0; i < nSeq; i++) {
		int seqlen = strlen(seqs[i]);
		if (seqlen != nPos) {
		fprintf(stderr, "Wrong number of characters for %s: expected %d but have %d instead.\n"
			"This sequence may be truncated, or another sequence may be too long.\n",
			names[i], nPos, seqlen);
		exit(1);
		}
	}
  	/* Replace "." with "-" and warn if we find any */
  	/* If nucleotide sequences, replace U with T and N with X */
  	bool findDot = false;
	for (i = 0; i < nSeq; i++) {
		char *p;
		for (p = seqs[i]; *p != '\0'; p++) {
			if (*p == '.') {
				findDot = true;
				*p = '-';
			}
			if (nCodes == 4 && *p == 'U')
				*p = 'T';
			if (nCodes == 4 && *p == 'N')
				*p = 'X';
		}
	}
	if (findDot)
		fprintf(stderr, "Warning! Found \".\" character(s). These are treated as gaps\n");

	if (ferror(fp)) {
		fprintf(stderr, "Error reading input file\n");
		exit(1);
	}

	alignment_t *align = (alignment_t*)mymalloc(sizeof(alignment_t));
	align->nSeq = nSeq;
	align->nPos = nPos;
	align->names = names;
	align->seqs = seqs;
	align->nSaved = nSaved;
	return(align);
}