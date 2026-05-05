#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

char *usage =
  "  FastTree protein_alignment > tree\n"
  "  FastTree < protein_alignment > tree\n"
  "  FastTree -out tree protein_alignment\n"
  "  FastTree -nt nucleotide_alignment > tree\n"
  "  FastTree -nt -gtr < nucleotide_alignment > tree\n"
  "  FastTree < nucleotide_alignment > tree\n"
  "FastTree accepts alignments in fasta or phylip interleaved formats\n"
  "\n"
  "Common options (must be before the alignment file):\n"
  "  -quiet to suppress reporting information\n"
  "  -nopr to suppress progress indicator\n"
  "  -log logfile -- save intermediate trees, settings, and model details\n"
  "  -fastest -- speed up the neighbor joining phase & reduce memory usage\n"
  "        (recommended for >50,000 sequences)\n"
  "  -n <number> to analyze multiple alignments (phylip format only)\n"
  "        (use for global bootstrap, with seqboot and CompareToBootstrap.pl)\n"
  "  -nosupport to not compute support values\n"
  "  -intree newick_file to set the starting tree(s)\n"
  "  -intree1 newick_file to use this starting tree for all the alignments\n"
  "        (for faster global bootstrap on huge alignments)\n"
  "  -pseudo to use pseudocounts (recommended for highly gapped sequences)\n"
  "  -gtr -- generalized time-reversible model (nucleotide alignments only)\n"
  "  -lg -- Le-Gascuel 2008 model (amino acid alignments only)\n"
  "  -wag -- Whelan-And-Goldman 2001 model (amino acid alignments only)\n"
  "  -quote -- allow spaces and other restricted characters (but not ' ) in\n"
  "           sequence names and quote names in the output tree (fasta input only;\n"
  "           FastTree will not be able to read these trees back in)\n"
  "  -noml to turn off maximum-likelihood\n"
  "  -nome to turn off minimum-evolution NNIs and SPRs\n"
  "        (recommended if running additional ML NNIs with -intree)\n"
  "  -nome -mllen with -intree to optimize branch lengths for a fixed topology\n"
  "  -cat # to specify the number of rate categories of sites (default 20)\n"
  "      or -nocat to use constant rates\n"
  "  -gamma -- after optimizing the tree under the CAT approximation,\n"
  "      rescale the lengths to optimize the Gamma20 likelihood\n"
  "  -constraints constraintAlignment to constrain the topology search\n"
  "       constraintAlignment should have 1s or 0s to indicates splits\n"
  "  -expert -- see more options\n"
  "For more information, see http://www.microbesonline.org/fasttree/\n";

char *expertUsage =
  "FastTree [-nt] [-n 100] [-quote] [-pseudo | -pseudo 1.0]\n"
  "           [-boot 1000 | -nosupport]\n"
  "           [-intree starting_trees_file | -intree1 starting_tree_file]\n"
  "           [-quiet | -nopr]\n"
  "           [-nni 10] [-spr 2] [-noml | -mllen | -mlnni 10]\n"
  "           [-mlacc 2] [-cat 20 | -nocat] [-gamma]\n"
  "           [-slow | -fastest] [-2nd | -no2nd] [-slownni] [-seed 1253] \n"
  "           [-top | -notop] [-topm 1.0 [-close 0.75] [-refresh 0.8]]\n"
  "           [-gtr] [-gtrrates ac ag at cg ct gt] [-gtrfreq A C G T]\n"
  "           [ -lg | -wag | -trans transitionmatrixfile ]\n"
  "           [-matrix Matrix | -nomatrix] [-nj | -bionj]\n"
  "           [ -constraints constraintAlignment [ -constraintWeight 100.0 ] ]\n"
  "           [-log logfile]\n"
  "         [ alignment_file ]\n"
  "        [ -out output_newick_file | > newick_tree]\n"
  "\n"
  "or\n"
  "\n"
  "FastTree [-nt] [-matrix Matrix | -nomatrix] [-rawdist] -makematrix [alignment]\n"
  "    [-n 100] > phylip_distance_matrix\n"
  "\n"
  "  FastTree supports fasta or phylip interleaved alignments\n"
  "  By default FastTree expects protein alignments,  use -nt for nucleotides\n"
  "  FastTree reads standard input if no alignment file is given\n"
  "\n"
  "Input/output options:\n"
  "  -n -- read in multiple alignments in. This only\n"
  "    works with phylip interleaved format. For example, you can\n"
  "    use it with the output from phylip's seqboot. If you use -n, FastTree\n"
  "    will write 1 tree per line to standard output.\n"
  "  -intree newickfile -- read the starting tree in from newickfile.\n"
  "     Any branch lengths in the starting trees are ignored.\n"
  "    -intree with -n will read a separate starting tree for each alignment.\n"
  "  -intree1 newickfile -- read the same starting tree for each alignment\n"
  "  -quiet -- do not write to standard error during normal operation (no progress\n"
  "     indicator, no options summary, no likelihood values, etc.)\n"
  "  -nopr -- do not write the progress indicator to stderr\n"
  "  -log logfile -- save intermediate trees so you can extract\n"
  "    the trees and restart long-running jobs if they crash\n"
  "    -log also reports the per-site rates (1 means slowest category)\n"
  "  -quote -- quote sequence names in the output and allow spaces, commas,\n"
  "    parentheses, and colons in them but not ' characters (fasta files only)\n"
  "\n"
  "Distances:\n"
  "  Default: For protein sequences, log-corrected distances and an\n"
  "     amino acid dissimilarity matrix derived from BLOSUM45\n"
  "  or for nucleotide sequences, Jukes-Cantor distances\n"
  "  To specify a different matrix, use -matrix FilePrefix or -nomatrix\n"
  "  Use -rawdist to turn the log-correction off\n"
  "  or to use %different instead of Jukes-Cantor\n"
  "  (These options affect minimum-evolution computations only;\n"
  "   use -trans to affect maximum-likelihoood computations)\n"
  "\n"
  "  -pseudo [weight] -- Use pseudocounts to estimate distances between\n"
  "      sequences with little or no overlap. (Off by default.) Recommended\n"
  "      if analyzing the alignment has sequences with little or no overlap.\n"
  "      If the weight is not specified, it is 1.0\n"
  "\n"
  "Topology refinement:\n"
  "  By default, FastTree tries to improve the tree with up to 4*log2(N)\n"
  "  rounds of minimum-evolution nearest-neighbor interchanges (NNI),\n"
  "  where N is the number of unique sequences, 2 rounds of\n"
  "  subtree-prune-regraft (SPR) moves (also min. evo.), and\n"
  "  up to 2*log(N) rounds of maximum-likelihood NNIs.\n"
  "  Use -nni to set the number of rounds of min. evo. NNIs,\n"
  "  and -spr to set the rounds of SPRs.\n"
  "  Use -noml to turn off both min-evo NNIs and SPRs (useful if refining\n"
  "       an approximately maximum-likelihood tree with further NNIs)\n"
  "  Use -sprlength set the maximum length of a SPR move (default 10)\n"
  "  Use -mlnni to set the number of rounds of maximum-likelihood NNIs\n"
  "  Use -mlacc 2 or -mlacc 3 to always optimize all 5 branches at each NNI,\n"
  "      and to optimize all 5 branches in 2 or 3 rounds\n"
  "  Use -mllen to optimize branch lengths without ML NNIs\n"
  "  Use -mllen -nome with -intree to optimize branch lengths on a fixed topology\n"
  "  Use -slownni to turn off heuristics to avoid constant subtrees (affects both\n"
  "       ML and ME NNIs)\n"
  "\n"
  "Maximum likelihood model options:\n"
  "  -lg -- Le-Gascuel 2008 model instead of (default) Jones-Taylor-Thorton 1992 model (a.a. only)\n"
  "  -wag -- Whelan-And-Goldman 2001 model instead of (default) Jones-Taylor-Thorton 1992 model (a.a. only)\n"
  "  -gtr -- generalized time-reversible instead of (default) Jukes-Cantor (nt only)\n"
  "  -cat # -- specify the number of rate categories of sites (default 20)\n"
  "  -nocat -- no CAT model (just 1 category)\n"
  " - trans filename -- use the transition matrix from filename\n"
  "      This is supported for amino acid alignments only\n"
  "      The file must be tab-delimited with columns in the order ARNDCQEGHILKMFPSTWYV*\n"
  "      The additional column named * is for the stationary distribution\n"
  "      Each row must have a row name in the same order ARNDCQEGHILKMFPSTWYV\n"
  "  -gamma -- after the final round of optimizing branch lengths with the CAT model,\n"
  "            report the likelihood under the discrete gamma model with the same\n"
  "            number of categories. FastTree uses the same branch lengths but\n"
  "            optimizes the gamma shape parameter and the scale of the lengths.\n"
  "            The final tree will have rescaled lengths. Used with -log, this\n"
  "            also generates per-site likelihoods for use with CONSEL, see\n"
  "            GammaLogToPaup.pl and documentation on the FastTree web site.\n"
  "\n"
  "Support value options:\n"
  "  By default, FastTree computes local support values by resampling the site\n"
  "  likelihoods 1,000 times and the Shimodaira Hasegawa test. If you specify -noml,\n"
  "  it will compute minimum-evolution bootstrap supports instead\n"
  "  In either case, the support values are proportions ranging from 0 to 1\n"
  "\n"
  "  Use -nosupport to turn off support values or -boot 100 to use just 100 resamples\n"
  "  Use -seed to initialize the random number generator\n"
  "\n"
  "Searching for the best join:\n"
  "  By default, FastTree combines the 'visible set' of fast neighbor-joining with\n"
  "      local hill-climbing as in relaxed neighbor-joining\n"
  "  -slow -- exhaustive search (like NJ or BIONJ, but different gap handling)\n"
  "      -slow takes half an hour instead of 8 seconds for 1,250 proteins\n"
  "  -fastest -- search the visible set (the top hit for each node) only\n"
  "      Unlike the original fast neighbor-joining, -fastest updates visible(C)\n"
  "      after joining A and B if join(AB,C) is better than join(C,visible(C))\n"
  "      -fastest also updates out-distances in a very lazy way,\n"
  "      -fastest sets -2nd on as well, use -fastest -no2nd to avoid this\n"
  "\n"
  "Top-hit heuristics:\n"
  "  By default, FastTree uses a top-hit list to speed up search\n"
  "  Use -notop (or -slow) to turn this feature off\n"
  "         and compare all leaves to each other,\n"
  "         and all new joined nodes to each other\n"
  "  -topm 1.0 -- set the top-hit list size to parameter*sqrt(N)\n"
  "         FastTree estimates the top m hits of a leaf from the\n"
  "         top 2*m hits of a 'close' neighbor, where close is\n"
  "         defined as d(seed,close) < 0.75 * d(seed, hit of rank 2*m),\n"
  "         and updates the top-hits as joins proceed\n"
  "  -close 0.75 -- modify the close heuristic, lower is more conservative\n"
  "  -refresh 0.8 -- compare a joined node to all other nodes if its\n"
  "         top-hit list is less than 80% of the desired length,\n"
  "         or if the age of the top-hit list is log2(m) or greater\n"
  "   -2nd or -no2nd to turn 2nd-level top hits heuristic on or off\n"
  "      This reduces memory usage and running time but may lead to\n"
  "      marginal reductions in tree quality.\n"
  "      (By default, -fastest turns on -2nd.)\n"
  "\n"
  "Join options:\n"
  "  -nj: regular (unweighted) neighbor-joining (default)\n"
  "  -bionj: weighted joins as in BIONJ\n"
  "          FastTree will also weight joins during NNIs\n"
  "\n"
  "Constrained topology search options:\n"
  "  -constraints alignmentfile -- an alignment with values of 0, 1, and -\n"
  "       Not all sequences need be present. A column of 0s and 1s defines a\n"
  "       constrained split. Some constraints may be violated\n"
  "       (see 'violating constraints:' in standard error).\n"
  "  -constraintWeight -- how strongly to weight the constraints. A value of 1\n"
  "       means a penalty of 1 in tree length for violating a constraint\n"
  "       Default: 100.0\n"
  "\n"
  "For more information, see http://www.microbesonline.org/fasttree/\n"
  "   or the comments in the source code\n";

void InitOptions(FastTreeOptions_t *opt) {
  opt->nAlign = 1;
  opt->matrixPrefix = NULL;
  opt->transitionFile = NULL;
  opt->make_matrix = false;
  opt->constraintsFile = NULL;
  opt->intreeFile = NULL;
  opt->intree1 = false;
  opt->nni = -1;
  opt->spr = 2;
  opt->maxSPRLength = 10;
  opt->MLnni = -1;
  opt->MLlen = false;
  opt->nBootstrap = 1000;
  opt->nRateCats = nDefaultRateCats;
  opt->logfile = NULL;
  opt->bUseGtr = false;
  opt->bUseLg = false;
  opt->bUseWag = false;
  opt->bUseGtrRates = false;
  int i;
  for (i = 0; i < 6; i++) opt->gtrrates[i] = 1.0;
  opt->bUseGtrFreq = false;
  for (i = 0; i < 4; i++) opt->gtrfreq[i] = 0.25;
  opt->bQuote = false;
  opt->fpOut = stdout;
  opt->fileName = NULL;
}

void ParseCommandLine(int argc, char **argv, FastTreeOptions_t *opt) {
  int iArg;

  if (isatty(STDIN_FILENO) && argc == 1) {
    fprintf(stderr,"Usage for FastTree version %s %s%s:\n%s",
	    FT_VERSION, SSE_STRING, OpenMPString(), usage);
#if (defined _WIN32 || defined WIN32 || defined WIN64 || defined _WIN64)
    fprintf(stderr, "Windows users: Please remember to run this inside a command shell\n");
    fprintf(stderr,"Hit return to continue\n");
    fgetc(stdin);
#endif
    exit(0);
  }    

  for (iArg = 1; iArg < argc; iArg++) {
    if (strcmp(argv[iArg],"-makematrix") == 0) {
      opt->make_matrix = true;
    } else if (strcmp(argv[iArg],"-logdist") == 0) {
      fprintf(stderr, "Warning: logdist is now on by default and obsolete\n");
    } else if (strcmp(argv[iArg],"-rawdist") == 0) {
      logdist = false;
    } else if (strcmp(argv[iArg],"-verbose") == 0 && iArg < argc-1) {
      verbose = atoi(argv[++iArg]);
    } else if (strcmp(argv[iArg],"-quiet") == 0) {
      verbose = 0;
      showProgress = 0;
    } else if (strcmp(argv[iArg],"-nopr") == 0) {
      showProgress = 0;
    } else if (strcmp(argv[iArg],"-slow") == 0) {
      slow = 1;
    } else if (strcmp(argv[iArg],"-fastest") == 0) {
      fastest = 1;
      tophitsRefresh = 0.5;
      useTopHits2nd = true;
    } else if (strcmp(argv[iArg],"-2nd") == 0) {
      useTopHits2nd = true;
    } else if (strcmp(argv[iArg],"-no2nd") == 0) {
      useTopHits2nd = false;
    } else if (strcmp(argv[iArg],"-slownni") == 0) {
      fastNNI = false;
    } else if (strcmp(argv[iArg], "-matrix") == 0 && iArg < argc-1) {
      iArg++;
      opt->matrixPrefix = argv[iArg];
    } else if (strcmp(argv[iArg], "-nomatrix") == 0) {
      useMatrix = false;
    } else if (strcmp(argv[iArg], "-n") == 0 && iArg < argc-1) {
      iArg++;
      opt->nAlign = atoi(argv[iArg]);
      if (opt->nAlign < 1) {
	fprintf(stderr, "-n argument for #input alignments must be > 0 not %s\n", argv[iArg]);
	exit(1);
      }
    } else if (strcmp(argv[iArg], "-quote") == 0) {
      opt->bQuote = true;
    } else if (strcmp(argv[iArg], "-nt") == 0) {
      nCodes = 4;
    } else if (strcmp(argv[iArg], "-intree") == 0 && iArg < argc-1) {
      iArg++;
      opt->intreeFile = argv[iArg];
    } else if (strcmp(argv[iArg], "-intree1") == 0 && iArg < argc-1) {
      iArg++;
      opt->intreeFile = argv[iArg];
      opt->intree1 = true;
    } else if (strcmp(argv[iArg], "-nj") == 0) {
      bionj = 0;
    } else if (strcmp(argv[iArg], "-bionj") == 0) {
      bionj = 1;
    } else if (strcmp(argv[iArg], "-boot") == 0 && iArg < argc-1) {
      iArg++;
      opt->nBootstrap = atoi(argv[iArg]);
    } else if (strcmp(argv[iArg], "-noboot") == 0 || strcmp(argv[iArg], "-nosupport") == 0) {
      opt->nBootstrap = 0;
    } else if (strcmp(argv[iArg], "-seed") == 0 && iArg < argc-1) {
      iArg++;
      long seed = atol(argv[iArg]);
      ran_start(seed);
    } else if (strcmp(argv[iArg],"-top") == 0) {
      if(tophitsMult < 0.01)
	tophitsMult = 1.0;
    } else if (strcmp(argv[iArg],"-notop") == 0) {
      tophitsMult = 0.0;
    } else if (strcmp(argv[iArg], "-topm") == 0 && iArg < argc-1) {
      iArg++;
      tophitsMult = atof(argv[iArg]);
    } else if (strcmp(argv[iArg], "-close") == 0 && iArg < argc-1) {
      iArg++;
      tophitsClose = atof(argv[iArg]);
      if (tophitsMult <= 0) {
	fprintf(stderr, "Cannot use -close unless -top is set above 0\n");
	exit(1);
      }
      if (tophitsClose <= 0 || tophitsClose >= 1) {
	fprintf(stderr, "-close argument must be between 0 and 1\n");
	exit(1);
      }
    } else if (strcmp(argv[iArg], "-refresh") == 0 && iArg < argc-1) {
      iArg++;
      tophitsRefresh = atof(argv[iArg]);
      if (tophitsMult <= 0) {
	fprintf(stderr, "Cannot use -refresh unless -top is set above 0\n");
	exit(1);
      }
      if (tophitsRefresh <= 0 || tophitsRefresh >= 1) {
	fprintf(stderr, "-refresh argument must be between 0 and 1\n");
	exit(1);
      }
    } else if (strcmp(argv[iArg],"-nni") == 0 && iArg < argc-1) {
      iArg++;
      opt->nni = atoi(argv[iArg]);
      if (opt->nni == 0)
	opt->spr = 0;
    } else if (strcmp(argv[iArg],"-spr") == 0 && iArg < argc-1) {
      iArg++;
      opt->spr = atoi(argv[iArg]);
    } else if (strcmp(argv[iArg],"-sprlength") == 0 && iArg < argc-1) {
      iArg++;
      opt->maxSPRLength = atoi(argv[iArg]);
    } else if (strcmp(argv[iArg],"-mlnni") == 0 && iArg < argc-1) {
      iArg++;
      opt->MLnni = atoi(argv[iArg]);
    } else if (strcmp(argv[iArg],"-noml") == 0) {
      opt->MLnni = 0;
    } else if (strcmp(argv[iArg],"-mllen") == 0) {
      opt->MLnni = 0;
      opt->MLlen = true;
    } else if (strcmp(argv[iArg],"-nome") == 0) {
      opt->spr = 0;
      opt->nni = 0;
    } else if (strcmp(argv[iArg],"-help") == 0) {
      fprintf(stderr,"FastTree %s %s%s:\n%s", FT_VERSION, SSE_STRING, OpenMPString(), usage);
      exit(0);
    } else if (strcmp(argv[iArg],"-expert") == 0) {
      fprintf(stderr, "Detailed usage for FastTree %s %s%s:\n%s",
	      FT_VERSION, SSE_STRING, OpenMPString(), expertUsage);
      exit(0);
    } else if (strcmp(argv[iArg],"-pseudo") == 0) {
      if (iArg < argc-1 && isdigit(argv[iArg+1][0])) {
	iArg++;
	pseudoWeight = atof(argv[iArg]);
	if (pseudoWeight < 0.0) {
	  fprintf(stderr,"Illegal argument to -pseudo: %s\n", argv[iArg]);
	  exit(1);
	}
      } else {
	pseudoWeight = 1.0;
      }
    } else if (strcmp(argv[iArg],"-constraints") == 0 && iArg < argc-1) {
      iArg++;
      opt->constraintsFile = argv[iArg];
    } else if (strcmp(argv[iArg],"-constraintWeight") == 0 && iArg < argc-1) {
      iArg++;
      constraintWeight = atof(argv[iArg]);
      if (constraintWeight <= 0.0) {
	fprintf(stderr, "Illegal argument to -constraintWeight (must be greater than zero): %s\n", argv[iArg]);
	exit(1);
      }
    } else if (strcmp(argv[iArg],"-mlacc") == 0 && iArg < argc-1) {
      iArg++;
      mlAccuracy = atoi(argv[iArg]);
      if (mlAccuracy < 1) {
	fprintf(stderr, "Illlegal -mlacc argument: %s\n", argv[iArg]);
	exit(1);
      }
    } else if (strcmp(argv[iArg],"-exactml") == 0 || strcmp(argv[iArg],"-mlexact") == 0) {
      fprintf(stderr,"-exactml is not required -- exact posteriors is the default now\n");
    } else if (strcmp(argv[iArg],"-approxml") == 0 || strcmp(argv[iArg],"-mlapprox") == 0) {
      exactML = false;
    } else if (strcmp(argv[iArg],"-cat") == 0 && iArg < argc-1) {
      iArg++;
      opt->nRateCats = atoi(argv[iArg]);
      if (opt->nRateCats < 1) {
	fprintf(stderr, "Illlegal argument to -ncat (must be greater than zero): %s\n", argv[iArg]);
	exit(1);
      }
    } else if (strcmp(argv[iArg],"-nocat") == 0) {
      opt->nRateCats = 1;
    } else if (strcmp(argv[iArg], "-lg") == 0) {
        opt->bUseLg = true;
    } else if (strcmp(argv[iArg], "-wag") == 0) {
        opt->bUseWag = true;
    } else if (strcmp(argv[iArg], "-gtr") == 0) {
      opt->bUseGtr = true;
    } else if (strcmp(argv[iArg], "-trans") == 0 && iArg < argc-1) {
      iArg++;
      opt->transitionFile = argv[iArg];
    } else if (strcmp(argv[iArg], "-gtrrates") == 0 && iArg < argc-6) {
      opt->bUseGtr = true;
      opt->bUseGtrRates = true;
      int i;
      for (i = 0; i < 6; i++) {
	opt->gtrrates[i] = atof(argv[++iArg]);
	if (opt->gtrrates[i] < 1e-5) {
	  fprintf(stderr, "Illegal or too small value of GTR rate: %s\n", argv[iArg]);
	  exit(1);
	}
      }
    } else if (strcmp(argv[iArg],"-gtrfreq") == 0 && iArg < argc-4) {
      opt->bUseGtr = true;
      opt->bUseGtrFreq = true;
      int i;
      double sum = 0;
      for (i = 0; i < 4; i++) {
	opt->gtrfreq[i] = atof(argv[++iArg]);
	sum += opt->gtrfreq[i];
	if (opt->gtrfreq[i] < 1e-5) {
	  fprintf(stderr, "Illegal or too small value of GTR frequency: %s\n", argv[iArg]);
	  exit(1);
	}
      }
      if (fabs(1.0-sum) > 0.01) {
	fprintf(stderr, "-gtrfreq values do not sum to 1\n");
	exit(1);
      }
      for (i = 0; i < 4; i++)
	opt->gtrfreq[i] /= sum;
    } else if (strcmp(argv[iArg],"-log") == 0 && iArg < argc-1) {
      iArg++;
      opt->logfile = argv[iArg];
    } else if (strcmp(argv[iArg],"-gamma") == 0) {
      gammaLogLk = true;
    } else if (strcmp(argv[iArg],"-out") == 0 && iArg < argc-1) {
      iArg++;
      opt->fpOut = fopen(argv[iArg],"w");
      if(opt->fpOut==NULL) {
	fprintf(stderr,"Cannot write to %s\n",argv[iArg]);
	exit(1);
      }
    } else if (argv[iArg][0] == '-') {
      fprintf(stderr, "Unknown or incorrect use of option %s\n%s", argv[iArg], usage);
      exit(1);
    } else
      break;
  }
  if(iArg < argc-1) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }

  opt->fileName = iArg == (argc-1) ?  argv[argc-1] : NULL;
}
