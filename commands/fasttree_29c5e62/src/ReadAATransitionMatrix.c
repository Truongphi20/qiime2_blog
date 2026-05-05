#include "hyper_parameters.h"
#include "support_functions.h"
#include "datastructs.h"


transition_matrix_t *ReadAATransitionMatrix(/*IN*/char *filename) {
  assert(nCodes==20);
  double stat[20];
  static double matrix[MAXCODES][MAXCODES];
  static char buf[BUFFER_SIZE];
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot read transition matrix file %s\n", filename);
    exit(1);
  }
  char expected[2*MAXCODES+20];
  int posE = 0;
  int i, j;
  for (i = 0; i < 20; i++) {
    expected[posE++] = codesStringAA[i];
    expected[posE++] = '\t';
  }
  expected[posE++] = '*';
  expected[posE++] = '\n';
  expected[posE++] = '\0';
  
  if (fgets(buf, sizeof(buf), fp) == NULL) {
    fprintf(stderr, "Error reading header line from transition matrix file\n");
    exit(1);
  }
  if (strcmp(buf, expected) != 0) {
    fprintf(stderr, "Invalid header line in transition matrix file, it must match:\n%s\n", expected);
    exit(1);
  }
  for (i = 0; i < 20; i++) {
    if (fgets(buf, sizeof(buf), fp) == NULL) {
      fprintf(stderr, "Error reading matrix line\n");
      exit(1);
    }
    char *field = strtok(buf,"\t\r\n");
    if (field == NULL || strlen(field) != 1 || field[0] != codesStringAA[i]) {
      fprintf(stderr, "Line for amino acid %c does not have the expected beginning\n", codesStringAA[i]);
      exit(1);
    }
    for (j = 0; j < 20; j++) {
      field = strtok(NULL, "\t\r\n");
      if (field == NULL) {
        fprintf(stderr, "Not enough fields for amino acid %c\n", codesStringAA[i]);
        exit(1);
      }
      matrix[i][j] = atof(field);
    }
    field = strtok(NULL, "\t\r\n");
    if (field == NULL) {
      fprintf(stderr, "Not enough fields for amino acid %c\n", codesStringAA[i]);
      exit(1);
    }
    stat[i] = atof(field);
  }

  double tol = 1e-5;
  /* Verify that stat is positive and sums to 1 */
  double statTot = 0;
  for (i = 0; i < 20; i++) {
    if (stat[i] < tol) {
      fprintf(stderr, "stationary frequency for amino acid %c must be positive\n", codesStringAA[i]);
      exit(1);
    }
    statTot += stat[i];
  }
  if (fabs(statTot - 1) > tol) {
    fprintf(stderr, "stationary frequencies must sum to 1 -- actual sum is %g\n", statTot);
    exit(1);
  }

  /* Verify that diagonals are negative and dot product of stat and diagonals is -1 */
  double totRate = 0;
  for (i = 0; i < 20; i++) {
    double diag = matrix[i][i];
    if (diag > -tol) {
      fprintf(stderr, "transition rate(%c,%c) must be negative\n",
              codesStringAA[i], codesStringAA[i]);
      exit(1);
    }
    totRate += stat[i] * diag;
  }
  if (fabs(totRate + 1) > tol) {
    fprintf(stderr, "Dot product of matrix diagonal and stationary frequencies must be -1 -- actual dot product is %g\n",
            totRate);
    exit(1);
  }

  /* Verify that each off-diagonal entry is nonnegative and that each column sums to 0 */
  for (j = 0; j < 20; j++) {
    double colSum = 0;
    for (i = 0; i < 20; i++) {
      double value = matrix[i][j];
      colSum += value;
      if (i != j && value < 0) {
        fprintf(stderr, "Off-diagonal matrix entry for (%c,%c) is negative\n",
                codesStringAA[i], codesStringAA[j]);
        exit(1);
      }
    }
    if (fabs(colSum) > tol) {
      fprintf(stderr, "Sum of column %c must be zero -- actual sum is %g\n",
              codesStringAA[j], colSum);
      exit(1);
    }
  }
  return CreateTransitionMatrix(matrix, stat);
}