#include "blosum50.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef enum {
  ALIGN_MATCH = 0,
  ALIGN_MISMATCH = 1,
  ALIGN_GAP = 2,
} ALIGN_ENUM;

typedef enum {
  BACKTRACE_UP = 0,
  BACKTRACE_LEFT = 1,
  BACKTRACE_DIAG = 2,
} BACKTRACE_ENUM;

static inline int
max2(int a, int b)
{
  return a > b ? a : b;
}

static inline int
max3(int a, int b, int c)
{
  return max2(max2(a, b), c);
}

typedef struct {
  char* seq;
  int len;
} seq_t;

typedef struct {
  seq_t* ref;
  seq_t* query;
  char* alignment;
} align_result_t;

typedef struct {
  int match;
  int mismatch;
  int gap;
  int T;
} align_option_t;

#define _init_matrix(name, type, reflen, querylen)                            \
  type** name = (type**)malloc(sizeof(type*) * (querylen + 1));               \
  for (int i = 0; i <= querylen; i++) {                                       \
    name[i] = (type*)malloc(sizeof(type) * (reflen + 1));                     \
    memset(name[i], 0, sizeof(type) * (reflen + 1));                          \
  }

#define free_matrxi(m, len)                                                   \
  for (int i = 0; i <= len; i++) {                                            \
    free(m[i]);                                                               \
  }                                                                           \
  free(m);

typedef struct {
  int prev_i;
  int prev_j;
  int score;
} cell_t;

static inline cell_t**
init_score_matrix(int refLen, int queryLen)
{
  _init_matrix(matrix, cell_t, refLen, queryLen);
  return matrix;
}

void
print_matrix(seq_t* ref, seq_t* query, cell_t** m)
{
  printf("\n");
  for (int i = 0; i <= query->len; i++) {
    if (i == 0) {
      printf("      ");
      for (int j = 0; j < ref->len; j++) {
        printf(" %c ", ref->seq[j]);
      }
      printf("\n");
    }
    for (int j = 0; j <= ref->len; j++) {
      if (j == 0) {
        if (i == 0) {
          printf("   ");
        } else {
          printf(" %c ", query->seq[i - 1]);
        }
      }
      printf("%2d ", m[i][j].score);
    }
    printf("\n");
  }
}

static inline align_result_t*
backtrace(seq_t* ref, seq_t* query, cell_t** matrix)
{
  align_result_t* result = (align_result_t*)malloc(sizeof(align_result_t));
  result->ref = ref;
  result->query = query;
  result->alignment = (char*)malloc(sizeof(char) * max2(ref->len, query->len));
  // TODO
  // backtrace core algorithm

  // free matrix
  free_matrxi(matrix, query->len);
  return result;
}

static inline int
make_score(char b1, char b2)
{
  // upper case
  b1 = b1 & 0xDF;
  b2 = b2 & 0xDF;
  int i1 = b1 - 'A';
  int i2 = b2 - 'A';
  return blosum50[i1][i2];
}

align_result_t*
repeat_align(seq_t* ref, seq_t* query, align_option_t* option)
{
  // cols -> query -> i
  // rows -> ref   -> j
  //
  //   ------------------> j(ref)
  //  |
  //  |
  //  |
  //  |
  //  i(query)

  cell_t** matrix = init_score_matrix(ref->len, query->len);
  // TODO
  // fill matrix and path
  // align core algorithm

  matrix[0][0].score = 0;
  for (int j = 1; j <= ref->len; j++) {
    for (int i = 0; i <= query->len; i++) {

      if (i == 0) {
        int mi, mj, ms = -9999999;
        for (int k = 0; k <= query->len; k++) {
          if (matrix[k][j - 1].score - option->T >= ms) {
            ms = matrix[k][j - 1].score - option->T;
            mi = k;
            mj = j - 1;
          }
        }
        if (ms >= matrix[i][j - 1].score) {
          matrix[i][j].score = ms;
          matrix[i][j].prev_i = mi;
          matrix[i][j].prev_j = mj;
        } else {
          matrix[i][j].score = matrix[i][j - 1].score;
          matrix[i][j].prev_i = i;
          matrix[i][j].prev_j = j - 1;
        }
        continue;
      }

      int score = make_score(ref->seq[j - 1], query->seq[i - 1]);
      int F_i_1_j_1 = matrix[i - 1][j - 1].score + score;
      int F_i_j_1 = matrix[i][j - 1].score + option->gap;
      int F_i_1_j = matrix[i - 1][j].score + option->gap;
      int F_i_0 = matrix[0][j].score;

      int max_F = max3(F_i_1_j_1, F_i_1_j, F_i_j_1);
      max_F = max2(max_F, F_i_0);

      if (F_i_0 == max_F) {
        matrix[i][j].score = F_i_0;
        matrix[i][j].prev_i = i;
        matrix[i][j].prev_j = 0;
        continue;
      }

      if (F_i_1_j_1 == max_F) {
        matrix[i][j].score = F_i_1_j_1;
        matrix[i][j].prev_i = i - 1;
        matrix[i][j].prev_j = j - 1;
        continue;
      }
      if (F_i_1_j == max_F) {
        matrix[i][j].score = F_i_1_j;
        matrix[i][j].prev_i = i - 1;
        matrix[i][j].prev_j = j;
        continue;
      }
      if (F_i_j_1 == max_F) {
        matrix[i][j].score = F_i_j_1;
        matrix[i][j].prev_i = i;
        matrix[i][j].prev_j = j - 1;
        continue;
      }
    }
    system("clear");
    print_matrix(ref, query, matrix);
    sleep(1);
  }
  return backtrace(ref, query, matrix);
}

void
print_align(align_result_t* result)
{
  // TODO
  // print alignment
}

int
main(int argc, char* argv[])
{
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <string1> <string2>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  char* s1 = argv[1];
  char* s2 = argv[2];
  int len1 = strlen(s1);
  int len2 = strlen(s2);

  printf("ref  :%3d %s\n", len1, s1);
  printf("query:%3d %s\n", len2, s2);

  align_option_t option = {
    .gap = -8,
    .T = 20,
  };

  seq_t ref = {
    .seq = s1,
    .len = len1,
  };

  seq_t query = {
    .seq = s2,
    .len = len2,
  };

  align_result_t* result = repeat_align(&ref, &query, &option);
  print_align(result);
  if (result != NULL) {
    free(result->alignment);
    free(result);
  }
}
