#include "blosum50.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef enum {
  BACKTRACE_UP = 0,
  BACKTRACE_LEFT,
  BACKTRACE_DIAG,
  BACKTRACE_JUMP,
} BACKTRACE_ENUM;

static char* backtrace_str[] = {
  "UP",
  "LEFT",
  "DIAG",
  "JUMP",
};

typedef enum {
  MATCH_MISMATCH = 0,
  MATCH_MATCH,
  MATCH_GAP,
  MATCH_JUMP,
} MATCH_ENUM;

#define max2(a, b) ((a) > (b) ? (a) : (b))
#define max3(a, b, c) max2(max2(a, b), c)

typedef struct {
  char* seq;
  int len;
} seq_t;

typedef struct {
  seq_t* ref;
  seq_t* query;
  int* query_alignment;
  MATCH_ENUM* alignment;
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
  int score;
  int break_point;
  BACKTRACE_ENUM backtrace;
  int prev_i;
  int prev_j;
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
      printf("     ");
      for (int j = 0; j < ref->len; j++) {
        printf("      %c", ref->seq[j]);
      }
      printf("\n");
    }
    for (int j = 0; j <= ref->len; j++) {
      if (j == 0) {
        if (i == 0) {
          printf("  ");
        } else {
          printf(" %c", query->seq[i - 1]);
        }
      }
      char c;
      switch (m[i][j].backtrace) {
      case BACKTRACE_UP:
        c = '^';
        break;
      case BACKTRACE_LEFT:
        c = '<';
        break;
      case BACKTRACE_DIAG:
        c = '\\';
        break;
      case BACKTRACE_JUMP:
        c = '*';
        break;
      }
      if (m[i][j].break_point) {
        // green
        printf("\033[31m%3d (%c)\033[0m", m[i][j].score, c);
      } else {
        printf("%3d (%c)", m[i][j].score, c);
      }
    }
    printf("\n");
  }
}

static inline int
get_previ(cell_t** matrix, seq_t* ref, seq_t* query, int T)
{
  // check ref.len - 1 column
  int max_score = -9999999;
  int max_i = 0;
  for (int i = 0; i <= query->len; i++) {
    if (matrix[i][ref->len - 1].score > max_score) {
      max_score = matrix[i][ref->len - 1].score;
      max_i = i;
    }
  }

  if (max_i < query->len && max_score > T) {
    if (matrix[max_i + 1][ref->len].backtrace == BACKTRACE_DIAG) {
      if (matrix[max_i + 1][ref->len].score
          > matrix[max_i][ref->len - 1].score) {
        return max_i + 1;
      }
    }
  }

  if (max_i == query->len && max_score > T) {
    return 0;
  }

  return 0;
}

static inline align_result_t*
backtrace(seq_t* ref, seq_t* query, cell_t** matrix, int T)
{
  align_result_t* result = (align_result_t*)malloc(sizeof(align_result_t));
  result->ref = ref;
  result->query = query;
  result->query_alignment = (int*)malloc(sizeof(int) * (ref->len + 1));
  result->alignment = (MATCH_ENUM*)malloc(sizeof(MATCH_ENUM) * (ref->len + 1));
  int previ = get_previ(matrix, ref, query, T);
  printf("previ: %d\n", previ);
  int prevj = ref->len;
  cell_t* cell = NULL;
  for (int i = ref->len - 1; i >= 0; i--) {
    cell = &matrix[previ][prevj];
#if 0
    printf("previ: %d, prevj: %d, b: %s, %d, i: %d %c:%c\n", previ, prevj,
           backtrace_str[cell->backtrace], cell->score, i, ref->seq[prevj],
           query->seq[previ]);
#endif
    switch (cell->backtrace) {
    case BACKTRACE_UP:
      previ = previ - 1;
      prevj = prevj;
      result->alignment[i] = MATCH_GAP;
      result->query_alignment[i] = previ;
      break;
    case BACKTRACE_LEFT:
      previ = previ;
      prevj = prevj - 1;
      result->alignment[i] = MATCH_GAP;
      result->query_alignment[i] = previ;
      break;
    case BACKTRACE_DIAG:
      previ = previ - 1;
      prevj = prevj - 1;
      if (ref->seq[prevj] == query->seq[previ]) {
        result->alignment[i] = MATCH_MATCH;
      } else {
        result->alignment[i] = MATCH_MISMATCH;
      }
      result->query_alignment[i] = previ;
      break;
    case BACKTRACE_JUMP:
      previ = cell->prev_i;
      prevj = cell->prev_j;
      result->alignment[i] = MATCH_JUMP;
      result->query_alignment[i] = previ;
      break;
    default:
      fprintf(stderr, "Unknown backtrace: %d\n", cell->backtrace);
      exit(EXIT_FAILURE);
    }
  }
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
          matrix[i][j].break_point = 1;
          matrix[i][j].backtrace = BACKTRACE_JUMP;
          matrix[i][j].prev_i = mi;
          matrix[i][j].prev_j = mj;
        } else {
          matrix[i][j].score = matrix[i][j - 1].score;
          matrix[i][j].backtrace = BACKTRACE_LEFT;
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

      if (F_i_1_j_1 == max_F) {
        matrix[i][j].score = F_i_1_j_1;
        matrix[i][j].backtrace = BACKTRACE_DIAG;
        continue;
      }

      if (F_i_0 == max_F) {
        matrix[i][j].score = F_i_0;
        matrix[i][j].backtrace = matrix[0][j].backtrace;
        matrix[i][j].prev_i = matrix[0][j].prev_i;
        matrix[i][j].prev_j = matrix[0][j].prev_j;
        continue;
      }

      if (F_i_1_j == max_F) {
        matrix[i][j].score = F_i_1_j;
        matrix[i][j].backtrace = BACKTRACE_UP;
        continue;
      }
      if (F_i_j_1 == max_F) {
        matrix[i][j].score = F_i_j_1;
        matrix[i][j].backtrace = BACKTRACE_LEFT;
        continue;
      }
    }
    printf("\033[2J\033[1;1H");
    print_matrix(ref, query, matrix);
    usleep(150000);
  }
  return backtrace(ref, query, matrix, option->T);
}

void
print_align(align_result_t* result)
{
  int* query_alignment = result->query_alignment;
  MATCH_ENUM* alignment = result->alignment;
  char* buffer = (char*)malloc(sizeof(char) * (result->ref->len + 1));
  char* query = (char*)malloc(sizeof(char) * (result->query->len + 1));
  for (int i = 0; i < result->ref->len; i++) {
    switch (alignment[i]) {
    case MATCH_MISMATCH:
      buffer[i] = '*';
      query[i] = result->query->seq[query_alignment[i]];
      break;
    case MATCH_MATCH:
      buffer[i] = '|';
      query[i] = result->query->seq[query_alignment[i]];
      break;
    case MATCH_GAP:
      buffer[i] = '-';
      query[i] = '-';
      break;
    case MATCH_JUMP:
      buffer[i] = '.';
      query[i] = '.';
      break;
    }
  }
  buffer[result->ref->len] = '\0';
  query[result->ref->len] = '\0';
  printf("%s\n", result->ref->seq);
  printf("%s\n", buffer);
  printf("%s\n", query);
  free(buffer);
  free(query);
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
    free(result->query_alignment);
    free(result->alignment);
    free(result);
  }
}
