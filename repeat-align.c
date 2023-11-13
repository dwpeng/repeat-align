#include "blosum50.h"
#include "list.h"
#include <assert.h>
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

// static const char* backtrace_str[] = {
//   "UP",
//   "LEFT",
//   "DIAG",
//   "JUMP",
// };

typedef enum {
  MATCH_MISMATCH = 0,
  MATCH_MATCH,
  MATCH_GAP_UP,
  MATCH_GAP_LEFT,
  MATCH_JUMP,
} MATCH_ENUM;

#define max2(a, b) ((a) > (b) ? (a) : (b))
#define max3(a, b, c) max2(max2(a, b), c)

typedef struct {
  char* seq;
  int len;
} seq_t;

typedef struct {
  int ref;
  int query;
  MATCH_ENUM match;
} match_t;

static inline match_t*
createMatch(int ref, int query, MATCH_ENUM match)
{
  match_t* m = (match_t*)malloc(sizeof(match_t));
  m->ref = ref;
  m->query = query;
  m->match = match;
  return m;
}

typedef struct {
  seq_t* ref;
  seq_t* query;
  List* alignment;
} align_result_t;

typedef enum {
  SEQ_PROTEIN = 0,
  SEQ_DNA,
} SEQ_ENUM;

typedef struct {
  int match;
  int mismatch;
  int gap;
  int T;
  SEQ_ENUM type;
} align_option_t;

#define _init_matrix(name, type, reflen, querylen)                            \
  type** name = (type**)malloc(sizeof(type*) * (querylen + 1));               \
  for (int i = 0; i <= querylen; i++) {                                       \
    name[i] = (type*)malloc(sizeof(type) * (reflen + 1));                     \
    memset(name[i], 0, sizeof(type) * (reflen + 1));                          \
  }

#define free_matrix(m, len)                                                   \
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
print_matrix(seq_t* ref, seq_t* query, cell_t** m, int flag)
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

      if (m[i][j].break_point == flag) {
        if (m[i][j].backtrace == BACKTRACE_JUMP) {
          printf("\033[32m%3d (%c)\033[0m", m[i][j].score, c);
        } else if (m[i][j].backtrace == BACKTRACE_DIAG) {
          printf("\033[31m%3d (%c)\033[0m", m[i][j].score, c);
        } else {
          printf("\033[33m%3d (%c)\033[0m", m[i][j].score, c);
        }
      } else {
        printf("%3d (%c)", m[i][j].score, c);
      }
    }
    printf("\n");
  }
}

static inline void
flash_print(seq_t* ref, seq_t* query, cell_t** m, int flag)
{
  printf("\033[2J\033[1;1H");
  print_matrix(ref, query, m, flag);
  usleep(150000);
}

static inline int
get_previ(cell_t** matrix, seq_t* query, int T, int col)
{
  int max_score = -9999999;
  int max_i = 0;
  for (int i = 0; i <= query->len; i++) {
    if (matrix[i][col - 1].score > max_score) {
      max_score = matrix[i][col - 1].score;
      max_i = i;
    }
  }

  if (max_i < query->len && max_score > T) {
    if (matrix[max_i + 1][col].backtrace == BACKTRACE_DIAG) {
      if (matrix[max_i + 1][col].score > matrix[max_i][col - 1].score) {
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
  result->alignment = createList();
  int previ = get_previ(matrix, query, T, ref->len);
  int prevj = ref->len;
  cell_t *cell = NULL, *cell_next = NULL;
  match_t* match = NULL;
  for (int i = ref->len; i >= 0; i--) {
    cell = &matrix[previ][prevj];
    cell->break_point = 2;
    flash_print(result->ref, result->query, matrix, 2);
    switch (cell->backtrace) {
    case BACKTRACE_UP:
      previ = previ - 1;
      prevj = prevj;
      match = createMatch(prevj, previ, MATCH_GAP_UP);
      pushHeadList(result->alignment, match);
      break;
    case BACKTRACE_LEFT:
      previ = previ;
      prevj = prevj - 1;
      match = createMatch(prevj, previ, MATCH_GAP_LEFT);
      pushHeadList(result->alignment, match);
      break;
    case BACKTRACE_DIAG:
      previ = previ - 1;
      prevj = prevj - 1;
      match = createMatch(prevj, previ, MATCH_MATCH);
      if (ref->seq[prevj] != query->seq[previ]) {
        match->match = MATCH_MISMATCH;
      }
      pushHeadList(result->alignment, match);
      break;
    case BACKTRACE_JUMP:
      previ = cell->prev_i;
      prevj = cell->prev_j;
      if (previ + 1 < query->len + 1 && prevj + 1 < ref->len + 1) {
        cell_next = &matrix[previ + 1][prevj + 1];
      }
      if (cell_next) {
        cell_next->break_point = 2;
        if (ref->seq[prevj] != query->seq[previ]) {
          match = createMatch(prevj, previ, MATCH_MISMATCH);
        } else {
          match = createMatch(prevj, previ, MATCH_MATCH);
        }
        cell_next = NULL;
      } else {
        match = createMatch(prevj, previ, MATCH_JUMP);
      }
      pushHeadList(result->alignment, match);
      break;
    default:
      fprintf(stderr, "Unknown backtrace: %d\n", cell->backtrace);
      exit(EXIT_FAILURE);
    }
  }
  match_t* m = popHeadList(result->alignment);
  if (m->match != MATCH_GAP_UP) {
    pushHeadList(result->alignment, m);
  }
  free_matrix(matrix, query->len);
  return result;
}

static inline int
make_score(align_option_t* option, char b1, char b2)
{
  // upper case
  b1 = b1 & 0xDF;
  b2 = b2 & 0xDF;
  int i1 = b1 - 'A';
  int i2 = b2 - 'A';
  if (option->type == SEQ_PROTEIN) {
    return blosum50[i1][i2];
  } else {
    return (b1 == b2) ? option->match : option->mismatch;
  }
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
      int score = make_score(option, ref->seq[j - 1], query->seq[i - 1]);
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
    flash_print(ref, query, matrix, 1);
  }
  return backtrace(ref, query, matrix, option->T);
}

void
print_align(align_result_t* result)
{
  List* alignment = result->alignment;
  char* ref_string = malloc(sizeof(char) * (alignment->size + 1));
  char* match_string = malloc(sizeof(char) * (alignment->size + 1));
  char* query_string = malloc(sizeof(char) * (alignment->size + 1));
  memset(ref_string, 0, sizeof(char) * (alignment->size + 1));
  memset(query_string, 0, sizeof(char) * (alignment->size + 1));
  int i = 0;
  for_each_list(alignment, node)
  {
    match_t* match = (match_t*)node->value;
    switch (match->match) {
    case MATCH_MATCH:
      ref_string[i] = result->ref->seq[match->ref];
      query_string[i] = result->query->seq[match->query];
      match_string[i] = '|';
      i++;
      break;
    case MATCH_MISMATCH:
      ref_string[i] = result->ref->seq[match->ref];
      query_string[i] = result->query->seq[match->query];
      match_string[i] = '*';
      i++;
      break;
    case MATCH_GAP_UP:
      ref_string[i] = '-';
      query_string[i] = result->query->seq[match->query];
      match_string[i] = '-';
      i++;
      break;
    case MATCH_GAP_LEFT:
      ref_string[i] = result->ref->seq[match->ref];
      query_string[i] = '-';
      match_string[i] = '-';
      i++;
      break;
    case MATCH_JUMP:
      ref_string[i] = result->ref->seq[match->ref];
      query_string[i] = '.';
      match_string[i] = '.';
      i++;
      break;
    }
  }
  ref_string[i] = '\0';
  query_string[i] = '\0';
  match_string[i] = '\0';
  printf("%s\n", ref_string);
  printf("%s\n", match_string);
  printf("%s\n", query_string);
  free(ref_string);
  free(query_string);
  free(match_string);
}

static inline SEQ_ENUM
detect_seq_type(char* seq)
{
  int len = strlen(seq);
  SEQ_ENUM t = SEQ_DNA;
  for (int i = 0; i < len; i++) {
    char c = seq[i];
    c &= 0xDF;
    if (c != 'A' || c != 'C' || c != 'G' || c != 'T') {
      t = SEQ_PROTEIN;
      break;
    }
  }
  return t;
}

static inline void
upper_string(char* s)
{
  int len = strlen(s);
  for (int i = 0; i < len; i++) {
    s[i] &= 0xDF;
  }
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
  upper_string(s1);
  upper_string(s2);
  int len1 = strlen(s1);
  int len2 = strlen(s2);

  printf("ref  :%3d %s\n", len1, s1);
  printf("query:%3d %s\n", len2, s2);

  seq_t ref = {
    .seq = s1,
    .len = len1,
  };

  seq_t query = {
    .seq = s2,
    .len = len2,
  };

  SEQ_ENUM t1 = detect_seq_type(s1);
  SEQ_ENUM t2 = detect_seq_type(s2);
  assert(t1 == t2);
  align_option_t option = {
    .gap = -4,
    .match = 4,
    .mismatch = -4,
    .T = 20,
    .type = t1,
  };

  align_result_t* result = repeat_align(&ref, &query, &option);
  print_align(result);
  if (result != NULL) {
    freeList(result->alignment);
  }
}
