#ifndef __array_h__
#define __array_h__

#include <stdio.h>
#include <stdlib.h>

typedef struct Node Node;
struct Node {
  void* value;
  Node* next;
};

typedef struct {
  Node* head;
  Node* tail;
  int size;
} List;

static inline List*
createList()
{
  List* list = (List*)malloc(sizeof(List));
  list->head = NULL;
  list->tail = NULL;
  list->size = 0;
  return list;
}

static inline void
pushHeadList(List* list, void* value)
{
  Node* node = (Node*)malloc(sizeof(Node));
  node->value = value;
  node->next = list->head;
  list->head = node;
  if (list->tail == NULL) {
    list->tail = node;
  }
  list->size++;
}

static inline void
pushTailList(List* list, void* value)
{
  Node* node = (Node*)malloc(sizeof(Node));
  node->value = value;
  node->next = NULL;
  if (list->tail == NULL) {
    list->head = node;
    list->tail = node;
  } else {
    list->tail->next = node;
    list->tail = node;
  }
  list->size++;
}

static inline void*
popHeadList(List* list)
{
  if (list->head == NULL) {
    return NULL;
  }
  Node* node = list->head;
  void* value = node->value;
  list->head = node->next;
  if (list->tail == node) {
    list->tail = NULL;
  }
  free(node);
  list->size--;
  return value;
}

static inline void*
popTailList(List* list)
{
  if (list->tail == NULL) {
    return NULL;
  }
  Node* node = list->tail;
  void* value = node->value;
  if (list->head == node) {
    list->head = NULL;
    list->tail = NULL;
  } else {
    Node* prev = list->head;
    while (prev->next != node) {
      prev = prev->next;
    }
    prev->next = NULL;
    list->tail = prev;
  }
  free(node);
  list->size--;
  return value;
}


static inline void
freeList(List* list)
{
  Node* node = list->head;
  while (node != NULL) {
    Node* next = node->next;
    free(node->value);
    free(node);
    node = next;
  }
  free(list);
}

#define for_each_list(list, node)                                             \
  for (Node* node = list->head; node != NULL; node = node->next)

#endif // __array_h__
