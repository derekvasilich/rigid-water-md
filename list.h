#ifndef LIST_H
#define LIST_H

struct node_t
{
     void *data;
     struct node_t *prev;
     struct node_t *next;
};

typedef void* dataptr;
typedef struct node_t node;
typedef node* nodeptr;

// create (allocate) a new node with NULL head and tail
nodeptr createNode(dataptr dat);

// add node new after head 
nodeptr addNode(nodeptr head, nodeptr new);

// remove node head from list it is contained in
nodeptr removeNode(nodeptr head);

// delete (deallocate) node head
dataptr deleteNode(nodeptr head);

// delete (deallocate) a list and all of its nodes (prev and next)
// it is assumed that deallocation of list data is performed
void deleteList(nodeptr head);

#endif // LIST_H
