#include <stdio.h>
#include <assert.h>

#include "list.h"

printNode(nodeptr head, char *name)
{
     printf("%s->0x%x\n\tdata=%d\n\tprev->0x%x\n\tnext->0x%x\n", name, head, head->data, head->prev, head->next);
//     printf("%d ", head->data);
}

printList(nodeptr head)
{
     nodeptr nxt;
     nxt = head;
     printf("[");
     while(nxt)
     {
      printNode(nxt, "");
      nxt = nxt->next;
     }
     printf("]\n");
}

/* test list
////////////////
int main (void)
{
     nodeptr head, new, new2;
     int dat;
     head = createNode(NULL); // create the head node
     printf("Created head with null data\n");
     printNode(head, "head");

     dat = 10;
     printf("Created new with integer data\n");
     new = createNode((void*)dat);
     printNode(new, "new");

     addNode(head, new);
     printf("Added new to head\n");
     printList(head);

     removeNode(new);
     printf("Removed new from head\n");
     printNode(new, "new");
     printList(head);

     addNode(head, new);
     printf("Added new to head\n");
     printList(head);

     dat = 12;
     printf("Created new2 with integer data\n");
     new2 = createNode((void*)dat);
     printNode(new2, "new2");

     addNode(head, new2);
     printf("Added new2 to head\n");
     printList(head);

     deleteNode(new);
     printf("Deleted new->0x%x\n", new);
     printList(head);

     dat = 10;
     printf("Created new with integer data\n");
     new = createNode((void*)dat);
     printNode(new, "new");

     addNode(new2, new);
     printf("Added new to new2\n");
     printList(head);

     printf("Deleting list at head\n");
     deleteList(head);
}
*/
/////////////////////////////

nodeptr createNode(dataptr dat)
{
     nodeptr head;

     head = (nodeptr)malloc(sizeof(node));
     if (!head)
      return NULL;
     head->next = NULL;
     head->prev = NULL;
     head->data = dat;

     return head;
}

nodeptr addNode(nodeptr head, nodeptr new)
{
//     assert(head); assert(new);

     new->next = head->next;
     new->prev = head;
     if (head->next)
      head->next->prev = new;
     head->next = new;

     return new;     
}

dataptr deleteNode(nodeptr head)
{
//     assert(head);

     dataptr dat = head->data;
     head = removeNode(head);
     printf("Freeing node at 0x%x\n", head);
     free(head); // may be called for NULL head
     return dat;
}

nodeptr removeNode(nodeptr head)
{
//     assert(head);

     nodeptr p, n;

     p = head->prev; // NULL prev represents list head
     n = head->next; // could be NULL indicating tail

     if (p) p->next = n;
     if (n) n->prev = p;

     head->prev = NULL;
     head->next = NULL;

     return head;
}

void deleteList(nodeptr head)
{
//     assert(head);

     nodeptr prv, nxt;
     while (head->prev)
      head=head->prev; // start at beginning
     prv = head;
     nxt = head->next;
     while(nxt)
     {
      deleteNode(prv);
      prv = nxt;
      nxt = nxt->next;
     }
     deleteNode(prv);
}
