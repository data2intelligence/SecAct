#include <R.h>

void hello(int *n)
{
  int i;
  for(i=0; i < *n; i++)
  {
    Rprintf("Hello, world!\n");
  }
}

