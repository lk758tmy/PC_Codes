#include <stdio.h>

typedef struct {
    double d;
    char   c;
} CS;

typedef struct {
    char   c1;
    double d;
    char   c2;
} CS1;

typedef struct {
	char c1,c2;
	double d;
} CS2;

int main()
{
  CS  a;
  CS1 b;
	CS2 e;
  char c;

  printf("sizeof(char)=%ld\n", sizeof(c));
  printf("sizeof(CS)=%ld\n", sizeof(CS));
  printf("offset(a.d)=%ld, offset(a.c)=%ld\n",
	(char *)&a.d - (char *)&a, (char *)&a.c - (char *)&a);

  printf( "sizeof(CS1)=%ld\n", sizeof(CS1));
  printf( "offset(b.c1)=%ld, offset(b.d)=%ld, offset(b.c2)=%ld\n",
	(char *)&b.c1 - (char *)&b, (char *)&b.d - (char *)&b,
	(char *)&b.c2 - (char *)&b);
	
	printf( "sizeof(CS2)=%ld\n", sizeof(CS2));
	printf( "offset(e.c1)=%ld, offset(e.c2)=%ld, offset(e.d)=%ld\n",
		(char *)&e.c1 - (char *)&e, (char *)&e.c2 - (char *)&e,
		(char *)&e.d - (char *)&e);

  return 0;
}
//由于数据对齐等机制(与编译器有关)：1.CS1的占地比CS2多！
//2.使用MPI数据类型时，位移不能想当然！["地址修正量"]
