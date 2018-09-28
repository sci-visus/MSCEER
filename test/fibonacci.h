/* file fibonacci.h */
#ifndef FIBONACCI_H
#define FIBONACCI_H

#include "testlibA.h"

typedef int int64;

class testclass {
public:
	int stuff;
	float fstuff;
};

testprint* make_testprint();
void testprint_printme(testprint* p);
int64 fibonacci_recursion( int64 n);
int64 fibonacci_iteration( int64 n, int64 nmax);
int64 fibonacci_full_iteration( int64 n, int64 nmax);
testclass* make_testclass();
int make1();
int make2();
int make3();
int make4();
#endif