/* file fibonacci.cpp */
 #include <iostream>
 #include "stdio.h"
 #include "stdlib.h"
#include "fibonacci.h"
 #include <cmath>

testclass* make_testclass() {
	return new testclass();
}

testprint* make_testprint() {
	return new testprint();
}
void testprint_printme(testprint* p) {
	p->printme();
}

int make1() { return 1; }
int make2() { return 2; }
int make3() { return 3; }
int make4() { return barf(3); }
int myval = 3;
int64 fibonacci_recursion( int64 n)
 {
 int64 a=1;
 if(n >= 2)
 {
 a = fibonacci_recursion(n-1) + fibonacci_recursion(n-2);
 }
 // prlongf( ", %d", a );
 return a;
 }

int64 fibonacci_iteration(int64 n, int64 nmax)
 {
 int64 a, b,i;
 a = 0;
 b = 1;
 for(int k=0; k<n; k++){
    i = (a + b) % nmax;
    a = b;
    b = i;
 }
 return a;
 }

int64 fibonacci_full_iteration( int64 n, int64 nmax)
 {
 int64 n0 = 0;
 for(int64 i=0; i<n; i++)
 {
    n0+= fibonacci_iteration(i, nmax);
 }
 return n0;
 }

int main( int argc, char** argv) {
 if ( argc > 3 || argc == 1 ) {
 printf("argc: %d\n", argc);
 printf("Fibonacci takes one postitive integer greater\
 than two as it's argument\n");
 return EXIT_SUCCESS;
 }

int64 n = atof( argv[1] );
 int64 nmax_exp = atof( argv[2]);
 int64 nmax = pow(10, nmax_exp);
 int64 ret=fibonacci_full_iteration(n, nmax);
 std::cout << ret << std::endl;
 return EXIT_SUCCESS;
 }