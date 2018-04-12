#define STEP 3072

#include "secp.c"
#define iter_64(a) do {int _i; for (_i = 0; _i < 64; _i++) { a(_i) }} while (0)



__kernel void hello_kernel() {
  // Get the index of the current element to be processed
    int i = get_global_id(0), j;
  
	float f1=0,f2=1,f3=1;
	for (j=0; j<1000; j++) {
		f1 += i*j/(i+j);
		f2 *= f1;
		f3 = pow(f2, j);
	}
}

__kernel void vector_add(__global const int *A, __global const int *B, __global int *C) {
 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
 
    // Do the operation
    C[i] = A[i] + B[i];

}