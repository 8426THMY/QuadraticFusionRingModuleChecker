#ifndef LINEAR_DIOPHANTINE_H
#define LINEAR_DIOPHANTINE_H


#include <stddef.h>

#include "util_math.h"


// Represents the system of linear Diophantine equations
//     a[i][0]x[0] + a[i][1]x[1] + ... + a[i][n-1]x[n-1] = b[i],
// in the variables x[j]. Currently, we only solve the trivial
// cases where we have one equation in one or two variables.
typedef struct linear_diophantine {
	integer_t *x;
	size_t num_vars;
	integer_t **a;
	integer_t *b;
	size_t num_eqs;
	// Store gcd(a[0], a[1], ..., a[n-1]).
	integer_t g;
} linear_diophantine;


int linear_diophantine_start(linear_diophantine *const dio);
int linear_diophantine_next(linear_diophantine *const dio);


#endif