#ifndef PERMUTATION_H
#define PERMUTATION_H


#include <stddef.h>

#include "util_math.h"


typedef struct permutation {
	integer_t *values;
	size_t length;
} permutation;

typedef struct permutation_array {
	permutation *perms;
	size_t num_perms;
} permutation_array;


int permutation_init(
	permutation *const perm, const size_t length,
	const integer_t *const input, const size_t input_length
);
void permutation_start(
	permutation *const perm, const size_t length,
	const integer_t *const input, const size_t input_length
);
int permutation_next(permutation *const perm);
void permutation_restart(permutation *const perm);
void permutation_print(const permutation *const perm);
void permutation_free(permutation *const perm);

int permutation_array_next(permutation_array *const perm_array);


#endif