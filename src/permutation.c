#include "permutation.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>


static void swap_values(
	integer_t *const values, const size_t i, const size_t j
);


// Initialize a permutation from an input array which is
// sorted in ascending order. We pad the permutation to
// the given length by inserting zeroes at the beginning.
int permutation_init(
	permutation *const perm, const size_t length,
	const integer_t *const input, const size_t input_length
){

	perm->length = length;
	if(length == 0){
		perm->values = NULL;
	}else{
		const size_t padding =
			(length > input_length) ? length - input_length : 0;

		perm->values = malloc(length*sizeof(*perm->values));
		if(perm->values == NULL){
			puts("Failed to allocate memory for permutation!");
			return 0;
		}
		memset(perm->values, 0, padding*sizeof(*perm->values));
		// By using the padding here, we truncate the
		// input if it's longer than the desired length.
		memcpy(
			&perm->values[padding], input, (length - padding)*sizeof(*input)
		);
	}

	return 1;
}

// Start a permutation, assuming its memory has already been allocated.
void permutation_start(
	permutation *const perm, const size_t length,
	const integer_t *const input, const size_t input_length
){

	perm->length = length;
	if(length == 0){
		perm->values = NULL;
	}else{
		const size_t padding =
			(length > input_length) ? length - input_length : 0;

		memset(perm->values, 0, padding*sizeof(*perm->values));
		// By using the padding here, we truncate the
		// input if it's longer than the desired length.
		memcpy(
			&perm->values[padding], input, (length - padding)*sizeof(*input)
		);
	}
}

int permutation_next(permutation *const perm){
	size_t i, j;

	// We change the indexing around from the standard
	// algorithm to avoid looping around to the largest
	// positive integer when i < 2.
	i = perm->length;
	while(i >= 2 && perm->values[i - 2] >= perm->values[i - 1]){
		--i;
	}
	if(i <= 1){
		return 0;
	}
	// At this point, i >= 2, so we can fix the indexing.
	i -= 2;

	j = perm->length - 1;
	while(j > i && perm->values[j] <= perm->values[i]){
		--j;
	}
	swap_values(perm->values, i, j);

	++i;
	j = perm->length - 1;
	while(i < j){
		swap_values(perm->values, i, j);
		++i;
		--j;
	}

	return 1;
}

// Reset a permutation to the first one in ascending order.
// This assumes that we're currently on the last permutation.
void permutation_restart(permutation *const perm){
	if(perm->length > 0){
		const size_t last = perm->length - 1;
		const size_t midpoint = last/2;

		size_t i;
		// If we're on the last permutation, we can return
		// to the first by simply reversing the values array.
		for(i = 0; i < midpoint; ++i){
			const integer_t temp = perm->values[i];
			perm->values[i] = perm->values[last-i];
			perm->values[last-i] = temp;
		}
	}
}

void permutation_print(const permutation *const perm){
	if(perm->length > 0){
		size_t i;

		printf("(%u", perm->values[0]);
		for(i = 1; i < perm->length; ++i){
			printf(", %u", perm->values[i]);
		}
		printf(")\n");
	}
}

void permutation_free(permutation *const perm){
	free(perm->values);
}


// Iteratively visit every possible set of permutations in the array.
int permutation_array_next(permutation_array *const perm_array){
	size_t i;
	for(i = 0; i < perm_array->num_perms - 1; ++i){
		if(permutation_next(&perm_array->perms[i])){
			return 1;
		}else{
			permutation_restart(&perm_array->perms[i]);
		}
	}
	return permutation_next(&perm_array->perms[perm_array->num_perms - 1]);
}


static void swap_values(
	integer_t *const values, const size_t i, const size_t j
){

	const integer_t temp = values[i];
	values[i] = values[j];
	values[j] = temp;
}