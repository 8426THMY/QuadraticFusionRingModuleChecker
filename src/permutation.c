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

void permutation_free(permutation *const perm){
	free(perm->values);
}


static void swap_values(
	integer_t *const values, const size_t i, const size_t j
){

	const integer_t temp = values[i];
	values[i] = values[j];
	values[j] = temp;
}