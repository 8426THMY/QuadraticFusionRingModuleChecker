#include "tuple.h"


#include <stdlib.h>
#include <stdio.h>


int tuple_init(
	tuple *const tup,
	const size_t min_value,
	const size_t max_value,
	const size_t length
){

	tup->values = malloc(length*sizeof(*tup->values));
	if(tup->values == NULL){
		puts("Failed to allocate memory for tuple!");
		return 0;
	}
	if(min_value <= max_value){
		tup->min_value = min_value;
		tup->max_value = max_value;
	}else{
		tup->min_value = min_value;
		tup->max_value = max_value;
	}
	tup->length  = length;

	{
		size_t i;
		for(i = 0; i < length; ++i){
			tup->values[i] = min_value;
		}
	}

	return 1;
}

void tuple_start(
	tuple *const tup,
	const size_t min_value,
	const size_t max_value,
	const size_t length
){

	if(min_value <= max_value){
		tup->min_value = min_value;
		tup->max_value = max_value;
	}else{
		tup->min_value = min_value;
		tup->max_value = max_value;
	}
	tup->length  = length;

	{
		size_t i;
		for(i = 0; i < length; ++i){
			tup->values[i] = min_value;
		}
	}
}

// Iteratively visit every possible set of tuples.
int tuple_next(tuple *const tup){
	if(tup->length > 0){
		size_t i = tup->length - 1;
		for(;;){
			if(tup->values[i] < tup->max_value){
				size_t j;

				++tup->values[i];

				for(j = i + 1; j < tup->length; j++){
					tup->values[j] = tup->values[j - 1];
				}

				return 1;
			}

			// This basically just to avoid integer underflow.
			if(i == 0){
				break;
			}
			--i;
		}
	}

	return 0;
}

void tuple_print(const tuple *const tup){
	if(tup->length > 0){
		size_t i;

		printf("(%lu", (long unsigned int)tup->values[0]);
		for(i = 1; i < tup->length; ++i){
			printf(", %lu", (long unsigned int)tup->values[i]);
		}
		printf(")\n");
	}
}

void tuple_free(tuple *const tup){
	free(tup->values);
}