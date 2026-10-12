#include "partition.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int partition_init(
	partition *const p, const size_t n, const size_t max_length
){

	p->max_length = max_length;
	if(max_length == 0){
		p->parts = NULL;
		p->length = 0;
	}else{
		p->parts = calloc(max_length, sizeof(*p->parts));
		if(p->parts == NULL){
			puts("Failed to allocate memory for partition!");
			return 0;
		}

		if(max_length == 1){
			p->parts[0] = n;
			p->length = 1;
		}else{
			p->parts[1] = n;
			p->length = 2;
			partition_next(p);
		}
	}

	return 1;
}

// Start a partition, assuming its memory has already been allocated.
void partition_start(
	partition *const p, const size_t n, const size_t max_length
){

	p->max_length = max_length;
	if(max_length == 0){
		p->length = 0;
	}else{
		if(max_length == 1){
			p->parts[0] = n;
			p->length = 1;
		}else{
			memset(p->parts, 0, max_length*sizeof(*p->parts));

			p->parts[1] = n;
			p->length = 2;
			partition_next(p);
		}
	}
}

// We use a modified version of the RuleAsc algorithm from Kelleher and
// O'Sullivan's Generating All Partitions: A Comparison Of Two Encodings.
// The function returns PARTITION_LAST when we've seen every partition.
int partition_next(partition *const p){
	if(p->length > 1){
		integer_t y = p->parts[--p->length] - 1;
		const integer_t x = p->parts[p->length - 1] + 1;

		while(x <= y && p->length < p->max_length){
			p->parts[p->length - 1] = x;
			y -= x;
			++p->length;
		}

		p->parts[p->length - 1] = x + y;
		return 1;
	}

	return 0;
}

// Reset a partition to the first one in ascending order.
void partition_restart(partition *const p){
	if(p->max_length > 1){
		integer_t n = 0;
		size_t i;
		for(i = 0; i < p->length; ++i){
			n += p->parts[i];
		}

		memset(p->parts, 0, p->max_length*sizeof(*p->parts));

		p->parts[1] = n;
		p->length = 2;
		partition_next(p);
	}
}

void partition_print(const partition *const p){
	if(p->length > 0){
		size_t i;

		printf("(%u", p->parts[0]);
		for(i = 1; i < p->length; ++i){
			printf(", %u", p->parts[i]);
		}
		printf(")\n");
	}
}

void partition_free(partition *const p){
	free(p->parts);
}


// Iteratively visit every possible set of partitions in the array.
int partition_array_next(partition_array *const part_array){
	size_t i;
	for(i = 0; i < part_array->num_parts - 1; ++i){
		if(partition_next(&part_array->parts[i])){
			return 1;
		}else{
			partition_restart(&part_array->parts[i]);
		}
	}
	return partition_next(&part_array->parts[part_array->num_parts - 1]);
}