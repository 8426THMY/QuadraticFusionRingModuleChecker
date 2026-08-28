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
		memset(p->parts, 0, max_length*sizeof(*p->parts));

		if(max_length == 1){
			p->parts[0] = n;
			p->length = 1;
		}else{
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

void partition_free(partition *const p){
	free(p->parts);
}