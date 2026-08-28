#ifndef PARTITION_H
#define PARTITION_H


#include <stddef.h>

#include "util_math.h"


typedef struct partition {
	integer_t *parts;
	size_t length;
	size_t max_length;
} partition;


int partition_init(
	partition *const p, const size_t n, const size_t max_length
);
void partition_start(
	partition *const p, const size_t n, const size_t max_length
);
int partition_next(partition *const p);
void partition_free(partition *const p);


#endif