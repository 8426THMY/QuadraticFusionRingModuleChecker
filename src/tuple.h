#ifndef TUPLE_H
#define TUPLE_H


#include <stddef.h>


// A tuple is a finite sequence of indices
//     (t_1, t_2, ..., t_n),
// where t_i >= t_j for all i >= j.
typedef struct tuple {
	size_t *values;
	size_t min_value;
	size_t max_value;
	size_t length;
} tuple;


int tuple_init(
	tuple *const tup,
	const size_t min_value,
	const size_t max_value,
	const size_t length
);
void tuple_start(
	tuple *const tup,
	const size_t min_value,
	const size_t max_value,
	const size_t length
);
int tuple_next(tuple *const tup);
void tuple_print(const tuple *const tup);
void tuple_free(tuple *const tup);


#endif