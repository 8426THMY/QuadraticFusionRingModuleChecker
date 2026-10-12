#ifndef MODULE_CHECKER_H
#define MODULE_CHECKER_H


#include <stddef.h>


#include "quadratic_ring.h"
#include "group.h"
#include "util_math.h"

#include "tuple.h"
#include "linear_diophantine.h"
#include "partition.h"
#include "permutation.h"

#include "pool_allocator.h"


// Store all of the data we use while searching
// for modules over a quadratic fusion ring.
typedef struct module_checker {
	const quadratic_ring *R;

	size_t min_orbits;
	size_t num_orbits;
	size_t max_orbits;

	size_t **psi_ids;
	group_elem_id *inner_factors;
	size_t num_inner_factors;

	tuple H_ids;
	integer_t ***a;

	linear_diophantine diag_eq;
	partition_array diag_sums;
	permutation_array diag_mults;

	linear_diophantine offdiag_eqs;
	partition_array offdiag_sums;
	permutation_array offdiag_mults;

	pool_allocator pool;
} module_checker;


int module_checker_init(
	module_checker *const mc,
	const quadratic_ring *const R
);
int module_checker_find_modules(const quadratic_ring *const R);
void module_checker_delete(module_checker *const mc);


#endif