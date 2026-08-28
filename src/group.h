#ifndef GROUP_H
#define GROUP_H


#include <stddef.h>
#include <stdio.h>

#include "util_math.h"


typedef integer_t *pgroup_elem;
typedef pgroup_elem *group_elem;
typedef size_t group_elem_id;

// A quotient group is stored as a bidirectional set
// of cosets and group elements. Given the index of
// one, we can use this structure to find the other.
typedef struct quotient_group {
	// For each coset, store the indices of
	// the group elements that live in it.
	group_elem_id **cosets;
	size_t num_cosets;
	// For each group element, store the
	// index of the coset that it lives in.
	size_t *coset_ids;
} quotient_group;

typedef struct subgroup {
	group_elem_id *elem_ids;
	size_t order;
	quotient_group quotient;
} subgroup;

typedef struct pgroup {
	integer_t *powers;
	size_t num_powers;
	integer_t prime;
	size_t order;
} pgroup;

typedef struct group {
	pgroup *pgroups;
	size_t num_pgroups;
	size_t num_factors;

	group_elem *elems;
	size_t order;

	// Group multiplication table. We store the
	// ID of the element gh in multiply[g][h].
	group_elem_id **multiply;
	// Group inverse table. We store the ID
	// of the element g^{-1} in inverse[g].
	group_elem_id *inverse;

	subgroup *subgroups;
	size_t num_subgroups;
} group;


group_elem_id group_elem_get_id(
	const group *const G,
	const group_elem g
);
group_elem_id group_elem_multiply(
	const group *const G,
	const group_elem g, const group_elem h
);
group_elem_id group_elem_invert(
	const group *const G,
	const group_elem g
);
void group_elem_print(const group *const G, const group_elem g);

int subgroup_init_quotient(const group *const G, subgroup *const H);
#if 0
void subgroup_init_zero(
	const group *const G, subgroup *const H,
	const unsigned int num_gens
);
void subgroup_init_direct_sum(
	const group *const G, subgroup *const H,
	const subgroup *const H1, const subgroup *const H2
);
#endif

int pgroup_load(pgroup *const Gp, FILE *const fp);
void pgroup_delete(pgroup *const Gp);

int group_load(group *const G, FILE *const fp);
#if 0
size_t group_count_subgroups(const group *const G);
void group_compute_subgroups(group *const G);
#endif
void group_delete(group *const G);


#endif