#include "module_checker.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


static void init_diag_eq(module_checker *const mc);
static int valid_diag_eq(const module_checker *const mc);

static void init_diag_sums(module_checker *const mc);

static int valid_diag_mults(const module_checker *const mc, const size_t i);
static int init_diag_mults(module_checker *const mc);
static int next_diag_mults(module_checker *const mc);

static void init_offdiag_eqs(module_checker *const mc);

static int init_offdiag_sums(module_checker *const mc);

static void update_lower_diag_mults(
	module_checker *const mc,
	const size_t i, const size_t j
);
static int valid_offdiag_mults(
	const module_checker *const mc,
	const size_t i, size_t j
);
static int init_offdiag_mults(module_checker *const mc);
static int next_offdiag_mults(module_checker *const mc);

static int check_module_matrix(const module_checker *const mc);
static int check_module_associativity(const module_checker *const mc);
static int construct_modules(module_checker *const mc);


// Allocate all the memory we'll need when looking for modules.
int module_checker_init(
	module_checker *const mc,
	const quadratic_ring *const R
){

	size_t i;
	size_t cur_id;

	const size_t min_orbits =
		(R->d == 0) ? 2 : 1;
	const size_t max_orbits =
		(R->d == 0) ? 2 : 2 + (R->m_sum*R->d)/R->Q->order;

	const size_t inner_factors_bytes =
		R->Q->quotient.num_cosets*sizeof(*mc->inner_factors);

	const size_t H_ids_bytes = max_orbits*sizeof(*mc->H_ids.values);
	const size_t a_row_bytes = max_orbits*sizeof(*mc->a);
	const size_t a_col_bytes = max_orbits*sizeof(**mc->a);
	const size_t a_elem_bytes = R->Q->quotient.num_cosets*sizeof(***mc->a);

	const size_t diag_eq_x_bytes = max_orbits*sizeof(*mc->diag_eq.x);
	const size_t diag_eq_a_bytes = sizeof(*mc->diag_eq.a);
	const size_t diag_eq_a_elem_bytes = max_orbits*sizeof(**mc->diag_eq.a);
	const size_t diag_eq_b_bytes = sizeof(*mc->diag_eq.b);
	const size_t diag_sums_bytes = max_orbits*sizeof(*mc->diag_sums.parts);
	const size_t diag_sums_elem_bytes =
		R->Q->quotient.num_cosets*sizeof(*mc->diag_sums.parts->parts);
	const size_t diag_mults_bytes = max_orbits*sizeof(*mc->diag_mults.perms);

	const size_t max_vars = ((max_orbits - 1)*max_orbits)/2;
	const size_t offdiag_eq_x_bytes = max_vars*sizeof(*mc->offdiag_eqs.x);
	const size_t offdiag_eq_a_bytes = max_orbits*sizeof(*mc->offdiag_eqs.a);
	const size_t offdiag_eq_a_elem_bytes =
		max_vars*sizeof(**mc->offdiag_eqs.a);
	const size_t offdiag_eq_b_bytes = max_orbits*sizeof(*mc->offdiag_eqs.b);
	const size_t offdiag_sums_bytes =
		max_vars*sizeof(*mc->offdiag_sums.parts);
	const size_t offdiag_sums_elem_bytes =
		R->Q->quotient.num_cosets*sizeof(*mc->offdiag_sums.parts->parts);
	const size_t offdiag_mults_bytes =
		max_vars*sizeof(*mc->offdiag_mults.perms);

	const size_t total_bytes =
		inner_factors_bytes +

		H_ids_bytes +
		a_row_bytes + max_orbits*a_col_bytes +
		max_orbits*max_orbits*a_elem_bytes +

		diag_eq_x_bytes + diag_eq_b_bytes +
		diag_eq_a_bytes + diag_eq_a_elem_bytes +
		diag_sums_bytes + max_orbits*diag_sums_elem_bytes +
		diag_mults_bytes +

		offdiag_eq_x_bytes + offdiag_eq_b_bytes +
		offdiag_eq_a_bytes + max_orbits*offdiag_eq_a_elem_bytes +
		offdiag_sums_bytes + max_vars*offdiag_sums_elem_bytes +
		offdiag_mults_bytes;

	if(!pool_init(&mc->pool, total_bytes)){
		puts("Failed to allocate memory for module checker!");
		return 0;
	}

	mc->R = R;

	mc->min_orbits = min_orbits;
	mc->max_orbits = max_orbits;

	// Have our data structures point to the new memory.
	mc->inner_factors     = pool_alloc(&mc->pool, inner_factors_bytes);
	mc->num_inner_factors =
		quadratic_ring_compute_inner_factors(R, mc->inner_factors);

	mc->H_ids.values = pool_alloc(&mc->pool, H_ids_bytes);
	mc->a            = pool_alloc(&mc->pool, a_row_bytes);
	for(i = 0; i < max_orbits; ++i){
		mc->a[i] = pool_alloc(&mc->pool, a_col_bytes);
	}
	for(i = 0; i < max_orbits; ++i){
		size_t j;
		for(j = 0; j < max_orbits; ++j){
			mc->a[i][j] = pool_alloc(&mc->pool, a_elem_bytes);
		}
	}

	mc->diag_eq.x        = pool_alloc(&mc->pool, diag_eq_x_bytes);
	mc->diag_eq.a        = pool_alloc(&mc->pool, diag_eq_a_bytes);
	mc->diag_eq.a[0]     = pool_alloc(&mc->pool, diag_eq_a_elem_bytes);
	mc->diag_eq.b        = pool_alloc(&mc->pool, diag_eq_b_bytes);
	mc->diag_sums.parts  = pool_alloc(&mc->pool, diag_sums_bytes);
	mc->diag_mults.perms = pool_alloc(&mc->pool, diag_mults_bytes);
	for(i = 0; i < max_orbits; ++i){
		mc->diag_sums.parts[i].parts = pool_alloc(
			&mc->pool, diag_sums_elem_bytes
		);
		mc->diag_mults.perms[i].values = mc->a[i][i];
	}

	mc->offdiag_eqs.x        = pool_alloc(&mc->pool, offdiag_eq_x_bytes);
	mc->offdiag_eqs.a        = pool_alloc(&mc->pool, offdiag_eq_a_bytes);
	for(i = 0; i < max_orbits; ++i){
		mc->offdiag_eqs.a[i] = pool_alloc(&mc->pool, offdiag_eq_a_elem_bytes);
	}
	mc->offdiag_eqs.b        = pool_alloc(&mc->pool, offdiag_eq_b_bytes);
	mc->offdiag_sums.parts   = pool_alloc(&mc->pool, offdiag_sums_bytes);
	mc->offdiag_mults.perms  = pool_alloc(&mc->pool, offdiag_mults_bytes);
	cur_id = 0;
	for(i = 0; i < max_orbits; ++i){
		size_t j;
		for(j = i + 1; j < max_orbits; ++j){
			mc->offdiag_sums.parts[cur_id].parts = pool_alloc(
				&mc->pool, offdiag_sums_elem_bytes
			);
			mc->offdiag_mults.perms[cur_id].values = mc->a[i][j];
			++cur_id;
		}
	}

	return 1;
}

int module_checker_find_modules(const quadratic_ring *const R){
	module_checker mc;
	if(!module_checker_init(&mc, R)){
		return 0;
	}

	// For a given number of G-orbits o, we start by constructing
	// o potentially Morita equivalent algebra objects as follows.
	//     1. For each i in {1, 2, ..., o}, choose a subgroup H_i.
	//
	//     2. For each subgroup H_i, choose a_{i,i,g\psi(H_i Q)H_i Q}.
	//
	//     3. Check that the squared dimensions of the corresponding
	//        module objects sum to the dimension of the fusion category.
	// After this, we have the simple objects of a potential module
	// category along with the action of G. All that remains is to
	// determine the action of the non-invertible objects, which we
	// do in the construct_modules() function.
	mc.num_orbits = mc.min_orbits;
	for(; mc.num_orbits <= mc.max_orbits; ++mc.num_orbits){
		tuple_start(&mc.H_ids, 0, R->G.num_subgroups - 1, mc.num_orbits);
		// Iterate over every combination of subgroups.
		do {
			printf("H = ");
			tuple_print(&mc.H_ids);

			init_diag_eq(&mc);
			if(!linear_diophantine_start(&mc.diag_eq)){
				continue;
			}
			// Iterate over every solution to the multiplicity equation.
			do {
				if(!valid_diag_eq(&mc)){
					continue;
				}
				init_diag_sums(&mc);
				// Iterate over every set of multiplicities
				// that sum to the solution we found above.
				do {
					if(!init_diag_mults(&mc)){
						continue;
					}
					// Iterate over every permutation of the multiplicities.
					do {
						if(!construct_modules(&mc)){
							module_checker_delete(&mc);
							return 0;
						}
					} while(next_diag_mults(&mc));
				} while(partition_array_next(&mc.diag_sums));
			} while(linear_diophantine_next(&mc.diag_eq));
		} while(tuple_next(&mc.H_ids));
	}

	module_checker_delete(&mc);
	return 1;
}

void module_checker_delete(module_checker *const mc){
	pool_delete(&mc->pool);
}


// The squared dimensions of the simple objects in the module category
// should sum to the dimension of the quadratic category. This gives
// a linear Diophantine in the squared sums of the a_{i,i}, which is
// slightly nicer to solve than brute force.
static void init_diag_eq(module_checker *const mc){
	const quadratic_ring *const R = mc->R;

	size_t i;
	// Store a[i] = |(\psi(H_i Q)H_i Q)/H_i|.
	for(i = 0; i < mc->num_orbits; ++i){
		const size_t H_id = mc->H_ids.values[i];
		const size_t psi_H_id = R->psi_ids[H_id][H_id];
		mc->diag_eq.a[0][i] =
			R->G.subgroups[psi_H_id].order/R->G.subgroups[H_id].order;
	}

	mc->diag_eq.num_vars = mc->num_orbits;
	// Compute the right-hand side. When d = FPdim(X_Q)
	// is irrational, we store its value as 0.
	mc->diag_eq.b[0] = (R->d == 0) ? R->m_sum :
		(2 - mc->num_orbits)*R->Q->order/R->d + R->m_sum;
	mc->diag_eq.num_eqs = 1;
}

// If we have multiple orbits corresponding to the same subgroup,
// multiple solutions to the Diophantine may correspond to the
// same module. To fix this, we ask that the multiplicity sums
// for identical subgroups are in ascending order.
static int valid_diag_eq(const module_checker *const mc){
	size_t i;
	for(i = 1; i < mc->num_orbits; ++i){
		// This is a pretty bad way of doing this, it'd be
		// better to handle the case of multiple identical
		// subgroups before going to the Diophantine.
		if(mc->H_ids.values[i] == mc->H_ids.values[i - 1]){
			if(mc->diag_eq.x[i] < mc->diag_eq.x[i - 1]){
				return 0;
			}
		}
	}
	return 1;
}

static void init_diag_sums(module_checker *const mc){
	const quadratic_ring *const R = mc->R;

	size_t i;
	// Initialize each partition.
	for(i = 0; i < mc->num_orbits; ++i){
		const size_t H_id = mc->H_ids.values[i];
		partition_start(
			&mc->diag_sums.parts[i], mc->diag_eq.x[i],
			R->G.subgroups[
				R->psi_ids[H_id][H_id]
			].quotient.num_cosets
		);
	}

	mc->diag_sums.num_parts = mc->num_orbits;
}

// Return whether the permutation a_{i,i} is a valid set
// of multiplicities for the algebra object candidate
//     A_i = H_i \oplus \bigoplus_{gQ \in G/Q}{X_{gQ}^{\oplus a_{i,i,g\psi(H_i Q)H_i Q}}}.
// There are two conditions that this function checks.
//     1. If A_i is the internal end of an object in a module we've
//        already seen before, we can skip it. This is the case if
//        a_{i,i} is an inner automorphism of a permutation that
//        precedes it in lexicographic order.
//
//     2. We require that for all gH_i Q \in G/(H_i Q),
//            a_{i,i,g\psi(H_i Q)H_i Q} = a_{i,i,\psi(g^{-1}H_i Q)H_i Q}
//                                      = a_{i,i,r(Q)\psi(g^{-1}H_i Q)H_i Q},
//        The first equality comes from symmetry of the fusion matrix,
//        while the second comes from associativity of the module.
static int valid_diag_mults(const module_checker *const mc, const size_t i){
	const quadratic_ring *const R = mc->R;
	const size_t Hi_id = mc->H_ids.values[i];
	const subgroup *const psi_Hii =
		&R->G.subgroups[R->psi_ids[Hi_id][Hi_id]];
	const subgroup *const Hi_Q =
		&R->G.subgroups[R->G.product_ids[Hi_id][R->Q_id]];
	const integer_t *a = mc->diag_mults.perms[i].values;

	size_t g1_coset;
	// Check the first condition.
	for(g1_coset = 0; g1_coset < mc->num_inner_factors; ++g1_coset){
		const group_elem_id g1 = mc->inner_factors[g1_coset];

		size_t g2_coset;
		for(g2_coset = 0; g2_coset < psi_Hii->quotient.num_cosets; ++g2_coset){
			// a_{i,i,g_1^{-1}\psi(g_1 Q)g_2\psi(H_i Q)H_i Q}
			const integer_t a_aut = a[psi_Hii->quotient.coset_ids[
				R->G.multiply[g1][psi_Hii->quotient.cosets[g2_coset][0]]
			]];
			// If the inner automorphism is smaller than a_{i,i}
			// in the first place that they differ, it precedes it
			// lexicographically. Thus, we've already seen it.
			if(a_aut < a[g2_coset]){
				return 0;
			}else if(a_aut > a[g2_coset]){
				break;
			}
		}
	}

	// Check the second condition.
	for(g1_coset = 0; g1_coset < Hi_Q->quotient.num_cosets; ++g1_coset){
		const group_elem_id g1 = Hi_Q->quotient.cosets[g1_coset][0];
		// \psi(g_1^{-1}H_i Q)
		const group_elem_id psi_g1_inv = R->psi[R->G.inverse[g1]];
		const integer_t ag1 = a[psi_Hii->quotient.coset_ids[g1]];
		if(
			a[psi_Hii->quotient.coset_ids[psi_g1_inv]] != ag1 ||
			a[psi_Hii->quotient.coset_ids[
				R->G.multiply[R->r][psi_g1_inv]
			]] != ag1
		){

			return 0;
		}
	}

	return 1;
}

static int init_diag_mults(module_checker *const mc){
	size_t i;
	// Initialize each permutation.
	for(i = 0; i < mc->diag_sums.num_parts; ++i){
		const partition *const cur_part = &mc->diag_sums.parts[i];
		permutation *const cur_perm = &mc->diag_mults.perms[i];
		permutation_start(
			cur_perm, cur_part->max_length,
			cur_part->parts, cur_part->length
		);
		// Make sure the starting permutations are valid.
		// We only really need to check the second condition.
		while(!valid_diag_mults(mc, i)){
			if(!permutation_next(cur_perm)){
				return 0;
			}
		}
	}

	mc->diag_mults.num_perms = mc->diag_sums.num_parts;

	return 1;
}

static int next_diag_mults(module_checker *const mc){
	size_t i;
	for(i = 0; i < mc->diag_mults.num_perms; ++i){
		permutation *const cur_perm = &mc->diag_mults.perms[i];
		if(permutation_next(cur_perm)){
			/** We're technically checking a lot of permutations **/
			/** multiple times here. We should compute a list of **/
			/** valid permutations for each orbit and then only  **/
			/** iterate over tuples of those.                    **/
			if(valid_diag_mults(mc, i)){
				return 1;
			}
			// If the current multiplicities aren't valid,
			// we need to increment this permutation again.
			--i;

		// Don't worry about restarting the last permutation.
		}else if(i < mc->diag_mults.num_perms - 1){
			permutation_restart(cur_perm);
		}
	}

	return 0;

	#if 0
	/** We should eventually do something like this to **/
	/** fix the issue above. In fact, it might be good **/
	/** to do this earlier in case we have tuples with **/
	/** with the same groups and sums appearing later. **/

	/** NOTE: This is outdated! There may be stuff we  **/
	/** don't check down here that we do check above.  **/
	/*
	/// Everything except the vectors' buffers can be allocated in the pool.
	vector *valid_perms = malloc(num_orbits*sizeof(*valid_perms));
	size_t *num_valid_perms = malloc(num_orbits*sizeof(*num_valid_perms));
	permutation perm;
	size_t *cur_perm;

	perm.values = malloc(R->G.order*sizeof(*cur_perm.values));

	size_t i;
	for(i = 0; i < part_array.num_parts; ++i){
		const partition *const cur_part = &part_array.parts[i];
		const size_t H_id = H_ids[i];
		const subgroup *const psi_H = &R->G.subgroups[mc->psi_ids[H_id][H_id]];

		permutation_start(
			&perm, psi_H->quotient.num_cosets,
			cur_part->parts, cur_part->length
		);
		vector_init(&valid_perms[i], perm.length*sizeof(*perm.values));
		num_valid_perms[i] = 0;

		// Store any permutations that
		do {
			if(valid_algebra_multiplicities(
				R, inner_factors, num_inner_factors, psi_H, perm.values
			)){

				vector_add(&valid_perms[i], perm.values);
				++num_valid_perms[i];
			}
		} while(permutation_next(&perm));

		// If we found no valid permutations for this
		// partition, go to the next tuple of groups.
		if(num_valid_perms[i] == 0){
			return 0;
		}
	}

	cur_perm = malloc(num_orbits*sizeof(*cur_perm));
	memset(cur_perm, 0, num_orbits*sizeof(*cur_perm));
	for(i = 0; i < num_orbits; ++i){
		if(cur_perm[i] < num_valid_perms[i] - 1){
			++cur_perm[i];
			return 1;
		}else{
			cur_perm[i] = 0;
		}
	}
	*/
	#endif
}

// By summing the fusion matrix compatibility equations
// over each gQ \in G/Q, we get a linear Diophantine in
// the squared sums of the a_{i,j}. This is a slightly
// nicer way of determining them versus brute force.
static void init_offdiag_eqs(module_checker *const mc){
	const quadratic_ring *const R = mc->R;

	size_t i;

	mc->offdiag_eqs.num_vars =
		((mc->num_orbits - 1)*mc->num_orbits)/2;

	for(i = 0; i < mc->num_orbits; ++i){
		const size_t Hi_id = mc->H_ids.values[i];
		const subgroup *Hi = &R->G.subgroups[Hi_id];
		const subgroup *psi_Hii =
			&R->G.subgroups[R->psi_ids[Hi_id][Hi_id]];

		// We store the variables in the order
		//     a_{0,1}, a_{0,2}, ..., a_{0,o-1},
		//     a_{1,2}, a_{1,3}, ..., a_{1,o-1},
		//     ...
		//     a_{o-3,o-2}, a_{o-3,o-1},
		//     a_{o-2,o-1}.
		// Therefore, for i < j, we have
		//     a_{i,j} = offdiag_eqs.x[o*i - (i*(i + 1))/2 + (j - (i + 1))].
		size_t cur_id = i - 1;
		size_t j;

		// The coefficients for most of the variables before
		// a_{i,i+1} are zero, so we just set them all here.
		memset(
			&mc->offdiag_eqs.a[i][0], 0,
			mc->num_orbits*i - (i*(i + 1))/2
		);
		// Set the non-zero coefficients for a_{j,i} when j < i.
		for(j = 0; j < i; ++j){
			const size_t Hj_id = mc->H_ids.values[j];
			const subgroup *Hj = &R->G.subgroups[Hj_id];
			const subgroup *psi_Hij =
				&R->G.subgroups[R->psi_ids[Hi_id][Hj_id]];

			// |(\psi(H_i Q)H_j Q)/Q||(\psi(H_i Q)H_j Q)/H_j|
			mc->offdiag_eqs.a[i][cur_id] =
				(psi_Hij->order/R->Q->order)*(psi_Hij->order/Hj->order);
			cur_id += mc->num_orbits - j - 2;
		}

		// After the previous loop, cur_id contains the
		// index of a_{i-1,o-1}, so we should increment it.
		++cur_id;
		// Set the coefficients for a_{i,j} when j > i.
		for(j = i + 1; j < mc->num_orbits; ++j){
			const size_t Hj_id = mc->H_ids.values[j];
			const subgroup *Hj = &R->G.subgroups[Hj_id];
			const subgroup *psi_Hij =
				&R->G.subgroups[R->psi_ids[Hi_id][Hj_id]];

			// |(\psi(H_i Q)H_j Q)/Q||(\psi(H_i Q)H_j Q)/H_j|
			mc->offdiag_eqs.a[i][cur_id] =
				(psi_Hij->order/R->Q->order)*(psi_Hij->order/Hj->order);
			++cur_id;
		}
		// The coefficients for the variables after a_{i,o-1}
		// are all zero, so we just set them all here.
		memset(
			&mc->offdiag_eqs.a[i][cur_id], 0,
			mc->offdiag_eqs.num_vars - cur_id
		);

		// Set the right-hand side of the equation.
		mc->offdiag_eqs.b[i] =
			Hi->order +
			(psi_Hii->order/R->Q->order)*mc->diag_eq.x[i]*(
				R->m_sum - mc->diag_eq.a[0][i]*mc->diag_eq.x[i]
			);
	}

	mc->offdiag_eqs.num_eqs = mc->num_orbits;
}

static int init_offdiag_sums(module_checker *const mc){
	const quadratic_ring *const R = mc->R;

	size_t cur_id = 0;
	size_t i;
	// Initialize each partition.
	for(i = 0; i < mc->num_orbits - 1; ++i){
		const size_t Hi_id = mc->H_ids.values[i];

		size_t j;
		for(j = i + 1; j < mc->num_orbits; ++j){
			const subgroup *const psi_Hij =
				&R->G.subgroups[R->psi_ids[Hi_id][mc->H_ids.values[j]]];

			const integer_t cur_sum =
				(integer_t)sqrt(mc->offdiag_eqs.x[cur_id]);
			// Make sure the sum is a square!
			if(cur_sum*cur_sum != mc->offdiag_eqs.x[cur_id]){
				return 0;
			}

			partition_start(
				&mc->offdiag_sums.parts[cur_id], cur_sum,
				psi_Hij->quotient.num_cosets
			);
			++cur_id;
		}
	}

	mc->offdiag_sums.num_parts = mc->offdiag_eqs.num_vars;

	return 1;
}

// We only solve the Diophantine for the a_{i,j} with i < j.
// To get the lower diagonal of the matrix, we use the fact that
//     a_{j,i,g\psi(H_j Q)H_i Q} = a_{i,j,r(Q)\psi(g^{-1}H_i Q)H_j Q}.
static void update_lower_diag_mults(
	module_checker *const mc,
	const size_t i, const size_t j
){

	const quadratic_ring *const R = mc->R;
	const size_t Hi_id = mc->H_ids.values[i];
	const size_t Hj_id = mc->H_ids.values[j];
	const subgroup *const psi_Hji =
		&R->G.subgroups[R->psi_ids[Hj_id][Hi_id]];
	const subgroup *const psi_Hij =
		&R->G.subgroups[R->psi_ids[Hi_id][Hj_id]];

	size_t g_coset;
	for(g_coset = 0; g_coset < psi_Hji->quotient.num_cosets; ++g_coset){
		const group_elem_id g = psi_Hji->quotient.cosets[g_coset][0];
		mc->a[j][i][g_coset] = mc->a[i][j][psi_Hij->quotient.coset_ids[
			// r(Q)psi(g^{-1})
			R->G.multiply[R->r][R->psi[R->G.inverse[g]]]
		]];
	}
}

// Return whether the permutation a_{i,j} is a valid set
// of multiplicities for the module object associated to
//     A_i = H_i \oplus \bigoplus_{gQ \in G/Q}{X_{gQ}^{\oplus a_{i,i,g\psi(H_i Q)H_i Q}}}.
// This time, we only need to check that this doesn't
// correspond to a module object that we've seen before.
static int valid_offdiag_mults(
	const module_checker *const mc,
	const size_t i, size_t j
){

	const quadratic_ring *const R = mc->R;
	const size_t Hi_id = mc->H_ids.values[i];
	const subgroup *const Hi = &R->G.subgroups[Hi_id];
	const subgroup *const psi_Hii =
		&R->G.subgroups[R->psi_ids[Hi_id][Hi_id]];
	integer_t ***const a = mc->a;

	for(; j < mc->num_orbits; ++j){
		const subgroup *const psi_Hij =
			&R->G.subgroups[R->psi_ids[Hi_id][mc->H_ids.values[j]]];
		size_t gi_coset;

		// We already checked this case when
		// constructing the algebra objects!
		if(i == j){
			continue;
		}

		// For every module object in this G-orbit, we need to check
		// whether the new multiplicities make it an object from a
		// module we've already seen before.
		for(gi_coset = 0; gi_coset < Hi->quotient.num_cosets; ++gi_coset){
			const group_elem_id gi = Hi->quotient.cosets[gi_coset][0];
			// \psi(g_i H_i Q)
			const group_elem_id psi_gi = R->psi[gi];
			// g_i^{-1}\psi(g_i H_i Q)
			const group_elem_id gi_inv_psi_gi =
				R->G.multiply[R->G.inverse[gi]][psi_gi];

			size_t gij_coset;
			// Check whether the multiplicities for the current module
			// object precede our new ones in lexicographic order.
			for(gij_coset = 0; gij_coset < psi_Hij->quotient.num_cosets; ++gij_coset){
				const group_elem_id gij = psi_Hij->quotient.cosets[gij_coset][0];
				// a_{i,j,\psi(g_i H_i Q)g_{ij}H_i Q}
				const integer_t aij_aut = a[i][j][
					psi_Hij->quotient.coset_ids[R->G.multiply[psi_gi][gij]]
				];

				// If the current module object's multiplicities do
				// come earlier, we may have seen this module before.
				// We'll need to make sure the internal end of the
				// module object matches the one we're looking at.
				if(aij_aut < a[i][j][gij_coset]){
					size_t gii_coset;

					// Early exit if the internal end is trivially the same.
					if(psi_Hii->quotient.coset_ids[gi_inv_psi_gi] == 0){
						return 0;
					}

					for(gii_coset = 0; gii_coset < psi_Hii->quotient.num_cosets; ++gii_coset){
						const integer_t aii_aut = a[i][i][
							psi_Hii->quotient.coset_ids[R->G.multiply[
								gi_inv_psi_gi
							][
								psi_Hii->quotient.cosets[gii_coset][0]
							]]
						];

						// If this module object has a different
						// internal end, keep searching.
						if(aii_aut != a[i][i][gii_coset]){
							break;
						}
					}
					// If the fusion matrix is a match, we
					// can safely ignore the new permutation.
					if(gii_coset == psi_Hii->quotient.num_cosets){
						return 0;
					}
					// Otherwise, there's no point continuing with
					// this permutation, so move on to the next one.
					break;

				// If they come after, keep searching.
				}else if(aij_aut > a[i][j][gij_coset]){
					break;
				}
			}
		}
	}

	return 1;
}

static int init_offdiag_mults(module_checker *const mc){
	size_t i;
	size_t cur_id = 0;
	for(i = 0; i < mc->num_orbits - 1; ++i){
		size_t j;
		for(j = i + 1; j < mc->num_orbits; ++j){
			const partition *const cur_part = &mc->offdiag_sums.parts[cur_id];
			permutation *const cur_perm = &mc->offdiag_mults.perms[cur_id];
			// Since we iterate over permutations in lexicographic
			// order, it's impossible for the first one to come
			// from one that appears earlier. Hence, we don't need
			// to check validity here.
			permutation_start(
				cur_perm, cur_part->max_length,
				cur_part->parts, cur_part->length
			);
			// Update the entries of the
			// matrix below the diagonal.
			update_lower_diag_mults(mc, i, j);

			++cur_id;
		}
		if(i == 0){
			break;
		}
	}

	mc->offdiag_mults.num_perms = mc->offdiag_sums.num_parts;

	return 1;
}

static int next_offdiag_mults(module_checker *const mc){
	size_t i;
	size_t cur_id = 0;
	for(i = 0; i < mc->num_orbits - 1; ++i){
		size_t j;
		for(j = i + 1; j < mc->num_orbits; ++j){
			permutation *const cur_perm = &mc->offdiag_mults.perms[cur_id];
			if(permutation_next(cur_perm)){
				// Update the entries of the
				// matrix below the diagonal.
				update_lower_diag_mults(mc, i, j);

				/** We're technically checking a lot of permutations   **/
				/** multiple times here. We should compute a list of   **/
				/** valid permutations for each variable and then only **/
				/** iterate over tuples of those.                      **/
				if(valid_offdiag_mults(mc, i, j)){
					return 1;
				}
				// If the current multiplicities aren't valid,
				// we need to increment this permutation again.
				--j;

			// If we're not yet on the last variable,
			// restart it and move to the next.
			}else if(cur_id < mc->offdiag_mults.num_perms - 1){
				permutation_restart(cur_perm);
				update_lower_diag_mults(mc, i, j);
				++cur_id;
			}
		}
	}

	// If we've reached the last permutation of the last
	// variable, we've completely exhausted all permutations.
	return 0;
}

// Make sure that the fusion matrices are compatible in
// the sense that F^M (F^M)^T = F^A. We only need to
// check the rows not coming from the pointed subcategory.
static int check_module_matrix(const module_checker *const mc){
	const quadratic_ring *const R = mc->R;
	integer_t ***const a = mc->a;

	size_t i;
	for(i = 0; i < mc->num_orbits; ++i){
		const size_t Hi_id = mc->H_ids.values[i];
		const subgroup *const Hi = &R->G.subgroups[Hi_id];
		const subgroup *const psi_Hii =
			&R->G.subgroups[R->psi_ids[Hi_id][Hi_id]];

		size_t gi_coset;
		for(gi_coset = 0; gi_coset < R->Q->quotient.num_cosets; ++gi_coset){
			const group_elem_id gi = R->Q->quotient.cosets[gi_coset][0];
			size_t gj_coset;
			size_t j;

			integer_t lhs;
			integer_t rhs;

			// Compute the left-hand side.
			lhs = group_subsets_intersection(
				// H_i
				Hi->elem_ids, Hi->order,
				// \psi(g_i Q)
				R->Q->quotient.cosets[R->Q->quotient.coset_ids[R->psi[gi]]],
				R->Q->order
			);
			for(gj_coset = 0; gj_coset < R->Q->quotient.num_cosets; ++gj_coset){
				const group_elem_id gj = R->Q->quotient.cosets[gj_coset][0];

				// a_{i,i,g_i g_j\psi(H_i Q)H_i Q}m_{g_j Q}
				lhs += a[i][i][psi_Hii->quotient.coset_ids[
					R->G.multiply[gi][gj]
				]]*R->m[gj];
			}

			// Compute the right-hand side.
			rhs = 0;
			for(j = 0; j < mc->num_orbits; ++j){
				const size_t Hj_id = mc->H_ids.values[j];
				const subgroup *const Hj = &R->G.subgroups[Hj_id];
				const subgroup *const psi_Hij =
					&R->G.subgroups[R->psi_ids[Hi_id][Hj_id]];

				for(gj_coset = 0; gj_coset < Hj->quotient.num_cosets; ++gj_coset){
					const group_elem_id gj = Hj->quotient.cosets[gj_coset][0];

					// a_{i,j,g_i g_j\psi(H_i Q)H_j Q}a_{i,j,g_j\psi(H_i Q)H_j Q}
					rhs += a[i][j][psi_Hij->quotient.coset_ids[
						R->G.multiply[gi][gj]
					]]*a[i][j][psi_Hij->quotient.coset_ids[gj]];
				}
			}

			if(lhs != rhs){
				return 0;
			}
		}
	}

	return 1;
}

// Check associativity of the new module. We only need to check
//     (X_Q \otimes X_Q) \otimes M_{i,gH_i} \cong X_Q \otimes (X_Q \otimes M_{i,gH_i}),
// for each gH_i \in G/H_i. This gives us an equation for each
// G-orbit j and module object M_{j,g'H_j} in the orbit.
static int check_module_associativity(const module_checker *const mc){
	const quadratic_ring *const R = mc->R;
	integer_t ***const a = mc->a;

	size_t i;
	for(i = 0; i < mc->num_orbits; ++i){
		const size_t Hi_id = mc->H_ids.values[i];
		const subgroup *const Hi = &R->G.subgroups[Hi_id];

		size_t gi_coset;
		for(gi_coset = 0; gi_coset < Hi->quotient.num_cosets; ++gi_coset){
			const group_elem_id gi = Hi->quotient.cosets[gi_coset][0];

			size_t j;
			for(j = 0; j < mc->num_orbits; ++j){
				const size_t Hj_id = mc->H_ids.values[j];
				const subgroup *const Hj = &R->G.subgroups[Hj_id];
				const subgroup *const psi_Hij =
					&R->G.subgroups[R->psi_ids[Hi_id][Hj_id]];

				size_t gj_coset;
				for(gj_coset = 0; gj_coset < Hj->quotient.num_cosets; ++gj_coset){
					const group_elem_id gj = Hj->quotient.cosets[gj_coset][0];
					size_t gk_coset;
					size_t k;

					integer_t lhs;
					integer_t rhs;

					// Compute the left-hand side.
					if(i == j){
						lhs = group_subsets_intersection(
							// g_i^{-1}g_j H_i
							Hi->quotient.cosets[Hi->quotient.coset_ids[
								R->G.multiply[R->G.inverse[gi]][gj]
							]], Hi->order,
							// r(Q)^{-1}
							R->Q->quotient.cosets[R->Q->quotient.coset_ids[
								R->G.inverse[R->r]
							]], R->Q->order
						);
					}else{
						lhs = 0;
					}
					for(gk_coset = 0; gk_coset < R->Q->quotient.num_cosets; ++gk_coset){
						const group_elem_id gk = R->Q->quotient.cosets[gk_coset][0];

						// m_{r(Q)g_k Q}a_{i,j,\psi(g_i H_i Q)g_j^{-1}g_k H_j Q}
						lhs += R->m[
							R->G.multiply[R->r][gk]
						]*a[i][j][psi_Hij->quotient.coset_ids[R->G.multiply[
							R->psi[gi]
						][
							R->G.multiply[R->G.inverse[gj]][gk]
						]]];
					}

					// Compute the right-hand side.
					rhs = 0;
					for(k = 0; k < mc->num_orbits; ++k){
						const size_t Hk_id = mc->H_ids.values[k];
						const subgroup *const Hk = &R->G.subgroups[Hk_id];
						const subgroup *const psi_Hik =
							&R->G.subgroups[R->psi_ids[Hi_id][Hk_id]];
						const subgroup *const psi_Hkj =
							&R->G.subgroups[R->psi_ids[Hk_id][Hj_id]];

						for(gk_coset = 0; gk_coset < Hk->quotient.num_cosets; ++gk_coset){
							const group_elem_id gk = Hk->quotient.cosets[gk_coset][0];

							// a_{i,k,\psi(g_i H_i Q)g_k^{-1}H_k Q}a_{k,j,g_j^{-1}\psi(g_k H_k Q)H_j Q}
							rhs += a[i][k][psi_Hik->quotient.coset_ids[
								R->G.multiply[R->psi[gi]][R->G.inverse[gk]]
							]]*a[k][j][psi_Hkj->quotient.coset_ids[
								R->G.multiply[R->G.inverse[gj]][R->psi[gk]]
							]];
						}
					}

					if(lhs != rhs){
						return 0;
					}
				}
			}
		}
	}

	return 1;
}

// Try and build modules that give a Morita equivalence between the algebras
//     A_i = H_i \oplus \bigoplus_{gQ \in G/Q}{X_{gQ}^{\oplus a_{i,i,g\psi(H_i Q)H_i Q}}}.
// The algebra objects themselves give us the simple objects and the action
// of G, so we must now determine the action of the non-invertible objects.
static int construct_modules(module_checker *const mc){
	init_offdiag_eqs(mc);
	if(!linear_diophantine_start(&mc->offdiag_eqs)){
		return 1;
	}
	do {
		// Solutions to the Diophantine system represent the
		// squared sums of the a_{i,j}. If any of them aren't
		// squares, we can continue on to the next possibility.
		if(!init_offdiag_sums(mc)){
			continue;
		}
		do {
			// Setup offdiag_mults.
			if(!init_offdiag_mults(mc)){
				continue;
			}
			do {
				// If everything was satisfied, we've found a new module!
				if(check_module_matrix(mc) && check_module_associativity(mc)){
					quadratic_module M;
					if(!quadratic_module_init(
						&M, mc->R, mc->H_ids.values, mc->a, mc->num_orbits
					)){

						return 0;
					}
					quadratic_module_print(&M, mc->R);
					quadratic_module_delete(&M);
					puts("Found module!");
					{
						size_t i;
						for(i = 0; i < mc->num_orbits; ++i){
							size_t j;
							const size_t Hi_id = mc->H_ids.values[i];
							const subgroup *const psi_Hii = &mc->R->G.subgroups[
								mc->R->psi_ids[Hi_id][Hi_id]
							];

							printf("\tH_%u + {%u", Hi_id, mc->a[i][i][0]);

							for(j = 1; j < mc->R->Q->quotient.num_cosets; ++j){
								printf(", %u", mc->a[i][i][
									psi_Hii->quotient.coset_ids[
										mc->R->Q->quotient.cosets[j][0]
									]
								]);
							}
							printf("}X\n");
						}
					}
				}
			} while(next_offdiag_mults(mc));
		} while(partition_array_next(&mc->offdiag_sums));
	} while(linear_diophantine_next(&mc->offdiag_eqs));

	return 1;
}