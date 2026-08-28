#include "hi_ring.h"


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "partition.h"
#include "permutation.h"


static void print_module_matrix(
	const group *const G,
	const subgroup *const H1, const integer_t *const a,
	const subgroup *const H2, const integer_t *const x
);
static size_t compute_inner_factors(
	const hi_ring *const R, group_elem_id **inner_factors
);
static int valid_algebra_multiplicities(
	const hi_ring *const R,
	const group_elem_id *const inner_factors, const size_t num_inner_factors,
	const subgroup *const H, const integer_t *const a
);
static int valid_module_multiplicities(
	const hi_ring *const R,
	const subgroup *const H1, const subgroup *const H2,
	const integer_t *const a1, const integer_t *const x1
);
static int valid_module_associativity(
	const hi_ring *const R,
	const subgroup *const H1, const subgroup *const H2,
	const integer_t *const a1, const integer_t *const a2,
	const integer_t *const x1, const integer_t *const x2
);
static int valid_module_rigidity(
	const hi_ring *const R,
	const subgroup *const H1, const subgroup *const H2,
	const integer_t *const a1, const integer_t *const a2,
	const integer_t *const x1, const integer_t *const x2
);
static void construct_modules(
	hi_ring *const R,
	integer_t *const scratchpad, const integer_t a1_sum,
	const subgroup *const H1, const integer_t *const a1,
	const subgroup *const H2, const integer_t *const a2
);


int hi_ring_load(hi_ring *const R, const char *const path){
	group_elem_id g;
	long unsigned int scan_val;

	FILE *const fp = fopen(path, "r");
	if(fp == NULL){
		printf("Failed to open file \"%s\"!\n", path);
		return 0;
	}

	// Initialize the structure. This is mostly
	// just so that if anything fails, we can
	// use our normal deletion function.
	R->psi = NULL;
	R->modules = NULL;
	R->num_modules = 0;

	// Load the group of invertibles and all of its subgroups.
	if(!group_load(&R->G, fp)){
		printf("Failed to load group from file \"%s\"!", path);
		fclose(fp);
		return 0;
	}
	// Load the index of the subgroup Q.
	fscanf(fp, "%lu", &scan_val);
	if(scan_val < R->G.num_subgroups){
		R->Q = &R->G.subgroups[scan_val];
	}else{
		printf(
			"Given subgroup index %lu, but only "
			"%lu subgroups were specified! Using "
			"the subgroup at index 0.\n",
			scan_val, (long unsigned int)R->G.num_subgroups
		);
		R->Q = &R->G.subgroups[0];
	}

	R->psi = malloc(R->G.order*sizeof(*R->psi));
	if(R->psi == NULL){
		puts("Failed to allocate memory for group automorphism!");
		fclose(fp);
		return 0;
	}
	// Load the group automorphism \psi.
	for(g = 0; g < R->G.order; ++g){
		fscanf(fp, "%lu", &scan_val);
		R->psi[g] = scan_val;
	}

	fclose(fp);
	return 1;
}

void hi_ring_print_module(const hi_ring *const R, const hi_module *const M){
	print_module_matrix(&R->G, M->H1, M->a1, M->H2, M->x1);
	printf("\n");
	print_module_matrix(&R->G, M->H2, M->a2, M->H1, M->x2);
	printf("\n\n");
}

int hi_ring_find_modules(hi_ring *const R){
	size_t H1_id;
	integer_t *scratchpad;
	partition p1, p2;
	permutation a1, a2;

	group_elem_id *inner_factors;
	const size_t num_inner_factors =
		compute_inner_factors(R, &inner_factors);
	if(num_inner_factors == -1){
		return 0;
	}

	// Allocate memory for all of our partitions and
	// permutations. Each one has at most |G| integers.
	scratchpad = malloc(7*R->G.order*sizeof(*scratchpad));
	if(scratchpad == NULL){
		if(inner_factors != NULL){
			free(inner_factors);
		}
		return 0;
	}
	// Set up the pointers for our partitions and permutations.
	p1.parts  = scratchpad;
	p2.parts  = p1.parts  + R->G.order;
	a1.values = p2.parts  + R->G.order;
	a2.values = a1.values + R->G.order;

	for(H1_id = 0; H1_id < R->G.num_subgroups; ++H1_id){
		size_t H2_id;
		const subgroup *const H1 = &R->G.subgroups[H1_id];

		printf("Starting subgroup H_1 = %lu\n", (long unsigned int)H1_id);

		// Start from the first subgroup so we
		// aren't checking every pair twice.
		for(H2_id = H1_id; H2_id < R->G.num_subgroups; ++H2_id){
			const subgroup *const H2 = &R->G.subgroups[H2_id];
			integer_t a1_sum;
			// Avoid checking every possibility twice when H_1 = H_2.
			const integer_t sum_max = (H1 == H2) ? R->Q->order/2 : R->Q->order;

			printf("Starting subgroup H_2 = %lu\n", (long unsigned int)H2_id);

			// To compute the multiplicities, we use the dimension identity,
			//     \sum_{gH_1 \in G/H_1}{a_1[gH_1]} + \sum_{gH_2 \in G/H_2}{a_2[gH_2]} = |Q|.
			// We start by choosing the total value of the first sum.
			for(a1_sum = 0; a1_sum <= sum_max; ++a1_sum){
				// Now choose a partition of this value. This is our
				// array of multiplicities a_1 (up to permutations).
				partition_start(
					&p1, a1_sum, H1->quotient.num_cosets
				);
				do {
					// The value of the second sum is |Q| - a1_sum,
					// so we next choose a partition of that too.
					partition_start(
						&p2, R->Q->order - a1_sum, H2->quotient.num_cosets
					);
					do {
						permutation_start(
							&a1, H1->quotient.num_cosets,
							p1.parts, p1.length
						);
						do {
							// Ignore this object if we've seen it before.
							if(!valid_algebra_multiplicities(
								R,
								inner_factors, num_inner_factors,
								H1, a1.values
							)){

								continue;
							}

							permutation_start(
								&a2, H2->quotient.num_cosets,
								p2.parts, p2.length
							);
							do {
								// Ignore this object if we've seen it before.
								if(!valid_algebra_multiplicities(
									R,
									inner_factors, num_inner_factors,
									H2, a2.values
								)){

									continue;
								}

								// This is a genuinely new candidate, so
								// check to see whether it gives any modules!
								construct_modules(
									R,
									a2.values + R->G.order, a1_sum,
									H1, a1.values, H2, a2.values
								);
							} while(permutation_next(&a2));
						} while(permutation_next(&a1));
					} while(partition_next(&p2));
				} while(partition_next(&p1));
			}
		}
	}

	free(scratchpad);
	if(inner_factors != NULL){
		free(inner_factors);
	}

	return 1;
}

void hi_ring_delete(hi_ring *const R){
	free(R->modules);
	free(R->psi);
	group_delete(&R->G);
}


static void print_module_matrix(
	const group *const G,
	const subgroup *const H1, const integer_t *const a,
	const subgroup *const H2, const integer_t *const x
){

	group_elem_id g;
	// Print the rows corresponding to invertible objects.
	for(g = 0; g < G->order; ++g){
		size_t i;
		for(i = 0; i < H1->quotient.num_cosets; ++i){
			printf("%u ", H1->quotient.coset_ids[g] == i);
		}
		for(i = 0; i < H2->quotient.num_cosets; ++i){
			printf("%u ", 0);
		}
		printf("\n");
	}
	// Print the rows corresponding to non-invertible objects.
	for(g = 0; g < G->order; ++g){
		size_t i;
		for(i = 0; i < H1->quotient.num_cosets; ++i){
			printf("%u ", a[H1->quotient.coset_ids[
				G->multiply[G->inverse[H1->quotient.cosets[i][0]]][g]
			]]);
		}
		for(i = 0; i < H2->quotient.num_cosets; ++i){
			printf("%u ", x[H2->quotient.coset_ids[
				G->multiply[G->inverse[g]][H2->quotient.cosets[i][0]]
			]]);
		}
		printf("\n");
	}
}

// We often need to compute whether an object of the form gX is an
// inner automorphism of some hX (that is, gX = khXk^{-1} for some
// k \in G). For many categories (especially crossed products), we
// can have k\psi(k^{-1}) = 1, so we need not check every element.
static size_t compute_inner_factors(
	const hi_ring *const R, group_elem_id **inner_factors
){

	size_t num_inner_factors = 0;
	group_elem_id k;

	// We start by assuming the worst case scenario (every element
	// gives a non-trivial inner factor), then reallocate to a
	// smaller size once we know how many we actually have.
	group_elem_id *temp_factors = malloc(R->G.order*sizeof(*temp_factors));
	if(temp_factors == NULL){
		puts("Failed to allocate memory for inner factors!");
		*inner_factors = NULL;
		return -1;
	}

	// The inner factors are essentially the elements
	// that aren't fixed by the automorphism \psi.
	for(k = 0; k < R->G.order; ++k){
		const group_elem_id psi_k = R->psi[k];
		if(psi_k != k){
			temp_factors[num_inner_factors] =
				R->G.multiply[k][R->G.inverse[psi_k]];
			++num_inner_factors;
		}
	}

	// Reallocate our array now that we know how big it is.
	if(num_inner_factors > 0){
		*inner_factors = realloc(
			temp_factors, num_inner_factors*sizeof(*temp_factors)
		);
		if(*inner_factors == NULL){
			puts("Failed to reallocate memory for inner factors!");
			free(temp_factors);
			return -1;
		}
	}else{
		free(temp_factors);
		*inner_factors = NULL;
	}

	return num_inner_factors;
}

// Return whether the permutation a is a valid set of
// multiplicities for the algebra object candidate
//     A = H \oplus \bigoplus_{g \in G}{gX^{\oplus a[gH]}}.
// There are two conditions that this function checks.
//     1. If a is an inner automorphism by H of a permutation
//        we've already visited, it must be a simple object
//        in a module category that we've already seen before.
//
//     2. We require a[\psi(g^{-1}h)H] = a[gH] for all g \in G
//        and h \in H. This comes from module associativity.
static int valid_algebra_multiplicities(
	const hi_ring *const R,
	const group_elem_id *const inner_factors, const size_t num_inner_factors,
	const subgroup *const H, const integer_t *const a
){

	group_elem_id g;

	size_t i;
	for(i = 0; i < num_inner_factors; ++i){
		const group_elem_id inner_i = inner_factors[i];
		size_t j;
		for(j = 0; j < H->quotient.num_cosets; ++j){
			// It doesn't matter which element of the coset we choose.
			// a[g\psi(g^{-1})hH]
			const integer_t a_aut = a[H->quotient.coset_ids[
				R->G.multiply[inner_i][H->quotient.cosets[j][0]]
			]];
			// If the first difference between a and its inner
			// automorphism by some element results in smaller
			// coordinate, then the automorphism must precede a
			// lexicographically. Thus, we've already seen it.
			if(a_aut < a[j]){
				return 0;
			}else if(a_aut > a[j]){
				break;
			}
		}
	}

	// Check that a[\psi(g^{-1}h)H] = a[gH] for all g \in G and h \in H.
	for(g = 0; g < R->G.order; ++g){
		const group_elem_id ginv = R->G.inverse[g];
		const integer_t ag = a[H->quotient.coset_ids[g]];
		size_t j;
		for(j = 0; j < H->order; ++j){
			// \psi(g^{-1}h)
			const group_elem_id psi_ginv_h = R->psi[
				R->G.multiply[ginv][H->elem_ids[j]]
			];
			if(a[H->quotient.coset_ids[psi_ginv_h]] != ag){
				return 0;
			}
		}
	}

	return 1;
}

// Return whether the permutation x_1 is a valid set of
// multiplicities for the module object associated to
//     A_1 = H_1 \oplus \bigoplus_{g \in G}{gX^{\oplus a_1[gH_1]}}.
// There are two conditions that this function checks.
//     1. If the corresponding module object with internal end
//        A_1 is an inner automorphism by G of an object we've
//        already found, we won't find anything new from it.
//
//     2. We require x_1[g\psi(h)H_2] = x_1[gH_2] for all g \in G
//        and h \in H_1. This comes from module associativity.
static int valid_module_multiplicities(
	const hi_ring *const R,
	const subgroup *const H1, const subgroup *const H2,
	const integer_t *const a1, const integer_t *const x1
){

	group_elem_id g;

	size_t i;
	for(i = 0; i < H1->quotient.num_cosets; ++i){
		size_t j;
		for(j = 0; j < H2->quotient.num_cosets; ++j){
			// It doesn't matter which element of the coset we choose.
			// \psi(g^{-1}H_1)
			const group_elem_id psi_ginv =
				R->psi[R->G.inverse[H1->quotient.cosets[i][0]]];
			// x1[\psi(g^{-1}H_1)hH_2]
			const integer_t x1_aut = x1[H2->quotient.coset_ids[
				R->G.multiply[psi_ginv][H2->quotient.cosets[j][0]]
			]];
			// If the set x from the fusion matrix of M_{gH} (x1_aut)
			// is a permutation that precedes the current set (x_1) in
			// lexicographic order, we may have seen this module.
			if(x1_aut < x1[j]){
				size_t k;
				for(k = 0; k < H1->quotient.num_cosets; ++k){
					// a1[\psi(g^{-1}H_1)hgH_1]
					const integer_t a1_aut = a1[H1->quotient.coset_ids[
						R->G.multiply[
							R->G.multiply[psi_ginv][H1->quotient.cosets[k][0]]
						][H1->quotient.cosets[i][0]]
					]];
					// If the resulting multiplicities a for M_{gH}
					// are different than what they should be, we
					// can't have seen this module before.
					if(a1_aut != a1[k]){
						break;
					}
				}
				// Otherwise, if the rest of the fusion matrix is a
				// match, we can safely ignore the new permutation x_1.
				if(k == H1->quotient.num_cosets){
					return 0;
				}

				break;
			}else if(x1_aut > x1[j]){
				break;
			}
		}
	}

	// Check that x_1[g\psi(h)H_2] = x_1[gH_2] for all g \in G and h \in H_1.
	for(g = 0; g < R->G.order; ++g){
		const integer_t x1g = x1[H2->quotient.coset_ids[g]];
		size_t j;
		for(j = 0; j < H1->order; ++j){
			// g\psi(h)
			const group_elem_id g_psi_h =
				R->G.multiply[g][R->psi[H1->elem_ids[j]]];
			if(x1[H2->quotient.coset_ids[g_psi_h]] != x1g){
				return 0;
			}
		}
	}

	return 1;
}

// For the corresponding module category to be associative, we need
//     (X \otimes X) \otimes M_{gH} \cong X \otimes (X \otimes M_{gH}.
// This gives us two equations for each orbit of the module category
// coming from equating the multiplicities of objects in the two orbits.
static int valid_module_associativity(
	const hi_ring *const R,
	const subgroup *const H1, const subgroup *const H2,
	const integer_t *const a1, const integer_t *const a2,
	const integer_t *const x1, const integer_t *const x2
){

	size_t i;
	for(i = 0; i < H1->quotient.num_cosets; ++i){
		size_t j;
		const group_elem_id gH1 = H1->quotient.cosets[i][0];
		const group_elem_id ginvH1 = R->G.inverse[H1->quotient.cosets[i][0]];
		integer_t lhs1 = (i == 0) ? 1 : 0;
		integer_t lhs2 = 0;
		integer_t rhs1 = 0;
		integer_t rhs2 = 0;

		for(j = 0; j < R->Q->order; ++j){
			const group_elem_id q = R->Q->elem_ids[j];
			// a1[g^{-1}qH_1]
			lhs1 += a1[H1->quotient.coset_ids[
				R->G.multiply[ginvH1][q]
			]];
			// x2[qgH_1]
			lhs2 += x2[H1->quotient.coset_ids[
				R->G.multiply[q][gH1]
			]];
		}
		for(j = 0; j < H1->quotient.num_cosets; ++j){
			const group_elem_id hH1 = H1->quotient.cosets[j][0];
			// a1[hH_1]a1[g^{-1}H_1\psi(h^{-1}H_1)]
			rhs1 += a1[j]*a1[H1->quotient.coset_ids[
				R->G.multiply[ginvH1][R->psi[R->G.inverse[hH1]]]
			]];
			// x2[hH_1]a1[g^{-1}H_1\psi(hH_1)]
			rhs2 += x2[j]*a1[H1->quotient.coset_ids[
				R->G.multiply[ginvH1][R->psi[hH1]]
			]];
		}
		for(j = 0; j < H2->quotient.num_cosets; ++j){
			const group_elem_id hH2 = H2->quotient.cosets[j][0];
			// x1[hH_2]x2[\psi(h^{-1}H_2)gH_1]
			rhs1 += x1[j]*x2[H1->quotient.coset_ids[
				R->G.multiply[R->psi[R->G.inverse[hH2]]][gH1]
			]];
			// a2[hH_2]x2[\psi(hH_2)gH_1]
			rhs2 += a2[j]*x2[H1->quotient.coset_ids[
				R->G.multiply[R->psi[hH2]][gH1]
			]];
		}

		if(lhs1 != rhs1 || lhs2 != rhs2){
			return 0;
		}
	}

	return 1;
}

// Once we've chosen candidates for x_1 and x_2, we have candidate
// fusion matrices F_{M_1} and F_{M_2} for the two module objects.
// We need to check that these agree with the fusion matrices for
// our two algebra objects in the sense that
//     F^{M_1}(F^{M_1})^T = F^{A_1},
//     F^{M_2}(F^{M_2})^T = F^{A_2}.
// This can be reduced two checking two equations. Because we need
// to loop over H_1 and H_2 for both, we check them simultaneously.
static int valid_module_rigidity(
	const hi_ring *const R,
	const subgroup *const H1, const subgroup *const H2,
	const integer_t *const a1, const integer_t *const a2,
	const integer_t *const x1, const integer_t *const x2
){

	group_elem_id g;
	for(g = 0; g < R->G.order; ++g){
		size_t j;
		// \psi(g^{-1})
		const group_elem_id psi_ginv = R->psi[R->G.inverse[g]];
		integer_t lhs1 = (H1->quotient.coset_ids[psi_ginv] == 0) ? 1 : 0;
		integer_t lhs2 = (H2->quotient.coset_ids[psi_ginv] == 0) ? 1 : 0;
		integer_t rhs1 = 0;
		integer_t rhs2 = 0;

		for(j = 0; j < R->Q->order; ++j){
			const group_elem_id psi_ginv_q =
				R->G.multiply[psi_ginv][R->psi[R->Q->elem_ids[j]]];
			// a1[\psi(g^{-1}q)H_1]
			lhs1 += a1[H1->quotient.coset_ids[psi_ginv_q]];
			// a2[\psi(g^{-1}q)H_2]
			lhs2 += a2[H2->quotient.coset_ids[psi_ginv_q]];
		}
		for(j = 0; j < H1->quotient.num_cosets; ++j){
			const size_t hgH1_coset = H1->quotient.coset_ids[
				R->G.multiply[H1->quotient.cosets[j][0]][g]
			];
			// a1[hgH1]a1[hH1]
			rhs1 += a1[hgH1_coset]*a1[j];
			// x2[hgH1]x2[hH1]
			rhs2 += x2[hgH1_coset]*x2[j];
		}
		for(j = 0; j < H2->quotient.num_cosets; ++j){
			const size_t hgH2_coset = H2->quotient.coset_ids[
				R->G.multiply[H2->quotient.cosets[j][0]][g]
			];
			// x1[hgH2]x1[hH2]
			rhs1 += x1[hgH2_coset]*x1[j];
			// a2[hgH2]a2[hH2]
			rhs2 += a2[hgH2_coset]*a2[j];
		}

		if(lhs1 != rhs1 || lhs2 != rhs2){
			return 0;
		}
	}

	return 1;
}

// Try and build modules that give a Morita equivalence between the algebras
//     A_1 = H_1 \oplus \bigoplus_{g \in G}{gX^{\oplus a_1[gH_1]}},
//     A_2 = H_2 \oplus \bigoplus_{g \in G}{gX^{\oplus a_2[gH_2]}}.
static void construct_modules(
	hi_ring *const R,
	integer_t *const scratchpad, const integer_t a1_sum,
	const subgroup *const H1, const integer_t *const a1,
	const subgroup *const H2, const integer_t *const a2
){

	integer_t x1_sum;
	integer_t x2_sum;
	partition p;
	permutation x1;

	// The algebra fusion matrix relies on integers x_1[gH_2] satisfying
	//     |H_2|(\sum_{gH_2 \in G/H_2}{x_1[gH_2]})^2 = |H_1|(1 + \sum_{gH_1 \in G/H_1}{a_1[gH_1]}\sum_{gH_2 \in G/H_2}{a_2[gH_2]}).
	// We have a similar set of constants x_2[gH_1] for our other candidate.
	// Start by computing the right-hand sides of the sums above.
	integer_t x1_square_sum = H1->order*(1 + a1_sum*(R->Q->order - a1_sum));
	integer_t x2_square_sum = H2->order*(1 + a1_sum*(R->Q->order - a1_sum));
	// Make sure they're evenly divisible by the subgroup orders.
	if((x1_square_sum % H2->order) != 0 || (x2_square_sum % H1->order) != 0){
		return;
	}
	x1_square_sum /= H2->order;
	x2_square_sum /= H1->order;

	x1_sum = sqrt(x1_square_sum);
	x2_sum = sqrt(x2_square_sum);
	// They need to be squares, too!
	if(x1_sum*x1_sum != x1_square_sum || x2_sum*x2_sum != x2_square_sum){
		return;
	}

	// Set up the pointers for our partition and permutation.
	p.parts   = scratchpad;
	x1.values = p.parts + R->G.order;

	// Now that we know what the elements of
	// x_1 sum to, we can look for partitions.
	// This builds a complete fusion matrix for A.
	partition_start(&p, x1_sum, H2->quotient.num_cosets);
	do {
		permutation_start(&x1, H2->quotient.num_cosets, p.parts, p.length);
		do {
			integer_t *x2;
			size_t i;

			if(!valid_module_multiplicities(R, H1, H2, a1, x1.values)){
				continue;
			}

			x2 = x1.values + R->G.order;
			// Compute the x2[gH_1] using x2[gH_1] = x1[\psi(g^{-1})H_2].
			// One of the conditions we check in the function
			// "valid_module_multiplicities()" ensures this makes sense.
			for(i = 0; i < H1->quotient.num_cosets; ++i){
				x2[i] = x1.values[H2->quotient.coset_ids[
					// \psi(g^{-1})
					R->psi[H1->quotient.cosets[i][0]]
				]];
			}

			if(
				valid_module_associativity(R, H1, H2, a1, a2, x1.values, x2) &&
				valid_module_associativity(R, H2, H1, a2, a1, x2, x1.values) &&
				valid_module_rigidity(R, H1, H2, a1, a2, x1.values, x2)
			){

				hi_module M;
				M.H1 = H1;
				M.a1 = a1;
				M.x1 = x1.values;
				M.H2 = H2;
				M.a2 = a2;
				M.x2 = x2;
				hi_ring_print_module(R, &M);
			}
		} while(permutation_next(&x1));
	} while(partition_next(&p));
}