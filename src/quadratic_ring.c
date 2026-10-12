#include "quadratic_ring.h"


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "partition.h"
#include "permutation.h"

#include "module_checker.h"


static int compute_psi_subgroups(quadratic_ring *const R);

static void print_module_matrix(
	const group *const G,
	const subgroup *const H1, const integer_t *const a,
	const subgroup *const H2, const integer_t *const x
);


int quadratic_ring_load(quadratic_ring *const R, const char *const path){
	long unsigned int scan_val;

	FILE *const fp = fopen(path, "r");
	if(fp == NULL){
		printf("Failed to open file \"%s\"!\n", path);
		return 0;
	}

	// Initialize the structure. This is mostly
	// just so that if anything fails, we can
	// use our normal deletion function.
	R->psi         = NULL;
	R->m           = NULL;
	R->modules     = NULL;
	R->psi_ids     = NULL;
	R->num_modules = 0;


	// Load the group of invertibles and all of its subgroups.
	if(!group_load(&R->G, fp)){
		printf("Failed to load group from file \"%s\"!\n", path);
		fclose(fp);
		return 0;
	}
	// Load the index of the subgroup Q.
	fscanf(fp, "%lu", &scan_val);
	if(scan_val < R->G.num_subgroups){
		R->Q_id = scan_val;
		R->Q = &R->G.subgroups[scan_val];
	}else{
		printf(
			"Given Q index %lu, but "
			"only has %lu subgroups!\n",
			scan_val, (long unsigned int)R->G.num_subgroups
		);
		return 0;
	}


	R->psi = malloc(R->G.order*sizeof(*R->psi));
	if(R->psi == NULL){
		puts("Failed to allocate memory for \\psi!");
		fclose(fp);
		return 0;
	}
	{
		group_elem_id g;
		// Load the group automorphism \psi.
		for(g = 0; g < R->G.order; ++g){
			fscanf(fp, "%lu", &scan_val);
			if(scan_val < R->G.order){
				R->psi[g] = scan_val;
			}else{
				printf(
					"Given \\psi group element ID %lu, "
					"but G only has order %lu!\n",
					scan_val, (long unsigned int)R->G.order
				);
				return 0;
			}
		}
		// Check that \psi is an bijection of G/Q.
		for(g = 0; g < R->G.order; ++g){
			const size_t g_coset = R->Q->quotient.coset_ids[g];
			const size_t psi_g_coset =
				R->Q->quotient.coset_ids[R->psi[g]];

			group_elem_id h;
			for(h = g + 1; h < R->G.order; ++h){
				const size_t h_coset = R->Q->quotient.coset_ids[h];
				const size_t psi_h_coset =
					R->Q->quotient.coset_ids[R->psi[h]];

				if(g_coset == h_coset && psi_g_coset != psi_h_coset){
					printf(
						"The provided \\psi is not a function on G/Q! "
						"Please check \\psi(%lu) and \\psi(%lu).\n",
						(long unsigned int)g, (long unsigned int)h
					);
					return 0;
				}
				if(g_coset != h_coset && psi_g_coset == psi_h_coset){
					printf(
						"The provided \\psi is not injective! "
						"Please check \\psi(%lu) and \\psi(%lu).\n",
						(long unsigned int)g, (long unsigned int)h
					);
					return 0;
				}
			}
		}
		// Check that \psi is a group homomorphism.
		// We do this after checking if it's a bijection
		// so the error reporting is easier to understand.
		for(g = 0; g < R->G.order; ++g){
			const group_elem_id psi_g = R->psi[g];

			group_elem_id h;
			for(h = g; h < R->G.order; ++h){
				const size_t psi_gh_coset = R->Q->quotient.coset_ids[
					R->psi[R->G.multiply[g][h]]
				];
				const size_t psi_g_psi_h_coset = R->Q->quotient.coset_ids[
					R->G.multiply[psi_g][R->psi[h]]
				];

				if(psi_gh_coset != psi_g_psi_h_coset){
					printf(
						"The function \\psi is not a group homomorphism! "
						"Please check \\psi(%lu) and \\psi(%lu).\n",
						(long unsigned int)g, (long unsigned int)h
					);
					return 0;
				}
			}
		}
	}

	// Load an element in r(Q).
	fscanf(fp, "%lu", &scan_val);
	if(scan_val < R->G.order){
		R->r = scan_val;
		if(
			R->Q->quotient.coset_ids[R->psi[R->r]] !=
			R->Q->quotient.coset_ids[R->r]
		){

			printf(
				"The element %lu of r(Q) given "
				"is not fixed by \\psi!\n",
				scan_val
			);
			return 0;
		}
	}else{
		printf(
			"Given r(Q) group element ID %lu, "
			"but G only has order %lu!\n",
			scan_val, (long unsigned int)R->G.order
		);
		return 0;
	}

	R->m = malloc(R->G.order*sizeof(*R->m));
	if(R->m == NULL){
		puts("Failed to allocate memory for m!");
		fclose(fp);
		return 0;
	}
	{
		group_elem_id g;
		// Load the multiplicities m.
		// We also compute their sum here.
		R->m_sum = 0;
		for(g = 0; g < R->G.order; ++g){
			fscanf(fp, "%lu", &scan_val);
			if(scan_val >= 0){
				R->m[g] = scan_val;
			}else{
				printf(
					"Invalid multiplicity m_%lu = %lu < 0!\n",
					(long unsigned int)g, scan_val
				);
				return 0;
			}
			R->m_sum += R->m[g];
		}
		// Check that m is a function of G/Q.
		for(g = 0; g < R->G.order; ++g){
			const size_t g_coset = R->Q->quotient.coset_ids[g];
			const integer_t m_g = R->m[g];

			group_elem_id h;
			for(h = g + 1; h < R->G.order; ++h){
				const size_t h_coset = R->Q->quotient.coset_ids[h];

				if(g_coset == h_coset && m_g != R->m[h]){
					printf(
						"The provided m is not a function on G/Q! "
						"Please check m(%lu) and m(%lu).\n",
						(long unsigned int)g, (long unsigned int)h
					);
					return 0;
				}
			}
		}
		// Check that m_{g^{-1}r(Q)} = m_{gQ} and
		// m_{g\psi(hQ)} = m_{ghQ}. We do this after
		// checking if it's a function so the error
		// reporting is easier to understand.
		for(g = 0; g < R->G.order; ++g){
			group_elem_id h;

			const integer_t m_ginv_r = R->m[
				R->G.multiply[R->G.inverse[g]][R->r]
			];
			if(m_ginv_r != R->m[g]){
				printf(
					"The function m doesn't satisfy "
					"m_{g^{-1}r(Q)} = m_{gQ}! Please "
					"check g = %lu.\n",
					(long unsigned int)g
				);
				return 0;
			}

			// Note that if h < g, then we can write
			//     m_{g\psi(hQ)} = m_{\psi(g\psi(hQ))}
			//                   = m_{h\psi(gQ)}
			//                   = m_{ghQ}.
			// Thus, we only need to check h >= g.
			for(h = g; h < R->G.order; ++h){
				const integer_t m_g_psi_h =
					R->m[R->G.multiply[g][R->psi[h]]];
				const integer_t m_gh = R->m[R->G.multiply[g][h]];

				if(m_g_psi_h != m_gh){
					printf(
						"The function m doesn't satisfy "
						"m_{g\\psi(hQ)} = m_{ghQ}! Please "
						"check g = %lu and h = %lu.\n",
						(long unsigned int)g, (long unsigned int)h
					);
					return 0;
				}
			}
		}
		// When we summed up the values of m, we had some repeats
		// since we store a value for each element of G rather
		// than of G/Q. We fix this by simply dividing by |Q|.
		R->m_sum /= R->Q->order;
	}

	// Compute d = FPdim(X_Q), or store d = 0 if it's irrational.
	// Note that FPdim(X_Q) is an integer if it's not irrational.
	{
		const integer_t disc = R->m_sum*R->m_sum + 4*R->Q->order;
		const integer_t sqrt_disc = (integer_t)sqrt(disc);
		if(sqrt_disc*sqrt_disc == disc){
			R->d = (R->m_sum + sqrt_disc)/2;
			if(2*R->d != R->m_sum + sqrt_disc){
				puts(
					"The given quadratic fusion ring is not "
					"categorifiable, as FPdim(X_Q) is neither "
					"irrational nor an integer!"
				);
				return 0;
			}
		}else{
			R->d = 0;
		}
	}

	fclose(fp);
	return compute_psi_subgroups(R);
}

// We often need to know when two non-invertible objects are
// inner automorphisms of each other. For gQ, g'Q \in G/Q,
//     gX_{g'Q}g^{-1} \cong X_{gg'\psi(g^{-1}Q)}.
// Thus, X_{g'Q} is an inner automorphism of X_{g''Q}
// if g''Q = g^{-1}\psi(gQ)g'Q for some gQ \in G/Q.
// We precompute the non-trivial factors g^{-1}\psi(gQ)
// and return how many there are.
size_t quadratic_ring_compute_inner_factors(
	const quadratic_ring *const R,
	group_elem_id *const inner_factors
){

	size_t num_inner_factors = 0;
	group_elem_id g_coset;

	// The inner factors are essentially the elements
	// that aren't fixed by the automorphism \psi.
	for(g_coset = 0; g_coset < R->Q->quotient.num_cosets; ++g_coset){
		const group_elem_id g = R->Q->quotient.cosets[g_coset][0];
		const group_elem_id psi_g = R->psi[g];
		if(R->Q->quotient.coset_ids[psi_g] != g_coset){
			inner_factors[num_inner_factors] =
				R->G.multiply[R->G.inverse[g]][psi_g];
			++num_inner_factors;
		}
	}

	return num_inner_factors;
}

void quadratic_ring_delete(quadratic_ring *const R){
	free(R->modules);
	free(R->psi_ids);
	free(R->m);
	free(R->psi);
	group_delete(&R->G);
}


int quadratic_module_init(
	quadratic_module *const M, const quadratic_ring *const R,
	const size_t *const H_ids, integer_t ***const a, const size_t num_orbits
){

	size_t i;
	const size_t a_row_bytes = num_orbits*sizeof(*a);
	const size_t a_col_bytes = num_orbits*sizeof(**a);
	size_t total_bytes =
		num_orbits*sizeof(*H_ids) +
		a_row_bytes + num_orbits*a_col_bytes;

	// Figure out how much memory we need for the module.
	for(i = 0; i < num_orbits; ++i){
		const size_t Hi_id = H_ids[i];
		size_t j;
		for(j = 0; j < num_orbits; ++j){
			const subgroup *const psi_Hij =
				&R->G.subgroups[R->psi_ids[Hi_id][H_ids[j]]];
			total_bytes += psi_Hij->quotient.num_cosets*sizeof(***a);
		}
	}


	M->H_ids = malloc(total_bytes);
	if(M->H_ids == NULL){
		printf("Failed to allocate memory for module!\n");
		return 0;
	}
	M->a = (integer_t ***)&M->H_ids[num_orbits];
	M->num_orbits = num_orbits;

	// Copy the subgroup IDs over.
	memcpy(M->H_ids, H_ids, num_orbits*sizeof(*H_ids));

	{
		integer_t **cur_col = (integer_t **)&M->a[num_orbits];
		integer_t *cur_elem = (integer_t *)&cur_col[num_orbits*num_orbits];
		// Now copy the multiplicities over.
		for(i = 0; i < num_orbits; ++i){
			const size_t Hi_id = H_ids[i];
			size_t j;

			M->a[i] = cur_col;
			for(j = 0; j < num_orbits; ++j){
				const subgroup *const psi_Hij =
					&R->G.subgroups[R->psi_ids[Hi_id][H_ids[j]]];

				M->a[i][j] = cur_elem;
				memcpy(
					M->a[i][j], a[i][j],
					psi_Hij->quotient.num_cosets*sizeof(***a)
				);
				cur_elem += psi_Hij->quotient.num_cosets;
			}
			cur_col += num_orbits;
		}
	}


	return 1;
}

void quadratic_module_print(
	const quadratic_module *const M,
	const quadratic_ring *const R
){

	#if 0
	print_module_matrix(&R->G, M->H1, M->a1, M->H2, M->x1);
	printf("\n");
	print_module_matrix(&R->G, M->H2, M->a2, M->H1, M->x2);
	printf("\n\n");
	#endif
}

void quadratic_module_delete(quadratic_module *const M){
	free(M->H_ids);
}


// For each pair of subgroups H_i and H_j,
// compute the ID of \psi(H_i Q)H_j Q.
static int compute_psi_subgroups(quadratic_ring *const R){
	size_t i;
	const size_t total_bytes = R->G.num_subgroups*(
		sizeof(*R->psi_ids) + R->G.num_subgroups*sizeof(**R->psi_ids)
	);
	group_elem_id *psi_elems;
	size_t *cur_offset;

	// Allocate memory for the subgroup IDs and psi_elems.
	// We reallocate later once we're done with the latter.
	R->psi_ids = malloc(
		total_bytes + R->G.order*R->Q->order*sizeof(*psi_elems)
	);
	if(R->psi_ids == NULL){
		puts("Unable to allocate memory for \\psi subgroups!");
		return 0;
	}
	psi_elems = (group_elem_id *)((char *)R->psi_ids + total_bytes);

	cur_offset = (size_t *)&R->psi_ids[R->G.num_subgroups];
	for(i = 0; i < R->G.num_subgroups; ++i){
		const subgroup *const H = &R->G.subgroups[i];
		size_t j;

		// Find every element in the set \psi(HQ). This may
		// generate duplicate elements, but that's okay!
		for(j = 0; j < H->order; ++j){
			const size_t psi_H_coset =
				R->Q->quotient.coset_ids[R->psi[H->elem_ids[j]]];
			size_t k;
			for(k = 0; k < R->Q->order; ++k){
				psi_elems[j*R->Q->order + k] =
					R->Q->quotient.cosets[psi_H_coset][k];
			}
		}

		R->psi_ids[i] = cur_offset;
		// Store the ID of the subgroup with these elements.
		// Note that \psi(H_i Q)Q = \psi(H_i Q), so we can
		// store it in the column corresponding to Q for now.
		R->psi_ids[i][R->Q_id] = subgroup_find_id(
			&R->G, psi_elems, H->order*R->Q->order
		);
		cur_offset += R->G.num_subgroups;
	}

	R->psi_ids = realloc(R->psi_ids, total_bytes);
	if(R->psi_ids == NULL){
		puts("Unable to allocate memory for \\psi subgroups!");
		return 0;
	}

	// Now that we have \psi(H_i Q), we can compute \psi(H_i Q)H_j Q.
	for(i = 0; i < R->G.num_subgroups; ++i){
		const size_t psi_H_id = R->psi_ids[i][R->Q_id];
		size_t j;
		for(j = 0; j < R->G.num_subgroups; ++j){
			R->psi_ids[i][j] = R->G.product_ids[psi_H_id][j];
		}
	}

	return 1;
}


static void print_module_matrix(
	const group *const G,
	const subgroup *const H1, const integer_t *const a,
	const subgroup *const H2, const integer_t *const x
){

	#if 0
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
	#endif
}