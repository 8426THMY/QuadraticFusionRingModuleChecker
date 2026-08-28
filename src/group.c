#include "group.h"


#include <stdlib.h>
#include <string.h>


static int init_group_lookups(group *const G);


// Let G = Z/p_1^{e_1}Z \times ... \times Z/p_n^{e_n}Z.
// The index of a group element g = (g_1, ..., g_n) is
//     g_id = \sum_{i=1}^{n}{g_i\prod_{j=1}^{i-1}{p_j^{e_j}}}.
group_elem_id group_elem_get_id(
	const group *const G,
	const group_elem g
){

	group_elem_id g_id = 0;
	size_t offset = 1;

	size_t i;
	for(i = 0; i < G->num_pgroups; ++i){
		const pgroup *const Gp = &G->pgroups[i];
		size_t j;
		for(j = 0; j < Gp->num_powers; ++j){
			g_id += g[i][j]*offset;
			offset *= Gp->powers[j];
		}
	}

	return g_id;
}

// Return the index of g + h.
group_elem_id group_elem_multiply(
	const group *const G,
	const group_elem g, const group_elem h
){

	group_elem_id sum_id = 0;
	size_t offset = 1;

	size_t i;
	for(i = 0; i < G->num_pgroups; ++i){
		const pgroup *const Gp = &G->pgroups[i];
		size_t j;
		for(j = 0; j < Gp->num_powers; ++j){
			sum_id += ((g[i][j] + h[i][j]) % Gp->powers[j])*offset;
			offset *= Gp->powers[j];
		}
	}

	return sum_id;
}

// Return the index of g^{-1}.
group_elem_id group_elem_invert(
	const group *const G,
	const group_elem g
){

	group_elem_id g_id = 0;
	size_t offset = 1;

	size_t i;
	for(i = 0; i < G->num_pgroups; ++i){
		const pgroup *const Gp = &G->pgroups[i];
		size_t j;
		for(j = 0; j < Gp->num_powers; ++j){
			g_id += g[i][j]*offset;
			offset *= Gp->powers[j];
		}
	}

	return g_id;
}

void group_elem_print(const group *const G, const group_elem g){
	size_t i;

	printf("(");
	for(i = 0; i < G->num_pgroups; ++i){
		const pgroup *const Gp = &G->pgroups[i];
		size_t j;
		for(j = 0; j < Gp->num_powers - 1; ++j){
			printf("%i, ", g[i][j]);
		}
		if(i == G->num_pgroups - 1){
			printf("%i", g[i][j]);
		}else{
			printf("%i, ", g[i][j]);
		}
	}
	printf(")");
}


int subgroup_init_quotient(const group *const G, subgroup *const H){
	size_t i;
	group_elem_id g = 0;
	quotient_group *const G_H = &H->quotient;
	char *start_bytes;
	const size_t coset_bytes = H->order*sizeof(**G_H->cosets);


	// Allocate memory for a temporary array of element IDs. This
	// lets us keep track of which elements we've sorted into cosets.
	group_elem_id *const elem_sorted = calloc(G->order, sizeof(*elem_sorted));
	if(elem_sorted == NULL){
		puts("Failed to allocate memory for temporary subgroup array!");
		return 0;
	}


	G_H->num_cosets = G->order/H->order;
	// Allocate memory for the array of cosets
	// and the array of element coset IDs.
	G_H->cosets = malloc(
		G_H->num_cosets*(sizeof(*G_H->cosets) + coset_bytes) +
		G->order*sizeof(*G_H->coset_ids)
	);
	if(G_H->cosets == NULL){
		puts("Failed to allocate memory for subgroup data!");
		free(elem_sorted);
		return 0;
	}
	// Since we allocated memory for the element coset IDs
	// above, we'll need to adjust the pointer before we use it.
	G_H->coset_ids = (size_t *)((char *)G_H->cosets +
		G_H->num_cosets*(sizeof(*G_H->cosets) + coset_bytes)
	);


	start_bytes = (char *)&G_H->cosets[G_H->num_cosets];
	for(i = 0; i < G_H->num_cosets; ++i){
		size_t j;

		// Set up the pointer to the coset's element array.
		G_H->cosets[i] = (group_elem_id *)(start_bytes + i*coset_bytes);
		// Find the next unsorted group element.
		while(elem_sorted[g]){
			++g;
		}
		// Multiply g with every element of H to fill out the coset.
		for(j = 0; j < H->order; ++j){
			const group_elem_id gh = G->multiply[g][H->elem_ids[j]];
			G_H->cosets[i][j] = gh;
			G_H->coset_ids[gh] = i;
			elem_sorted[gh] = 1;

		}
	}


	return 1;
}

#if 0
void subgroup_init_zero(
	const group *const G, subgroup *const H,
	const size_t num_gens
){
	size_t i;
	// A group element is a matrix of integers. We need to
	// allocate memory for the row pointers and their entries.
	const size_t elem_bytes = (
		G->num_pgroups*sizeof(pgroup_elem) +
		G->num_factors*sizeof(integer_t)
	);
	char *start_bytes;

	H->gens = malloc(num_gens*(sizeof(group_elem) + elem_bytes));
	H->num_gens = num_gens;

	start_bytes = (char *)(H->gens + num_gens);
	// Initialize each generator to zero.
	for(i = 0; i < num_gens; ++i){
		integer_t *pgroup_offset;
		size_t j;

		// Initialize the pgroup_elem array
		// pointers for each group_elem.
		H->gens[i] = (pgroup_elem *)(start_bytes + i*elem_bytes);

		pgroup_offset = (integer_t *)(H->gens[i] + G->num_pgroups);
		// Initialize each group element to identity.
		memset(pgroup_offset, 0, G->num_factors*sizeof(integer_t));

		for(j = 0; j < G->num_pgroups; ++j){
			// Initialize the integer array
			// pointers for each pgroup_elem.
			H->gens[i][j] = pgroup_offset;
			pgroup_offset += G->pgroups[j].num_powers;
		}
	}
}

// Initialize H = H1 \oplus H2. We assume
// that H1 and H2 intersect trivially!
void subgroup_init_direct_sum(
	const group *const G, subgroup *const H,
	const subgroup *const H1, const subgroup *const H2
){

	size_t i;

	subgroup_init_zero(G, H, H1->num_gens + H2->num_gens);

	// Copy the generators from H1.
	for(i = 0; i < H1->num_gens; ++i){
		memcpy(H->gens[i][0], H1->gens[i][0], G->num_factors*sizeof(integer_t));
	}
	// Copy the generators from H2.
	for(i = 0; i < H2->num_gens; ++i){
		memcpy(H->gens[H1->num_gens + i][0], H2->gens[i][0], G->num_factors*sizeof(integer_t));
	}
}
#endif


// A p-group is specified in the following format:
//     prime num_powers power1 power2 ... powern
// Prime powers should be given in ascending order!
int pgroup_load(pgroup *const Gp, FILE *const fp){
	size_t i;
	long unsigned int scan_val;

	// ISO C doesn't support long long ints with fscanf,
	// so just read to a long unsigned int and copy it over.
	fscanf(fp, "%i %lu", &Gp->prime, &scan_val);
	Gp->num_powers = scan_val;
	// Allocate and initialize the array of powers for this prime.
	Gp->powers = malloc(Gp->num_powers*sizeof(*Gp->powers));
	if(Gp->powers == NULL){
		puts("Failed to allocate memory for pgroup powers!");
		return 0;
	}

	Gp->order = 1;
	for(i = 0; i < Gp->num_powers; ++i){
		fscanf(fp, "%i", &Gp->powers[i]);
		Gp->order *= Gp->powers[i];
	}

	return 1;
}

void pgroup_delete(pgroup *const Gp){
	free(Gp->powers);
}


// A group is specified in the following format:
//     num_pgroups
//     pgroup1
//     pgroup2
//     ...
//     pgroupn
// These p-groups should be p-Sylow subgroups,
// that is, we shouldn't repeat any primes.
int group_load(group *const G, FILE *const fp){
	size_t i;
	long unsigned int scan_val;

	// Initialize the structure. This is mostly
	// just so that if anything fails, we can
	// use our normal deletion function.
	G->pgroups = NULL;
	G->num_pgroups = 0;
	G->num_factors = 0;
	G->elems = NULL;
	G->order = 0;
	G->multiply = NULL;
	G->inverse = NULL;
	G->subgroups = NULL;
	G->num_subgroups = 0;


	// ISO C doesn't support long long ints with fscanf,
	// so just read to a long unsigned int and copy it over.
	fscanf(fp, "%lu", &scan_val);
	G->num_pgroups = scan_val;
	// Allocate the array of factors.
	G->pgroups = malloc(G->num_pgroups*sizeof(*G->pgroups));
	if(G->pgroups == NULL){
		puts("Failed to allocate memory for group factors!");
		G->num_pgroups = 0;
		return 0;
	}

	G->num_factors = 0;
	G->order = 1;
	// Load the p-groups for each prime.
	for(i = 0; i < G->num_pgroups; ++i){
		pgroup *const Gp = &G->pgroups[i];
		if(!pgroup_load(Gp, fp) || Gp->num_powers == 0){
			G->num_pgroups = i;
			return 0;
		}
		G->num_factors += Gp->num_powers;
		G->order *= Gp->order;
	}


	// Allocate memory for each group element and create
	// lookup tables for multiplication and inversion.
	if(!init_group_lookups(G)){
		return 0;
	}


	/// This is temporary! We should eventually find the subgroups
	/// ourselves rather than have them specified in the input file.
	fscanf(fp, "%lu", &scan_val);
	G->num_subgroups = scan_val;
	// Allocate the array of subgroups.
	G->subgroups = malloc(G->num_subgroups*sizeof(*G->subgroups));
	if(G->subgroups == NULL){
		puts("Failed to allocate memory for subgroups!");
		G->num_subgroups = 0;
		return 0;
	}

	// Load each subgroup.
	for(i = 0; i < G->num_subgroups; ++i){
		size_t j;
		subgroup *const H = &G->subgroups[i];

		fscanf(fp, "%lu", &scan_val);
		H->order = scan_val;
		// Allocate the array of subgroup element IDs.
		H->elem_ids = malloc(H->order*sizeof(*H->elem_ids));
		if(H->elem_ids == NULL){
			puts("Failed to allocate memory for subgroup elements!");
			G->num_subgroups = i;
			return 0;
		}

		// Initialize the array of element IDs.
		for(j = 0; j < H->order; ++j){
			fscanf(fp, "%lu", &scan_val);
			H->elem_ids[j] = scan_val;
		}

		// Initialize the quotient group G/H.
		if(!subgroup_init_quotient(G, H)){
			free(H->elem_ids);
			G->num_subgroups = i;
			return 0;
		}
	}


	return 1;
}

#if 0
size_t group_count_subgroups(const group *const G){
	size_t num_subgroups = 0;

	size_t i;
	for(i = 0; i < G->num_pgroups; ++i){
		pgroup *const Gp = &G->pgroups[i];
		size_t j;
		for(j = 0; j < Gp->num_powers; ++j){
			size_t cur_subgroups = 1;

			integer_t cur_power;
			// Count up the number of non-trivial subgroups this factor has.
			for(cur_power = Gp->prime; cur_power < Gp->powers[j]; cur_power *= Gp->prime){
				++cur_subgroups;
			}

			// We need to count all of the subgroups for this factor
			// and all of the subgroups we obtain by taking direct
			// products of them with the ones we've already counted.
			num_subgroups += cur_subgroups + cur_subgroups*num_subgroups;
		}
	}

	// Add one for the trivial subgroup!
	return num_subgroups + 1;
}

#error "These functions miss a whole bunch of subgroups (for example, the one generated by (1, 2) in Z/3Z \times Z/3Z)."
void group_compute_subgroups(group *const G){
	size_t cur_num_subgroups = 0;
	size_t i;

	G->num_subgroups = group_count_subgroups(G);
	G->subgroups = malloc(G->num_subgroups*sizeof(*G->subgroups));

	for(i = 0; i < G->num_pgroups; ++i){
		pgroup *const Gp = &G->pgroups[i];
		size_t j;
		// We start by adding the subgroups of each factor.
		for(j = 0; j < Gp->num_powers; ++j){
			// Make sure we don't add generators that
			// live in the same prime power factor.
			const size_t prev_num_subgroups = cur_num_subgroups;

			integer_t cur_power;
			// Note that we don't consider cur_power = 0, as we
			// would otherwise be adding a copy of the trivial
			// subgroup for every factor! We handle this after
			// determining all of the other subgroups.
			for(cur_power = 1; cur_power < Gp->powers[j]; cur_power *= Gp->prime){
				size_t k;

				subgroup *const cur_subgroup = &G->subgroups[cur_num_subgroups];
				subgroup_init_zero(G, cur_subgroup, 1);
				cur_subgroup->gens[0][i][j] = cur_power;
				++cur_num_subgroups;

				// We also have to add the direct product of this new
				// subgroup with each subgroup from the previous factors.
				for(k = 0; k < prev_num_subgroups; ++k){
					subgroup_init_direct_sum(
						G, &G->subgroups[cur_num_subgroups],
						&G->subgroups[k], cur_subgroup
					);
					++cur_num_subgroups;
				}
			}
		}
	}

	// Add the trivial subgroup at last!
	subgroup_init_zero(G, &G->subgroups[cur_num_subgroups], 1);
}
#endif

void group_delete(group *const G){
	size_t i;
	for(i = 0; i < G->num_pgroups; ++i){
		pgroup_delete(&G->pgroups[i]);
	}
	free(G->pgroups);
	free(G->elems);
}


static int init_group_lookups(group *const G){
	group_elem_id i;

	char *cur_offset;
	// A group element consists of an array of
	// pgroup_elems and an integer for each factor.
	const size_t elem_bytes = (
		G->num_pgroups*sizeof(pgroup_elem) +
		G->num_factors*sizeof(integer_t)
	);
	const size_t row_bytes = G->order*sizeof(**G->multiply);
	G->elems = malloc(G->order*(
		sizeof(group_elem) + elem_bytes +
		sizeof(*G->multiply) + row_bytes +
		sizeof(*G->inverse)
	));
	if(G->elems == NULL){
		puts("Failed to allocate memory for group lookups!");
		return 0;
	}
	// Set the pointers for the multiplication and inverse tables.
	G->multiply = (group_elem_id **)((char *)G->elems +
		G->order*(sizeof(group_elem) + elem_bytes)
	);
	G->inverse = (group_elem_id *)((char *)G->multiply +
		G->order*(sizeof(*G->multiply) + row_bytes)
	);

	cur_offset = (char *)&G->elems[G->order];
	// Initialize each group element.
	for(i = 0; i < G->order; ++i){
		integer_t *pgroup_offset;
		size_t j;
		size_t truncated_id = i;

		// Initialize the pgroup_elem array
		// pointers for each group_elem.
		G->elems[i] = (pgroup_elem *)cur_offset;
		cur_offset += elem_bytes;

		pgroup_offset = (integer_t *)&G->elems[i][G->num_pgroups];
		for(j = 0; j < G->num_pgroups; ++j){
			const pgroup *const Gp = &G->pgroups[j];
			size_t k;

			// Initialize the integer array
			// pointers for each pgroup_elem.
			G->elems[i][j] = pgroup_offset;
			pgroup_offset += Gp->num_powers;

			// Let G = Z/p_1^{e_1}Z \times ... \times Z/p_n^{e_n}Z.
			// The index of a group element g = (g_1, ..., g_n) is
			//     g_id = \sum_{i=1}^{n}{g_i\prod_{j=1}^{i-1}{p_j^{e_j}}}.
			// We can reconstruct the group element with this index using
			//     g_k = ((g_id - \sum_{i=1}^{k-1}{g_i\prod_{j=1}^{i-1}{p_j^{e_j}}})/\prod_{j=1}^{k-1}{p_j^{e_j}}) % p_k^{e_k}.
			for(k = 0; k < Gp->num_powers; ++k){
				G->elems[i][j][k] = truncated_id % Gp->powers[k];
				truncated_id = (truncated_id - G->elems[i][j][k])/Gp->powers[k];
			}
		}
	}

	cur_offset = (char *)&G->multiply[G->order];
	// Initialize the multiplication and inverse tables.
	for(i = 0; i < G->order; ++i){
		const group_elem g = G->elems[i];
		group_elem_id j;

		// Initialize the row pointers for each group element.
		G->multiply[i] = (group_elem_id *)cur_offset;
		cur_offset += row_bytes;

		for(j = 0; j < G->order; ++j){
			G->multiply[i][j] = group_elem_multiply(G, g, G->elems[j]);
		}

		G->inverse[i] = group_elem_invert(G, g);
	}

	return 1;
}