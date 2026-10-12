#ifndef QUADRATIC_RING_H
#define QUADRATIC_RING_H


#include <stddef.h>

#include "group.h"
#include "util_math.h"


// A simple module over a quadratic fusion ring consists
// of two pieces of data: a set of subgroups {H_i}_i,
// which tells us the group action of the ring, and a
// set of multiplicities {a_{i,j}}_{i,j}, which tells
// us the action of the non-invertible elements. The
// algebra object corresponding to a group orbit is
// given by inner automorphisms of the objects
//     A_i = H_i \oplus \bigoplus_{gQ \in G/Q}{X_{gQ}^{a_{i,i,g\psi(H_i Q)H_i Q}}}.
typedef struct quadratic_module {
	size_t *H_ids;
	integer_t ***a;
	size_t num_orbits;
} quadratic_module;

// A quadratic fusion ring has simple objects
// G \sqcup {X_{gQ}}_{gQ \in G/Q} and fusion rules
//     g \otimes X_Q = X_{gQ},
//     X_Q \otimes g' = X_{\psi(g'Q)},
//     X_Q^{*} = X_{r(Q)},
//     X_Q \otimes X_Q^{*} = r(Q) \oplus \bigoplus_{g \in G}{X_{gQ}^{m_{gQ}}.
// Here, Q is a subgroup of G, \psi is a group
// automorphism of G/Q, r is a bijection of sets
// G/Q and each m_{gQ} is a non-negative integer.
// We assume that G is a finite Abelian group.
typedef struct quadratic_ring {
	group G;
	size_t Q_id;
	const subgroup *Q;
	// This is technically a group automorphism G/Q -> G/Q,
	// but we store it as a map G -> G so it's easier to use
	// with how we store our groups. An element g in G is
	// sent to any element in the coset \psi(gQ).
	group_elem_id *psi;
	// This is technically a bijection G/Q -> G/Q,
	// but we only need to store an element in the
	// coset r(Q) since r(gQ) = r(Q)\psi(g^{-1}Q).
	integer_t r;
	// This is technically a function G/Q -> Z_{>=0},
	// but we store it as a function G -> Z_{>=0} so
	// it's easier to use with how we store groups.
	integer_t *m;
	// Store the sum of every value of m.
	integer_t m_sum;
	// Store d = FPdim(X_Q), or d = 0 if it's irrational.
	integer_t d;

	// Stores the ID of the subgroup \psi(H_i Q)H_j Q
	// for each pair of subgroups H_i and H_j.
	size_t **psi_ids;

	quadratic_module *modules;
	size_t num_modules;
} quadratic_ring;


int quadratic_ring_load(quadratic_ring *const R, const char *const path);
size_t quadratic_ring_compute_inner_factors(
	const quadratic_ring *const R,
	group_elem_id *const inner_factors
);
void quadratic_ring_delete(quadratic_ring *const R);

int quadratic_module_init(
	quadratic_module *const M, const quadratic_ring *const R,
	const size_t *const H_ids, integer_t ***const a, const size_t num_orbits
);
void quadratic_module_print(
	const quadratic_module *const M,
	const quadratic_ring *const R
);
void quadratic_module_delete(quadratic_module *const M);


#endif