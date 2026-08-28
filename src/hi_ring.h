#ifndef HI_RING_H
#define HI_RING_H


#include <stddef.h>

#include "group.h"
#include "util_math.h"


// A simple module over a Haagerup-Izumi-like has two
// group orbits. The corresponding algebra objects of
// a particular gorup orbit are given by inner
// automorphisms of objects of the form
//     H \oplus \bigoplus_{gH \in G/H}{X^{\oplus a[gH]}},
// for some subgroup H of G and integers a[gH]. We also
// have an array of integers x indexed by elements of a
// second group H' which tell us how X acts on the orbit.
typedef struct hi_module {
	// Data of the action on the first group orbit.
	// The integers a_1 are indexed by elements of
	// G/H_1, whereas the x_1 are indexed by G/H_2.
	const subgroup *H1;
	const integer_t *a1;
	const integer_t *x1;
	// Data of the action on the second group orbit.
	// The integers a_2 are indexed by elements of
	// G/H_2, whereas the x_2 are indexed by G/H_1.
	const subgroup *H2;
	const integer_t *a2;
	const integer_t *x2;
} hi_module;

// A Haagerup-Izumi-like fusion ring has simple
// objects  G \sqcup G{X} and fusion rules
//     g \otimes X = gX,
//     X \otimes g = \psi(g)X,
//     X \otimes X = 1 \oplus \bigoplus_{q \in Q}{qX}.
// Here, \psi is an involution of G and Q is a subgroup.
typedef struct hi_ring {
	group G;
	const subgroup *Q;
	group_elem_id *psi;

	hi_module *modules;
	size_t num_modules;
} hi_ring;


int hi_ring_load(hi_ring *const R, const char *const path);
void hi_ring_print_module(const hi_ring *const R, const hi_module *const M);
int hi_ring_find_modules(hi_ring *const R);
void hi_ring_delete(hi_ring *const R);


#endif