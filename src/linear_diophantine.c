#include "linear_diophantine.h"


#include <stdio.h>


int linear_diophantine_start(linear_diophantine *const dio){
	if(dio->num_eqs == 1){
		// The one or two variable cases are much easier to handle.
		if(dio->num_vars == 1){
			dio->x[0] = dio->b[0]/dio->a[0][0];
			dio->g = dio->a[0][0];
			return (dio->x[0]*dio->a[0][0] == dio->b[0]);
		}else if(dio->num_vars == 2){
			integer_t t;

			// We would like to solve the equation
			//     a[0]x'[0] + a[1]x'[1] = b,
			// in terms of non-negative integers x[0] and x[1].
			// We start by instead solving the equation
			//     a[0]x''[0] + a[1]x''[1] = gcd(a[0], a[1]),
			// using the extended Euclidean algorithm.
			dio->g = xgcd(
				dio->a[0][0], dio->a[0][1],
				&dio->x[0], &dio->x[1]
			);

			t = dio->b[0]/dio->g;
			// A solution to the original equation only exists if the
			// gcd of the subgroup orders divides the right-hand side.
			if(t*dio->g != dio->b[0]){
				return 0;
			}
			// In this case, we obtain a solution by rescaling
			//     x'[0] = (b/gcd(a[0], a[1]))x''[0],
			//     x'[1] = (b/gcd(a[0], a[1]))x''[1].
			dio->x[0] *= t;
			dio->x[1] *= t;

			// We now have integers x'[0] and x'[1] satisfying
			//     a[0]x'[0] + a[1]x'[1] = b.
			// The remaining solutions are given by
			//     x[0] = x'[0] + (a[1]/gcd(a[0], a[1]))t,
			//     x[1] = x'[1] - (a[0]/gcd(a[0], a[1]))t,
			// for some integer t. We'll just return the
			// non-negative solution that minimizes x[0].
			// If x[0] is negative, we need to do some
			// tricks to make sure this still rounds down.
			t = (dio->x[0] >= 0) ?
				(dio->x[0]*dio->g)/dio->a[0][1] :
				(dio->x[0]*dio->g - (dio->a[0][1] - 1))/dio->a[0][1];
			dio->x[0] -= (dio->a[0][1]/dio->g)*t;
			dio->x[1] += (dio->a[0][0]/dio->g)*t;

			// Note that minimizing x[1] will maximize x[0].
			// If the latter is still negative, there must
			// be no non-negative solutions.
			return (dio->x[1] >= 0);
		}

	// We only solve two equations if there's a single variable.
	}else if(dio->num_eqs == 2 && dio->num_vars == 1){
		dio->x[0] = dio->b[0]/dio->a[0][0];
		dio->g = dio->a[0][0];
		// Make sure we get the same solution from
		// both equations and that it's an integer.
		return (
			dio->b[0]*dio->a[1][0] == dio->b[1]*dio->a[0][0] &&
			dio->x[0]*dio->a[0][0] == dio->b[0]
		);
	}

	// Outside of the cases we handle above, I don't know if
	// there's a better way than just brute force. While we can
	// find solutions using the extended Euclidean algorithm or
	// Hermite normal form, finding the non-negative integer
	// solutions seems hard.
	return 0;
}

int linear_diophantine_next(linear_diophantine *const dio){
	if(dio->num_eqs == 1 && dio->num_vars == 2){
		dio->x[0] += dio->a[0][1]/dio->g;
		dio->x[1] -= dio->a[0][0]/dio->g;
		return (dio->x[1] >= 0);
	}

	return 0;
}