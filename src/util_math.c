#include "util_math.h"


// Compute the greatest common divisor of a and b.
integer_t gcd(integer_t a, integer_t b){
	while(a != 0){
		const integer_t c = a;
		a = b % a;
		b = c;
	}
	return b;
}

// Compute a^p in a very naive way. We assume that p >= 0.
integer_t powi(const integer_t a, integer_t p){
	integer_t power = 1;
	for(; p > 0; --p){
		power *= a;
	}
	return power;
}