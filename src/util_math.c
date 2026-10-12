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

// Compute the greatest common divisor of a and b
// as well as x and y such that ax + by = gcd(a, b).
integer_t xgcd(
	integer_t a, integer_t b,
	integer_t *const x, integer_t *const y
){

	integer_t x_new = 1;
	integer_t y_new = 0;
	*x = 0;
	*y = 1;

	while(a != 0){
		const integer_t q = b/a;
		integer_t temp;

		temp = a;
		a = b - q*a;
		b = temp;

		temp = x_new;
		x_new = *x - q*x_new;
		*x = temp;

		temp = y_new;
		y_new = *y - q*y_new;
		*y = temp;
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