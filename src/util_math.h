#ifndef UTIL_MATH_H
#define UTIL_MATH_H


typedef int integer_t;
typedef double number_t;


integer_t gcd(integer_t a, integer_t b);
integer_t xgcd(
	integer_t a, integer_t b,
	integer_t *const x, integer_t *const y
);
integer_t powi(const integer_t a, const integer_t p);


#endif