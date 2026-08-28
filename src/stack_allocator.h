#ifndef STACK_ALLOCATOR_H
#define STACK_ALLOCATOR_H


#include <stddef.h>


// A stack allocator for the heap. Most of the memory we allocate
// during loops has an upper bound, so by allocating it all once,
// we can avoid a lot of unnecessary system calls.
typedef struct stack_allocator {
	void *top;
	size_t free_bytes;
	void *bottom;
} stack_allocator;

// After each block, we store its size.
typedef size_t stack_footer;


int stack_init(stack_allocator *const stack, const size_t num_bytes);
void *stack_alloc(stack_allocator *const stack, const size_t num_bytes);
void stack_free_last(stack_allocator *const stack);
void stack_delete(stack_allocator *const stack);


#endif