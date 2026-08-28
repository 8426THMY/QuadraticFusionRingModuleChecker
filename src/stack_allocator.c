#include "stack_allocator.h"


#include <stdlib.h>


int stack_init(stack_allocator *const stack, const size_t num_bytes){
	stack->bottom = malloc(num_bytes);
	if(stack->bottom == NULL){
		return 0;
	}

	stack->top = stack->bottom;
	stack->free_bytes = num_bytes;

	return 1;
}

void *stack_alloc(stack_allocator *const stack, const size_t num_bytes){
	void *cur_block;
	const size_t total_bytes = num_bytes + sizeof(stack_footer);
	if(num_bytes == 0 || total_bytes > stack->free_bytes){
		return NULL;
	}

	cur_block = stack->top;
	// Store the block's size in its footer.
	*((stack_footer *)((char *)cur_block + num_bytes)) = num_bytes;

	stack->top = (char *)stack->top + total_bytes;
	stack->free_bytes -= total_bytes;

	return cur_block;
}

void stack_free_last(stack_allocator *const stack){
	if(stack->top != stack->bottom){
		const stack_footer free_size =
			// Read the last block's size from its footer.
			*((stack_footer *)((char *)stack->top - sizeof(stack_footer))) +
			// We also need to include the size of the footer!
			sizeof(stack_footer);
		stack->top = (char *)stack->top - free_size;
		stack->free_bytes += free_size;
	}
}

void stack_delete(stack_allocator *const stack){
	free(stack->bottom);
}