#include "pool_allocator.h"


#include <stdlib.h>
#include <stdio.h>


int pool_init(pool_allocator *const pool, const size_t num_bytes){
	pool->bottom = malloc(num_bytes);
	if(pool->bottom == NULL){
		return 0;
	}

	pool->top = pool->bottom;
	pool->free_bytes = num_bytes;

	return 1;
}

void *pool_alloc(pool_allocator *const pool, const size_t num_bytes){
	void *cur_block;
	if(num_bytes == 0 || num_bytes > pool->free_bytes){
		puts("Pool out of memory!");
		return NULL;
	}

	cur_block = pool->top;

	pool->top = (char *)pool->top + num_bytes;
	pool->free_bytes -= num_bytes;

	return cur_block;
}

void pool_delete(pool_allocator *const pool){
	free(pool->bottom);
}