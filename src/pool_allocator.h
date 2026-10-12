#ifndef POOL_ALLOCATOR_H
#define POOL_ALLOCATOR_H


#include <stddef.h>


// A pool allocator for the heap. This just makes it a little
// easier to manage objects whose memory usage we know in advance.
typedef struct pool_allocator {
	void *top;
	size_t free_bytes;
	void *bottom;
} pool_allocator;


int pool_init(pool_allocator *const pool, const size_t num_bytes);
void *pool_alloc(pool_allocator *const pool, const size_t num_bytes);
void pool_delete(pool_allocator *const pool);


#endif