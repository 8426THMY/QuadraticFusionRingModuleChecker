#include <stdio.h>

#include "hi_ring.h"


int main(int argc, char **argv){
	int i;
	for(i = 1; i < argc; ++i){
		hi_ring R;
		if(hi_ring_load(&R, argv[i])){
			hi_ring_find_modules(&R);
			printf("Done!");
		}else{
			printf("Failed to load fusion ring from file \"%s\"!", argv[i]);
		}
		hi_ring_delete(&R);
	}

	return 0;
}