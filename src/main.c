#include <stdio.h>

#include "quadratic_ring.h"
#include "module_checker.h"


int main(int argc, char **argv){
	int i;
	for(i = 1; i < argc; ++i){
		quadratic_ring R;
		if(quadratic_ring_load(&R, argv[i])){
			module_checker_find_modules(&R);
			printf("Done!");
		}else{
			printf("Failed to load fusion ring from file \"%s\"!\n", argv[i]);
		}
		quadratic_ring_delete(&R);
	}

	return 0;
}