#ifndef BIOMEMORY_H
#define BIOMEMORY_H

#include "bioerror.h"


// define BioFree
#define BioFree(x) 	{ if((x) != NULL) { free(x); (x) = NULL; } }

// define
#define BioMalloc(x)	bio_malloc(__FILE__, __LINE__, (x))
#define BioRealloc(x,y)	bio_realloc(__FILE__, __LINE__, (x), (y))


// allocate memory for 2D array
#define Bio_Allocate_Memory_for_2D_Array(dest, type, dim1, dim2) \
{ \
	dest = (type**) BioMalloc(sizeof(type*) * (dim1));\
	size_t i; \
	for(i=0; i<(dim1); i++) \
		{ (dest)[i] = (type*) BioMalloc(sizeof(type) * (dim2)); }\
}
// init 2D array
#define Bio_Init_2D_Array(dest, init_value, dim1, dim2) \
{ \
	size_t i, j; \
	for(i=0; i<(dim1); i++)\
		for(j=0; j<(dim2); j++)\
			{ (dest)[i][j] = (init_value); }\
}
// free memory for 2D array
#define Bio_Free_Memory_for_2D_Array(dest, dim1) \
{ \
	size_t i; \
	for(i=0; i<(dim1); i++) { free((dest)[i]); (dest)[i] = NULL; }\
	free(dest); dest = NULL;\
}

// init 1D array
#define Bio_Init_1D_Array(dest, init_value, dim1) \
{ \
	size_t i; \
	for(i=0; i<(dim1); i++) { (dest)[i] = (init_value); }\
}







// extern function
// malloc and realloc
extern inline void *bio_malloc(char *file, int line, size_t size);
extern inline void *bio_realloc(char *file, int line, void *p, size_t size);
// free
extern void bio_free_2D_array(void **p, int dim1);
extern void bio_free_3D_array(void ***p, int dim1, int dim2);





#endif	// end BIOMEMORY_H

