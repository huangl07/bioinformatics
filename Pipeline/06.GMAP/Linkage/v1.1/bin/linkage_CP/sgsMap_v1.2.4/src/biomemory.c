/* 
 * 
 * memory handling
 * 
 */



#include "biomemory.h"



/* Function:	bio_malloc()
 * 
 * Purpose:	allocates size bytes and returns a pointer to the allocated memory.
 * 		The memory is not initialized.
 * 		calling for bio_detail_die() if the request fails
 * 
 * Args:	file - the name of file: __FILE__
 * 		line - the line of error: __LINE__
 * 		size - the size of the allocated memory
 *           
 * Return:	return a pointer to the allocated memory 
 * 		that is suitably aligned for any kind of variable.
 * 
 */
inline void *bio_malloc(char *file, int line, size_t size) {
	void *ptr;
	if ((ptr = malloc (size)) == NULL)
		BioDie("malloc of %ld bytes failed -->> file %s line %d", size, file, line);
	return ptr;
}



/* Function:	bio_realloc()
 * 
 * Purpose:	changes the size of the memory block pointed to by p to size bytes.
 * 		calling for bio_detail_die() if the request fails
 * 
 * Args:	file - the name of file: __FILE__
 * 		line - the line of error: __LINE__
 * 		p - a pointer to the allocated memory
 * 		size - the size of the allocated memory
 *           
 * Return:	return a pointer to the allocated memory 
 * 		that is suitably aligned for any kind of variable.
 * 
 */
inline void *bio_realloc(char *file, int line, void *p, size_t size) {
	void *ptr;
	if ((ptr = realloc(p, size)) == NULL)
		BioDie("realloc of %ld bytes failed -->> file %s line %d", size, file, line);
	return ptr;
}



/* Function:	bio_free_2D_array(), bio_free_3D_array()
 *
 * Purpose:	Convenience functions for free'ing 2D
 * 		and 3D pointer arrays. Tolerates any of the
 * 		pointers being NULL, to allow "sparse" 
 * 		arrays.
 *
 * Args:	p     - array to be freed
 * 		dim1  - n for first dimension
 * 		dim2  - n for second dimension
 *
 * 		e.g. a 2d array is indexed p[0..dim1-1][]
 * 		a 3D array is indexed p[0..dim1-1][0..dim2-1][]
 *           
 * Returns:	void
 * 
 * Note:	After free, set it NULL
 * 		e.g. free(p); p = NULL;
 * 		bio_free_2D_array(p); p = NULL;
 * 		bio_free_3D_array(p); p = NULL;
 * 
 */
void bio_free_2D_array(void **p, int dim1) {
	int i;

	if (p != NULL) {
		for (i = 0; i < dim1; i++)
			if (p[i] != NULL) free(p[i]);
		free(p);
	}
}
void bio_free_3D_array(void ***p, int dim1, int dim2) {
	int i, j;

	if (p != NULL) {
		for (i = 0; i < dim1; i++)
			if (p[i] != NULL) {
				for (j = 0; j < dim2; j++)
					if (p[i][j] != NULL) free(p[i][j]);
				free(p[i]);
			}
		free(p);
	}
}




