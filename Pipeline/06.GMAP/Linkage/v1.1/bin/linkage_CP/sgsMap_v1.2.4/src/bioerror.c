/* 
 * 
 * error information handling and exit
 * 
 */



#include "bioerror.h"



/* Function:	bio_die()
 * 
 * Purpose:	Print an error message and die. The arguments
 * 		are formatted exactly like arguments to printf().
 *           
 * Return:	None. Exits the program.
 * 
 */          
void bio_die(char *file, int line, char *format, ...) {
	// print position: file and line
	fprintf(stderr, "\nFATAL ERROR (file %s line %d): ", file, line);
	/* format the error mesg */
	va_list  argp;
	va_start(argp, format);
	vfprintf(stderr, format, argp);
	va_end(argp);
	fprintf(stderr, "\n");
	fflush(stderr);
	/* exit  */
	exit(1);
}



