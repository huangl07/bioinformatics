#ifndef BIOIO_H
#define BIOIO_H



#include "bioerror.h"
#include "biomemory.h"
#include "biostring.h"


// define the max length of path
#ifndef MAX_PATH_LENGTH
#define MAX_PATH_LENGTH 512
#endif


#ifndef BIO_FCLOSE
#define BIO_FCLOSE
#define BioFCLOSE(x)	{ fflush((x)); fclose((x)); }
#endif


// define BioIO
typedef struct {
	size_t io_length;
	FILE **fp;
	char **io_buffer;
}BioIO;



// extern function
// processing line
extern inline char *bio_readline(FILE *fp);
extern size_t reading_file_names(const char *listfile, char ***mstr);
extern long get_file_lines(const char *filename);
// io processing
extern BioIO* init_BioIO(const char **path, const char **mode, const size_t n);
extern void free_BioIO(BioIO *bioio);



#endif	// end BIOIO

