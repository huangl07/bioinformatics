/* bioio.c
 *
 * IO handling for the kmerErrorCorrect library
 *
 */


#include "bioio.h"




/* Function:	bio_readline()
 *
 * Purpose:	Reading one line from file.
 *                
 * Args:	s - string to abstract
 *
 * Returns:	A pointer to the string of the line.
 * 		NULL on failure.
 *
 */
inline char *bio_readline(FILE *fp) {
	char *s, *line;
	int len, pos, n;

	line = (char *) BioMalloc(sizeof(char) * 128);
	n   = 128;

	/* Simple case 1. We're sitting at EOF, or there's an error.
	 *                fgets() returns NULL, so we return NULL.
	 */
	if (fgets(line, n, fp) == NULL) return NULL;

	/* Simple case 2. fgets() got a string, and it reached EOF.
	 *                return success status, so caller can use
	 *                the last line; on the next call we'll
	 *                return the 0 for the EOF.
	 */
	if (feof(fp)) return line;

	/* Simple case 3. We got a complete string, with \n,
	 *                and don't need to extend the buffer.
	 */
	len = strlen(line);
	if (line[len-1] == '\n') return line;

	/* The case we're waiting for. We have an incomplete string,
	 * and we have to extend the buffer one or more times. Make
	 * sure we overwrite the previous fgets's \0 (hence +(n-1)
	 * in first step, rather than 128, and reads of 129, not 128).
	 */
	pos = n - 1;
	while (1) {
		n  += 128;
		line = (char *) BioRealloc(line, sizeof(char) * n);
		s = line + pos;
		if (fgets(s, 129, fp) == NULL) return line;
		len = strlen(s);
		if (s[len-1] == '\n') return line;
		pos += 128;
	}
}




/* Function:	reading_file_names()
 *
 * Purpose:	Reading the content of file into mstr.
 * 		here the empty line will be ignore
 * 
 * Args:	infile - the name of the file
 * 		mstr - the address of char **;
 * 		one line is a record, char *;
 *
 * Returns:	the sum number of records
 *
 */
size_t reading_file_names(const char *infile, char ***mstr) {
	// open listfile
	FILE *fp;
        if( (fp=fopen(infile, "r")) == NULL )
                BioDie("open file %s is failed!", infile);

	// get the number of non-empty line in file
	char *line;
	int sum = 0;
	while ( (line=bio_readline(fp)) != NULL ){
		if( ! bio_is_blank_line(line) ) sum++;
		free(line); line=NULL;
		continue;
	}

	// init seqFileNames
	*mstr = (char **) BioMalloc(sizeof(char *) * sum);

	// reading file
	rewind(fp);	// back to the beginning of the file
	int count = 0;
	while ( (line=bio_readline(fp)) != NULL ){
		if( bio_is_blank_line(line) ) continue;
		bio_chomp(line);
		(*mstr)[count] = line;
		count++;
	}

	// over
	BioFCLOSE(fp);
	return sum;
}





/* Function:	get_file_lines()
 *
 * Purpose:	count the sum number of line in file.
 * 		here the empty line will be not ignored
 * 
 * Args:	infile - the name of the file
 *
 * Returns:	the sum number of line
 * 
 */
long get_file_lines(const char *infile){
//long get_file_lines(char *infile){
	char ch;
	long sum=0;

	// open listfile
	const char *mode[] = {"r"};      // this is a pointer array; [] > *
	const char *path[] = {infile};
	BioIO *bioio = init_BioIO(&path[0], &mode[0], 1);

	// iteration get '\n'
	while( (ch=fgetc(bioio->fp[0])) != EOF )
		if(ch=='\n') sum++;

	// check the end of file
	// check if the last line has '\n'
	// e.g. "aaaa\nEOF" or "aaaaaEOF"
	// get file size of end check 
	long size = ftell(bioio->fp[0]);
	if( size!=0 ){	// file is empty
		fseek(bioio->fp[0], size-1, SEEK_SET);
		ch = fgetc(bioio->fp[0]);
		if( ch != '\n' ) sum++;
	}

	// over
	free_BioIO(bioio); bioio = NULL;
	return sum;
}




/* Function:	init_BioIO()
 *
 * Purpose:	open multiple file handle and set io buffer
 * 
 * Args:	path - the file array
 * 		mode - the mode of file handle
 * 		n - the number of file
 *
 * Returns:	BioIO
 * 
 */
BioIO* init_BioIO(const char **path, const char **mode, const size_t n){
	// allocate memory for BioIO
	BioIO *bioio = (BioIO *) BioMalloc(sizeof(BioIO));
	bioio->io_length = n;

	// create io handle
	// allocate memory for fq
	bioio->fp = (FILE **) BioMalloc(sizeof(FILE *) * n);
	// iteration open file
	int i;
	for( i=0; i<n; i++ )
		if( (bioio->fp[i]=fopen(path[i], mode[i])) == NULL )
			BioDie("open or create file %s is failed!", path[i]);

	// set io buffer
	// allocate memory for io_buffer
	bioio->io_buffer = (char **) BioMalloc(sizeof(char *) * n);
	for(i=0; i<n; i++) bioio->io_buffer[i] = (char *) BioMalloc(IO_BUFFER_SIZE * SIZE_M * sizeof(char));
	// set io buffer
	for(i=0; i<n; i++) setvbuf(bioio->fp[i], bioio->io_buffer[i], _IOFBF, IO_BUFFER_SIZE * SIZE_M);

	// over and return
	return bioio;
}





/* Function:	free_BioIO()
 *
 * Purpose:	free BioIO
 * 
 * Args:	bioio - the pointer to BioIO
 *
 * Returns:	none
 * 
 * NOTE:	here bioio = NULL; is not usable;
 * 		So if you want to set bioio to NULL, you can:
 * 		free_BioIO(bioio);
 * 		bioio = NULL;
 */
void free_BioIO(BioIO *bioio){
	// NOTE: must close fp firstly, then free io buffer
	// close fp
	int i;
	for(i=0; i<bioio->io_length; i++) BioFCLOSE(bioio->fp[i]);
	// free fp
	free(bioio->fp); bioio->fp = NULL;

	// free io buffer
	for(i=0; i<bioio->io_length; i++) {
		free(bioio->io_buffer[i]);
		bioio->io_buffer[i] = NULL; 
	}
	free(bioio->io_buffer); bioio->io_buffer = NULL;

	// free bioio;
	// NOTE: here bioio = NULL; is not usable;
	// So if you want to set bioio to NULL, you can:
	// free_BioIO(bioio);
	// bioio = NULL;
	free(bioio);
}






