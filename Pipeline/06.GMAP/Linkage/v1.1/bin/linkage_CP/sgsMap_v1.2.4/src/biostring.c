/* biostring.c
 *
 * string handling for the kmerErrorCorrect library
 *
 */


#include "biostring.h"



/* Function:	bio_is_blank_line()
 *
 * Purpose:	Returns TRUE if string consists solely of whitespace.
 *
 * Args:	s - string to check
 *
 * Return:	TRUE or FALSE <===> 1 or 0
 *
 */
inline int bio_is_blank_line(const char *s) {
	for (; *s != '\0'; s++)
		if (! isspace((int) *s)) return FALSE;
	return TRUE;
}



/* Function:	bio_is_contain_alnum()
 *
 * Purpose:	Returns TRUE if string consists one or more alpha or numeric.
 *
 * Args:	s - string to check
 *
 * Return:	TRUE or FALSE <===> 1 or 0
 *
 */
inline int bio_is_contain_alnum(const char *s) {
	for (; *s != '\0'; s++)
		if (isalnum((int) *s)) return TRUE;
	return FALSE;
}



/* Function:	bio_chomp()
 *
 * Purpose:	Chop trailing whitespace off of a string.
 *
 * Args:	s - string to chop
 *
 * Return:	void
 *
 */
inline void bio_chomp(char *s) {
	int   i;

	i = strlen(s) - 1;		         /* set i at last char in string     */
	while (i >= 0 && isspace((int) s[i])) i--;   /* i now at last non-whitespace char, or -1 */
	s[i+1] = '\0';
}



/* Function:	bio_strcat()
 *
 * Purpose:	Contact two strings, and return a pointer which point to new string.
 *		The new string have new memory which allocated by molloc.
 *
 * Args:	s1 and s2 - string to contact
 *
 * Return:	A pointer to new string. 
 * 		NULL on failure.
 *
 */
inline char *bio_strcat(const char *s1, const char *s2) {
	int len1, len2;
	char *new;

	len1 = (s1 == NULL) ? 0 : strlen(s1);
	len2 = (s2 == NULL) ? 0 : strlen(s2);

	new = (char *) BioMalloc(sizeof(char) * (len1+len2+1));

	if (len1 > 0) memcpy(new, s1, len1);
	if (len2 > 0) memcpy(new+len1, s2, len2);
	new[len1+len2] = '\0';

	return new;
}



/* Function:	bio_substr()
 *
 * Purpose:	Abstract a substring from a string according to start and end.
 *                
 * Args:	s - string to abstract
 *
 * Returns:	A pointer to new string.
 * 		NULL on failure.
 *
 */
inline char *bio_substr(const char *s, int start, int end) {
	char *new;
	int len, range;

	if (s == NULL) return NULL;

	len = strlen(s);
	start = (start < 0) ? 0 : start;
	end = (end > (len-1)) ? (len-1) : end;

	range = end - start;
	if (range < 0) BioDie("substr error: end=%d < start=%d\n", end, start);

	new = (char *) BioMalloc(sizeof(char) * (range+2));	// 2 = a char and '\0'
	memcpy(new, s+start, range+1);
	new[range+1] = '\0';

	return new;
}



/* Function:	bio_strcpy()
 *
 * Purpose:	copies  the  string  pointed to by s
 *                
 * Args:	s - string to cpy
 *
 * Returns:	A pointer to new string.
 * 		NULL on failure.
 *
 */
inline char *bio_strcpy(const char *s) {
	char *new;
	if (s == NULL) return NULL;
	if ((new = (char *) BioMalloc (strlen(s) +1)) == NULL) return NULL;
	strcpy(new, s);
	return new;
}



/* Function:	bio_strlower(); bio_strupper()
 *
 * Purpose:	lower and upper string pointed to by s
 *                
 * Args:	s - a pointer to a string
 *
 * Returns:	void
 *
 */
inline void bio_strlower(char *s) {
	for (; *s != '\0'; s++)
		*s = tolower((int) *s);
}
inline void bio_strupper(char *s) {
	for (; *s != '\0'; s++)
		*s = toupper((int) *s);
}




/* Function:	bio_rand_str()
 *
 * Purpose:	generate a random string, its length is len.
 *
 * Args:	len - the length of string
 *
 * Return:	(char *)
 *
 * Note:	set seed before calling it
 * 		srand((unsigned int)time((time_t *)NULL));
 *		otherwise you will got the same strings
 * 
 */
char* bio_rand_str(size_t len){
	char str[62] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"; 
	int i;
	// allocated memory
	char *rstr;
	rstr = (char *)	BioMalloc(sizeof(char) * (len+1));
	// iteration
	for(i=0; i<len; i++)
		rstr[i] = str[(rand()%62)];
	// add '\0'
	rstr[len] = '\0';
	// over
	return rstr;
}




/* Function:	bio_revstr()
 * 
 * Purpose:	Returns a reversed version of s
 * 
 * Args:	s - string to reverse.
 *           
 * Return:	a pointer to a string
 * 
 */                
char *bio_revstr(const char *s) {
	int  len, pos;
	len = strlen(s);

	char *new;
	new = (char*) BioMalloc( sizeof(char) * (len+1) );

	// reverse
	for (pos = 0; pos < len; pos++)
		new[pos] = s[len-pos-1];
	new[len] = '\0';

	// return
	return new;
}

/* Function:	bio_revstr_in_place()
 * 
 * Purpose:	Returns a reversed version of s
 * 
 * Args:	s - string to reverse.
 *           
 * Return:	a pointer to a string
 * 
 */                
inline void bio_revstr_in_place(char *s) {
	int  len, i, j;
	char ch;
	len = strlen(s);

	// reverse
	for (i=0, j=len-1; i<=j; i++, j--){
		ch = s[i];
		s[i] = s[j];
		s[j] = ch;
	}
}




/* Function:	bio_strsplit()
 * 
 * Purpose:	split a string according to pattern
 * 
 * Args:	s - string to split.
 * 		pattern - a regular expersion
 * 		mstr - a pointer to "char **"
 * 		which is used to store the string frags
 * 		iscase - is ignore case
 * 
 * Return:	the number of the string frags
 * 
 */                
size_t bio_strsplit(const char *s, const char *pattern, char ***mstr, const int iscase) {
	// match split pattern
	size_t matchnum;
	char **tmpstr;
	regmatch_t *reg;
	if (iscase)	matchnum = bio_regex_gmsi(s, pattern, &reg, 1, &tmpstr);
	else		matchnum = bio_regex_gms(s, pattern, &reg, 1, &tmpstr);

	// abstract frags
	*mstr = (char **) BioMalloc(sizeof(char *) * (matchnum+1));
	int pos, i, len;
	pos = 0;
	for (i=0; i<matchnum; i++) {
		if(pos < reg[i].rm_so)		(*mstr)[i] = bio_substr(s, pos, reg[i].rm_so - 1);
		else if(pos == reg[i].rm_so)	(*mstr)[i] = bio_strcpy("\0");
		else				BioDie("start pos (%d) > reg[i].rm_so (%d)", pos, reg[i].rm_so);
		pos = reg[i].rm_eo;	// change the current pos
	}
	// abstract last frag
	len = strlen(s);
	if(pos < len)		(*mstr)[i] = bio_substr(s, pos, len-1);
	else if(pos == len)	(*mstr)[i] = bio_strcpy("\0");
	else			BioDie("start pos (%d) > len (%d)", pos, len);

	// over
	bio_free_2D_array((void **)tmpstr, matchnum); tmpstr = NULL;
	free(reg); reg = NULL;
	return i+1;	// return the number of frags
}


/* Function:	bio_random_string()
 *
 * Purpose:	generate random string from "s" string
 * 
 * Args:	r - random number generator "gls library"
 * 		s - string to sample
 * 		size - the length of random string
 *
 * Returns:	A pointer to new random string.
 * 		NULL on failure.
 *
 */
char *bio_random_string(const gsl_rng * r, const char *s, size_t size) {
	if (s == NULL) BioDie("s is NULL\n");
	size_t len = strlen(s);
	if (len < 1) BioDie("len < 1\n");

	char *old, *new;
	old = bio_strcpy(s);
	new = (char *) BioMalloc (sizeof(char) * (size+1));
	gsl_ran_sample (r, new, size, old, len , sizeof(char));
	new[size] = '\0';

	return new;
}



/* Function:	bio_random_string_prob()
 *
 * Purpose:	generate random string from "s" string
 * 		according to probability distribution "p"
 * 
 * Args:	r - random number generator "gls library"
 * 		s - string to sample
 * 		size - the length of random string
 * 		p - a probability distribution ( sum(p) = 1  )
 * 		K - the length of p, NOTE: K must <= strlen(s)
 *
 * Returns:	A pointer to new random string.
 * 		NULL on failure.
 *
 */
char *bio_random_string_prob(const gsl_rng * r, const char *s, size_t size, const double p[], size_t K) {
	if (s == NULL) BioDie("s is NULL");
	size_t len = strlen(s);
	if (len < 1) BioDie("len < 1");
	if (len < K) BioDie("len < K");

	// check sum(p)
	double sum=0.0;
	int i;
	for(i=0; i<K; i++) sum += p[i];
	if( fabs(sum - 1.0) > 10e-6 ) BioDie("sum(p[]) != 1");

	// The Multinomial Distribution
	unsigned int *n;
	n = (unsigned int *) BioMalloc(sizeof(unsigned int) * K);
	gsl_ran_multinomial (r, K, size, p, n);

	// sample
	return bio_random_string_count(r, s, n, K);
}



/* Function:	bio_random_string_count()
 *
 * Purpose:	generate random string from "s" string
 * 		according to probability distribution "p"
 * 
 * Args:	r - random number generator "gls library"
 * 		s - string to sample
 * 		size - the length of random string
 * 		p - a probability distribution ( sum(p) = 1  )
 * 		K - the length of p, NOTE: K must <= strlen(s)
 *
 * Returns:	A pointer to new random string.
 * 		NULL on failure.
 *
 */
char *bio_random_string_count(const gsl_rng * r, const char *s, const unsigned int count[], size_t K) {
	if (s == NULL) BioDie("s is NULL");
	size_t len = strlen(s);
	if (len < 1) BioDie("len < 1");
	if (len < K) BioDie("len < K");

	// calculate size
	int i,j;
	size_t size = 0;
	for(i=0; i<K; i++) size += count[i];

	// generate indexs
	int *indexs, n=0;
	indexs = (int *) BioMalloc(sizeof(int) * size);
	for(i=0; i<K; i++) 
		for(j=0; j<count[i]; j++)
			indexs[n++] = i;

	// choose random indexs
	gsl_ran_shuffle (r, indexs, size, sizeof(int));
	
	// sample
	char *new;
	new = (char *) BioMalloc (sizeof(char) * (size+1));
	for(i=0; i<size; i++)
		new[i] = s[indexs[i]];
	new[size] = '\0';

	// free
	free(indexs);
	indexs = NULL;

	return new;
}



//#define BIOSTRING_TEST
#ifdef BIOSTRING_TEST

int main(void) {
	// test bio_strsplit()
	char **str, *s, *tmp;
	size_t n, i;
	s = bio_strcpy(",hello,wor\tld,liumin,");
	printf("s: %s\n", s);
	n = bio_strsplit(s, "[,\t]", &str, FALSE);
	for(i=0; i<n; i++) {
		printf("%d\t%s\n", i, str[i]);
	}
	bio_free_2D_array((void **)str, 2); str = NULL;

//	tmp = join("\\@", "ok", "bmk", "bioinformatics", "bio", "", "ok");
//	fprintf(stderr, "tmp: %s\n", tmp);

	return 0;
}

#endif  // end BIOSTRING_TEST
