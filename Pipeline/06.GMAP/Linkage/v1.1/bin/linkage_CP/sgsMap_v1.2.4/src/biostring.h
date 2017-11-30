#ifndef BIOSTRING_H
#define BIOSTRING_H



#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include "bio.h"
#include "bioerror.h"
#include "biomemory.h"
#include "bioregex.h"

// define 
//#define join(x, ...)  bio_strjoin((x), ## __VA_ARGS__);



// extern function
// check
extern inline int bio_is_blank_line(const char *s);
extern inline int bio_is_contain_alnum(const char *s);
// chomp
extern inline void bio_chomp(char *s);
// strcat, substr and strcpy
extern inline char *bio_strcat(const char *s1, const char *s2);
extern inline char *bio_substr(const char *s, int start, int end);
extern inline char *bio_strcpy(const char *s);
// lower and upper string
extern inline void bio_strlower(char *s);
extern inline void bio_strupper(char *s);
// reverse string
extern char *bio_revstr(const char *s);
extern inline void bio_revstr_in_place(char *s);
// split string
extern size_t bio_strsplit(const char *s, const char *pattern, char ***mstr, const int iscase);
// join string
//extern inline char *bio_strjoin(const char *cc, const char *fs, ...);
// generate rand string
extern inline char* bio_rand_str(size_t len);
// generate random string
extern char *bio_random_string(const gsl_rng * r, const char *s, size_t size);
extern char *bio_random_string_prob(const gsl_rng * r, const char *s, size_t size, const double p[], size_t K);
extern char *bio_random_string_count(const gsl_rng * r, const char *s, const unsigned int count[], size_t K);



#endif	// end BIOSTRING_H

