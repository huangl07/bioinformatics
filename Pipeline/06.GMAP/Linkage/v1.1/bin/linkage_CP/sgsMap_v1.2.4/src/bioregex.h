#ifndef BIOREGEX_H
#define BIOREGEX_H


#include <regex.h>
#include "bio.h"
#include "bioerror.h"
#include "biomemory.h"
#include "biostring.h"



// define
// regcomp
#define BioRegcomp(x,y,z)  bio_regcomp(__FILE__, __LINE__, (x), (y), (z))



// extern function
// regcomp
extern void bio_regcomp(char *file, int line, regex_t *preg, const char *regex, int cflags);
// m, gm, ogm
extern int bio_regex_m(const char *s, const char *pattern);
extern size_t bio_regex_mt(const char *s, const char *pattern);
extern size_t bio_regex_omt(const char *s, const char *pattern);
// mi, gmi, ogmi
extern int bio_regex_mi(const char *s, const char *pattern);
extern size_t bio_regex_mti(const char *s, const char *pattern);
extern size_t bio_regex_omti(const char *s, const char *pattern);
// ms, gms, ogms
extern int bio_regex_ms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr);
extern size_t bio_regex_gms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr);
extern size_t bio_regex_ogms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr);
// msi, gmsi, ogmsi
extern int bio_regex_msi(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr);
extern size_t bio_regex_gmsi(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr);
extern size_t bio_regex_ogmsi(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr);
// s, gs, ogs
extern char* bio_regex_s(const char *s, const char *pattern, const char *target);
extern char* bio_regex_gs(const char *s, const char *pattern, const char *target);
extern char* bio_regex_ogs(const char *s, const char *pattern, const char *target);
// si, gsi, ogsi
extern char* bio_regex_si(const char *s, const char *pattern, const char *target);
extern char* bio_regex_gsi(const char *s, const char *pattern, const char *target);
extern char* bio_regex_ogsi(const char *s, const char *pattern, const char *target);





#endif	// end BIOREGEX_H

