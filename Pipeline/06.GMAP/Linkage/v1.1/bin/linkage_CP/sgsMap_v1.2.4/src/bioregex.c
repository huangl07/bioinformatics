/*****************************************************************************
 * Regular Expressions
 * From: http://en.wikipedia.org/wiki/Regular_expression
******************************************************************************/

/*****************************************************************************
 * 
 * POSIX Basic Regular Expressions
 * 
 * Traditional Unix regular expression syntax followed common conventions 
 * but often differed from tool to tool. The IEEE POSIX Basic Regular Expressions (BRE) 
 * standard (released alongside an alternative flavor called Extended Regular Expressions 
 * or ERE) was designed mostly for backward compatibility with the traditional 
 * (Simple Regular Expression) syntax but provided a common standard 
 * which has since been adopted as the default syntax of many Unix regular expression tools, 
 * though there is often some variation or additional features. 
 * Many such tools also provide support for ERE syntax with command line arguments.
 * 
 * In the BRE syntax, most characters are treated as literals — they match only themselves 
 * (e.g., a matches "a"). The exceptions, listed below, are called metacharacters or metasequences.
 * 
 * --------------------------------------------------------------------------------------------
 * Metacharacter	Description 
 * --------------------------------------------------------------------------------------------
 * 	.		Matches any single character (many applications exclude newlines, 
 * 			and exactly which characters are considered newlines is flavor-, character-encoding-, 
 * 			and platform-specific, but it is safe to assume that the line feed character is included). 
 * 			Within POSIX bracket expressions, the dot character matches a literal dot. 
 * 			For example, a.c matches "abc", etc., but [a.c] matches only "a", ".", or "c". 
 * 
 * 	[ ]		A bracket expression. Matches a single character that is contained within the brackets. 
 * 			For example, [abc] matches "a", "b", or "c". 
 * 			[a-z] specifies a range which matches any lowercase letter from "a" to "z". 
 * 			These forms can be mixed: [abcx-z] matches "a", "b", "c", "x", "y", or "z", as does [a-cx-z]. 
 * 			The - character is treated as a literal character if it is the last 
 * 			or the first (after the ^) character within the brackets: [abc-], [-abc]. 
 * 			Note that backslash escapes are not allowed. 
 * 			The ] character can be included in a bracket expression if it is the first (after the ^) character: []abc].
 * 
 * 	[^ ]		Matches a single character that is not contained within the brackets.
 * 			For example, [^abc] matches any character other than "a", "b", or "c". 
 * 			[^a-z] matches any single character that is not a lowercase letter 
 * 			from "a" to "z". As above, literal characters and ranges can be mixed.
 * 
 * 	^		Matches the starting position within the string. 
 * 			In line-based tools, it matches the starting position of any line. 
 * 
 * 	$		Matches the ending position of the string or the position just 
 * 			before a string-ending newline. In line-based tools, 
 * 			it matches the ending position of any line. 
 * 
 * 	BRE: \( \)
 * 	ERE: ( )	Defines a marked subexpression. 
 * 			The string matched within the parentheses can be recalled later 
 * 			(see the next entry, \n). 
 * 			A marked subexpression is also called a block or capturing group. 
 * 
 * 	\n		Matches what the nth marked subexpression matched, 
 * 			where n is a digit from 1 to 9. This construct is theoretically irregular 
 * 			and was not adopted in the POSIX ERE syntax. 
 * 			Some tools allow referencing more than nine capturing groups. 
 * 
 * 	*		Matches the preceding element zero or more times. For example,
 * 			ab*c matches "ac", "abc", "abbbc", etc. 
 * 			[xyz]* matches "", "x", "y", "z", "zx", "zyx", "xyzzy", and so on. 
 * 			\(ab\)* matches "", "ab", "abab", "ababab", and so on. 
 * 
 * 	BRE: \{m,n\}
 * 	ERE: {m,n}	Matches the preceding element at least m and not more than n times. 
 * 			For example, a\{3,5\} matches only "aaa", "aaaa", and "aaaaa". 
 * 			This is not found in a few older instances of regular expressions. 
 * --------------------------------------------------------------------------------------------
 * 
 * Examples:
 * 
 * .at matches any three-character string ending with "at", including "hat", "cat", and "bat". 
 * [hc]at matches "hat" and "cat". 
 * [^b]at matches all strings matched by .at except "bat". 
 * ^[hc]at matches "hat" and "cat", but only at the beginning of the string or line. 
 * [hc]at$ matches "hat" and "cat", but only at the end of the string or line. 
 * \[.\] matches any single character surrounded by "[" and "]" since the brackets are escaped, for example: "[a]" and "[b]". 
 * 
******************************************************************************/


/*****************************************************************************
 *
 * POSIX Extended Regular Expressions
 *
 * The meaning of metacharacters escaped with a backslash is reversed for 
 * some characters in the POSIX Extended Regular Expression (ERE) syntax. 
 * With this syntax, a backslash causes the metacharacter to be treated as a literal character. 
 * So, for example, \( \) is now ( ) and \{ \} is now { }. 
 * Additionally, support is removed for \n backreferences 
 * and the following metacharacters are added:
 * 
 * ------------------------------------------------------------------------------
 * Metacharacter	Description 
 * ------------------------------------------------------------------------------
 * 	?		Matches the preceding element zero or one time. 
 * 			For example, ba? matches "b" or "ba". 
 * 	+		Matches the preceding element one or more times. 
 * 			For example, ba+ matches "ba", "baa", "baaa", and so on. 
 * 	|		The choice (aka alternation or set union) operator matches 
 * 			either the expression before or the expression after the operator. 
 * 			For example, abc|def matches "abc" or "def". 
 * ------------------------------------------------------------------------------
 * 
 * Examples:
 * 
 * [hc]+at matches "hat", "cat", "hhat", "chat", "hcat", "ccchat", and so on, but not "at". 
 * [hc]?at matches "hat", "cat", and "at". 
 * [hc]*at matches "hat", "cat", "hhat", "chat", "hcat", "ccchat", "at", and so on. 
 * cat|dog matches "cat" or "dog". 
 * POSIX Extended Regular Expressions can often be used with modern Unix
 * utilities by including the command line flag -E.
 * 
 * 
******************************************************************************/

/*****************************************************************************
 *
 * POSIX character classes
 *
 * POSIX character classesSince many ranges of characters depend on
 * the chosen locale setting (i.e., in some settings letters are organized as
 * abc...zABC...Z, while in some others as aAbBcC...zZ), 
 * the POSIX standard defines some classes or categories of characters
 * as shown in the following table:
 * 
 * -------------------------------------------------------------------------------------------------------
 * POSIX	Non-standard	Perl	ASCII					Description 
 * -------------------------------------------------------------------------------------------------------
 * [:alnum:]				[A-Za-z0-9]				Alphanumeric characters 
 * 		[:word:]	\w	[A-Za-z0-9_]				Alphanumeric characters plus "_" 
 * 				\W	[^A-Za-z0-9_]				Non-word characters 
 * [:alpha:]				[A-Za-z]				Alphabetic characters 
 * [:blank:]				[ \t]					Space and tab 
 * 				\b	[(?<=\W)(?=\w)|(?<=\w)(?=\W)]		Word boundaries 
 * [:cntrl:]				[\x00-\x1F\x7F]				Control characters 
 * [:digit:]			\d	[0-9]					Digits 
 * 				\D	[^0-9]					Non-digits 
 * [:graph:]				[\x21-\x7E]				Visible characters 
 * [:lower:]				[a-z]					Lowercase letters 
 * [:print:]				[\x20-\x7E]				Visible characters and the space character 
 * [:punct:]				[\]\[!"#$%&'()*+,./:;<=>?@\^_`{|}~-]	Punctuation characters 
 * [:space:]			\s	[ \t\r\n\v\f]				Whitespace characters 
 * 				\S	[^ \t\r\n\v\f]				Non-whitespace characters 
 * [:upper:]				[A-Z]					Uppercase letters 
 * [:xdigit:]				[A-Fa-f0-9]				Hexadecimal digits 
 * --------------------------------------------------------------------------------------------------------
 * 
 * NOTE: POSIX character classes can only be used within bracket expressions.
 * For example, [[:upper:]ab] matches the uppercase letters and lowercase "a" and "b".
 *
******************************************************************************/




#include "bioregex.h"


// static function
// abstract strings
static void _bio_abstract_match_str(const char *s, regmatch_t **pm, const size_t nmatch, const size_t matchnum, char ***mstr);
// m, gm, ogm
static int _bio_regex_m(const char *s, const char *pattern, const int iscase);
static size_t _bio_regex_mt(const char *s, const char *pattern, const int iscase);
static size_t _bio_regex_omt(const char *s, const char *pattern, const int iscase);
// ms, gms, ogms
static int _bio_regex_ms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, const int iscase, char ***mstr);
static size_t _bio_regex_gms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, const int iscase, char ***mstr);
static size_t _bio_regex_ogms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, const int iscase, char ***mstr);
// s, gs, ogs
static char* _bio_regex_s(const char *s, const char *pattern, const char *target, const int iscase);
static char* _bio_regex_gs(const char *s, const char *pattern, const char *target, const int iscase);
static char* _bio_regex_ogs(const char *s, const char *pattern, const char *target, const int iscase);




/* Description
 * 
 * o 	- overlap 
 * g 	- global
 * s 	- substitution		// NOTE order
 * m 	- match
 * t 	- times
 * s	- string		// NOTE order
 * i	- ignore the case
 * 
 * priority: og > s|m > t|ss > i
 * 
 * m		match or isMatch
 * mt		match times
 * omt		match times can overlap
 * ms		match only first, and return string
 * gms		match all, and return string
 * ogms		match all can overlap, and return string
 * s		substitution first
 * gs		substitution all
 * ogs		substitution all can overlap
 * mi		ignore case, match or isMatch
 * mti		ignore case, match times
 * omti		ignore case, match times can overlap
 * msi		ignore case, match only first, and return string
 * gmsi		ignore case, match all, and return string
 * ogmsi	ignore case, match all can overlap, and return string
 * si		ignore case, substitution first
 * gsi		ignore case, substitution all
 * ogsi		ignore case, substitution all can overlap
 * 
 * 
 * NOTE: the site information of match is stored into "regmatch_t **pm"
 * 
 */

/* Args for static function: 
 * _m, _gm, _ogm, _ms, _gms, _oms, _s, _gs and _ogs 
 * 
 * s		- string to check
 * pattern	- regular expression
 * iscase	- differentiate case or not: 0 NO; otherwise YES
 * pm		- a pointer to regmatch_t pointer 
 * 			which will be used to store the site of matchs
 * nmatch	- the number of matchs will be stored 
 * 			(unused structure elements will contain the value -1)
 * target	- a pointer to the string which is used to replace the pattern
 * 
 */

// m, gm, ogm
int bio_regex_m(const char *s, const char *pattern) {
	return _bio_regex_m(s, pattern, FALSE);
}
size_t bio_regex_mt(const char *s, const char *pattern) {
	return _bio_regex_mt(s, pattern, FALSE);
}
size_t bio_regex_omt(const char *s, const char *pattern) {
	return _bio_regex_omt(s, pattern, FALSE);
}
// mi, gmi, ogmi
int bio_regex_mi(const char *s, const char *pattern) {
	return _bio_regex_m(s, pattern, TRUE);
}
size_t bio_regex_mti(const char *s, const char *pattern) {
	return _bio_regex_mt(s, pattern, TRUE);
}
size_t bio_regex_omti(const char *s, const char *pattern) {
	return _bio_regex_omt(s, pattern, TRUE);
}
// ms, gms, ogms
int bio_regex_ms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr) {
	return _bio_regex_ms(s, pattern, pm, nmatch, FALSE, mstr);
}
size_t bio_regex_gms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr) {
	return _bio_regex_gms(s, pattern, pm, nmatch, FALSE, mstr);
}
size_t bio_regex_ogms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr) {
	return _bio_regex_ogms(s, pattern, pm, nmatch, FALSE, mstr);
}
// msi, gmsi, ogmsi
int bio_regex_msi(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr) {
	return _bio_regex_ms(s, pattern, pm, nmatch, TRUE, mstr);
}
size_t bio_regex_gmsi(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr) {
	return _bio_regex_gms(s, pattern, pm, nmatch, TRUE, mstr);
}
size_t bio_regex_ogmsi(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, char ***mstr) {
	return _bio_regex_ogms(s, pattern, pm, nmatch, TRUE, mstr);
}
// s, gs, ogs
char* bio_regex_s(const char *s, const char *pattern, const char *target) {
	return _bio_regex_s(s, pattern, target, FALSE);
}
char* bio_regex_gs(const char *s, const char *pattern, const char *target) {
	return _bio_regex_gs(s, pattern, target, FALSE);
}
char* bio_regex_ogs(const char *s, const char *pattern, const char *target) {
	return _bio_regex_ogs(s, pattern, target, FALSE);
}
// si, gsi, ogsi
char* bio_regex_si(const char *s, const char *pattern, const char *target) {
	return _bio_regex_s(s, pattern, target, TRUE);
}
char* bio_regex_gsi(const char *s, const char *pattern, const char *target) {
	return _bio_regex_gs(s, pattern, target, TRUE);
}
char* bio_regex_ogsi(const char *s, const char *pattern, const char *target) {
	return _bio_regex_ogs(s, pattern, target, TRUE);
}



/* Function:	bio_regcomp()
 *
 * Purpose:	Safe version of regcomp.
 *
 * Args:	file - the name of file call bio_regcomp
 *		line - the number of line call bio_regcomp
 *		preg, regex and cflags - exactly like arguments to regcomp
 *
 * Return:	(void) or die
 *
 */
void bio_regcomp(char *file, int line, regex_t *preg, const char *regex, int cflags) {
	if (regcomp(preg, regex, cflags) != 0)
		BioDie("regcomp failed: file %s line %d", file, line);
}



/* Function:	_bio_regex_m()
 *
 * Purpose:	Returns TRUE if match, otherwise FALSE.
 *
 * Return:	TRUE or FALSE <===> 1 or 0 <===> match or nomatch
 *
 */
int _bio_regex_m(const char *s, const char *pattern, const int iscase) {
	regmatch_t	pm[1];
	regex_t		preg;

	if (iscase) BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE|REG_NOSUB|REG_ICASE);
	else BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE|REG_NOSUB);

	if (regexec(&preg, s, 1, pm, REG_NOTEOL)) {
		regfree(&preg);
		return FALSE;
	}
	regfree(&preg);
	return TRUE;
}



/* Function:	_bio_regex_mt()
 *
 * Purpose:	Return the times of matchs
 *
 * Return:	(size_t) Return the times of matchs
 *
 */
size_t _bio_regex_mt(const char *s, const char *pattern, const int iscase) {
	regmatch_t	pm[1];
	regex_t		preg;
	size_t		matchnum = 0;
	size_t		n = 0;

	if (iscase) BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE|REG_ICASE);
	else BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE);

	while (regexec(&preg, s+n, 1, pm, REG_NOTEOL) == 0) {
		n += pm[0].rm_eo;
		matchnum++;
	}
	regfree(&preg);
	return matchnum;
}



/* Function:	_bio_regex_omt()
 *
 * Purpose:	Return the times of matchs, can overlap
 *
 * Return:	(size_t) Return the times of matchs
 *
 */
size_t _bio_regex_omt(const char *s, const char *pattern, const int iscase) {
	regmatch_t	pm[1];
	regex_t		preg;
	size_t		matchnum = 0;
	size_t		n = 0;

	if (iscase) BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE|REG_ICASE);
	else BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE);

	while (regexec(&preg, s+n, 1, pm, REG_NOTEOL) == 0) {
		n += pm[0].rm_so + 1;
		matchnum++;
	}
	regfree(&preg);
	return matchnum;
}


/* Function:	_bio_abstract_match_str()
 *
 * Purpose:	Abstract string according to match information ("regmatch_t **pm"), 
 * 		and the matched strings will be stored into char ***.
 *
 * Return:	void
 *
 */
void _bio_abstract_match_str(const char *s, regmatch_t **pm, const size_t nmatch, const size_t matchnum, char ***mstr) {
	/* allocate memory for matched strings */
	*mstr = (char **) BioMalloc(sizeof(char *) * nmatch * matchnum);

	size_t i;
	for (i = 0; i < nmatch*matchnum; i++)
		(*mstr)[i] = ((*pm)[i].rm_so == -1) ? NULL : bio_substr(s, (*pm)[i].rm_so, (*pm)[i].rm_eo - 1);

}


/* Function:	_bio_regex_ms()
 *
 * Purpose:	match only once, and return the matched string 
 * 		which will be stored into char ***.
 *		NOTE: here only match once
 *
 * Return:	TRUE or FALSE <===> 1 or 0 <===> match or nomatch
 *
 */
int _bio_regex_ms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, const int iscase, char ***mstr) {
	*pm = (regmatch_t *) BioMalloc(sizeof(regmatch_t) * nmatch);
	regex_t		preg;

	if (iscase) BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE|REG_ICASE);
	else BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE);

	if (regexec(&preg, s, nmatch, *pm, REG_NOTEOL)) {
		regfree(&preg);
		mstr = NULL;	/* null mstr */
		return FALSE;
	}
	regfree(&preg);

	/* abstract matched strings */
	if(mstr != NULL) _bio_abstract_match_str(s, pm, nmatch, 1, mstr);

	return TRUE;
}




/* Function:	_bio_regex_gms()
 *
 * Purpose:	match all, and return the matched string 
 * 		which will be stored into char ***.
 *		NOTE: here all matchs will be found
 * 		NOTE: there is no overlap between two matchs
 *		for example: s="ababac",pattern="aba", here only match once time(s=0,e=3)
 *
 * Return:	Return the times of matchs
 *
 * Example:	regmatch_t *reg;
 *		size_t nmatch = 2, matchnum, i;
 *		char str[] = "abc def ABC";
 *		
 *		matchnum = _bio_regex_gms(str, "(ab)c", &reg, nmatch, 0);
 *		
 *		for (i = 0; i < nmatch*matchnum; i++) {
 *			printf("s=%d,e=%d,", reg[i].rm_so, reg[i].rm_eo);
 *			printf("str=%s\n", bio_substr(str, reg[i].rm_so, reg[i].rm_eo - 1));
 *		}
 *
 */
size_t _bio_regex_gms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, const int iscase, char ***mstr) {
	size_t matchnum = _bio_regex_mt(s, pattern, iscase);
	if (matchnum == 0) {
		*pm = NULL;
		mstr = NULL;	/* null mstr */
		return 0;
	}

	*pm = (regmatch_t *) BioMalloc(sizeof(regmatch_t) * nmatch * matchnum);
	regex_t		preg;

	if (iscase) BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE|REG_ICASE);
	else BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE);

	size_t		n = 0, matchcount = 0, i;
	while ((matchcount < matchnum) && regexec(&preg, s+n, nmatch, (*pm)+matchcount*nmatch, REG_NOTEOL) == 0) {
		for (i = 0; i < nmatch; i++) {
			if ((*pm)[i+matchcount*nmatch].rm_so == -1) break;
			(*pm)[i+matchcount*nmatch].rm_so += n;
			(*pm)[i+matchcount*nmatch].rm_eo += n;
		}
		n = (*pm)[0+matchcount*nmatch].rm_eo;
		matchcount++;
	}
	regfree(&preg);

	/* abstract matched strings */
	if(mstr != NULL) _bio_abstract_match_str(s, pm, nmatch, matchnum, mstr);

	return matchnum;
}



/* Function:	_bio_regex_ogms()
 *
 * Purpose:	match all, and return the matched string 
 * 		which will be stored into char ***.
 *		NOTE: here all matchs will be found
 * 		NOTE: there may be overlap between two matchs
 *		for example: s="ababac",pattern="aba", here match twice(s=0,e=3; s=2,e=5)
 *
 * Return:	Return the times of matchs
 *
 * Example:	regmatch_t *reg;
 *		size_t nmatch = 2, matchnum, i;
 *		char str[] = "ababac";
 *		
 *		//matchnum = _bio_regex_gms(str, "(ab)c", &reg, nmatch, 0);	// compare
 *		matchnum = _bio_regex_ogms(str, "(ab)c", &reg, nmatch, 0);	// compare
 *		
 *		for (i = 0; i < nmatch*matchnum; i++) {
 *			printf("s=%d,e=%d,", reg[i].rm_so, reg[i].rm_eo);
 *			printf("str=%s\n", bio_substr(str, reg[i].rm_so, reg[i].rm_eo - 1));
 *		}
 *
 */
size_t _bio_regex_ogms(const char *s, const char *pattern, regmatch_t **pm, const size_t nmatch, const int iscase, char ***mstr) {
	size_t omatchnum = _bio_regex_omt(s, pattern, iscase);
	if (omatchnum == 0) {
		*pm = NULL;
		mstr = NULL;	/* null mstr */
		return 0;
	}

	*pm = (regmatch_t *) BioMalloc(sizeof(regmatch_t) * nmatch * omatchnum);
	regex_t		preg;

	if (iscase) BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE|REG_ICASE);
	else BioRegcomp(&preg, pattern, REG_EXTENDED|REG_NEWLINE);

	size_t		n = 0, matchcount = 0, i;
	while ((matchcount < omatchnum) && regexec(&preg, s+n, nmatch, (*pm)+matchcount*nmatch, REG_NOTEOL) == 0) {
		for (i = 0; i < nmatch; i++) {
			if ((*pm)[i+matchcount*nmatch].rm_so == -1) break;
			(*pm)[i+matchcount*nmatch].rm_so += n;
			(*pm)[i+matchcount*nmatch].rm_eo += n;
		}
		n = (*pm)[0+matchcount*nmatch].rm_so + 1;
		matchcount++;
	}
	regfree(&preg);

	/* abstract matched strings */
	if(mstr != NULL) _bio_abstract_match_str(s, pm, nmatch, omatchnum, mstr);

	return omatchnum;
}





/*
 * Function:	_bio_regex_s()
 *
 * Purpose:	Replace first pattern with target.
 *		NOTE: here only replace once
 *
 * Return:	A pointer to the new string. If unmatch, it is similar to copy string
 *
 */
char* _bio_regex_s(const char *s, const char *pattern, const char *target, const int iscase) {
	regmatch_t	*reg;

	/* unmatch: copy origin string */
	if (!_bio_regex_ms(s, pattern, &reg, 1, iscase, NULL))	// unmatch copy origin string
		return strdup(s);

	/* match: substr */
	size_t		s_len = strlen(s);
	size_t		target_len = strlen(target);
	size_t		match_len = reg[0].rm_eo - reg[0].rm_so;
	char		*new = BioMalloc(s_len - match_len + target_len + 1);

	strncpy(new, s, reg[0].rm_so);
	strncpy(new + reg[0].rm_so, target, target_len);
	strncpy(new + reg[0].rm_so + target_len, s + reg[0].rm_eo, s_len - reg[0].rm_eo);
	new[s_len - match_len + target_len + 1] = '\0';
	return new;
}



/* Function:	_bio_regex_gs()
 *
 * Purpose:	Replace all matched pattern with target.
 *		NOTE: here replace all matched pattern. 
 *		NOTE: there is no overlap between two matched pattern
 *		for example: "ababac","aba","liua" ====> "liuabac"
 *							----###-----
 *
 * Return:	A pointer to the new string. If unmatch, it is similar to copy string
 *
 * Args:	s   - string to check
 *		pattern - regular expression
 *		target - a pointer to the string which is used to replace the pattern
 *		iscase - differentiate case or not: 0 NO; otherwise YES
 *
 */
char* _bio_regex_gs(const char *s, const char *pattern, const char *target, const int iscase) {
	regmatch_t	*reg;
	size_t		matchnum = _bio_regex_gms(s, pattern, &reg, 1, iscase, NULL);

	/* unmatch: copy origin string */
	if (matchnum == 0) return strdup(s);

	/* calculate the sum length of matched pattern */
	size_t		match_len = 0, i;
	for (i = 0; i < matchnum; i++)
		match_len += reg[i].rm_eo - reg[i].rm_so;

	/* calculate the size of memory for new string */
	size_t		s_len = strlen(s);
	size_t		target_len = strlen(target);
	char		*new = BioMalloc(s_len - match_len + target_len*matchnum + 1);

	size_t		n = 0;
	/* header */
	strncpy(new+n, s, reg[0].rm_so);	// copy the fragment of s
	n += reg[0].rm_so;
	strcpy(new + n, target);		// copy the target
	n += target_len;
	/* body */
	for (i = 1; i < matchnum; i++) {
		strncpy(new + n, s + reg[i-1].rm_eo, reg[i].rm_so - reg[i-1].rm_eo);	// copy the fragment of s
		n += reg[i].rm_so - reg[i-1].rm_eo;
		strcpy(new + n, target);						// copy the target
		n += target_len;
	}
	/* tailer */
	strncpy(new + n, s + reg[matchnum-1].rm_eo, s_len - reg[matchnum-1].rm_eo);
	n += s_len - reg[matchnum-1].rm_eo;
	new[n] = '\0';

	/* check */
	if (n != s_len - match_len + target_len*matchnum)
		BioDie("_bio_regex_gs length error!\n");

	return new;
}



/* Function:	_bio_regex_ogs()
 *
 * Purpose:	Like _bio_regex_gs, replace all matched pattern with target.
 *		NOTE: after replacing, there will not exist matched pattern
 *		for example: "ababac","aba","liua" ====> "liuabac"	// _bio_regex_gs
 *							----###-----
 *		for example: "ababac","aba","liua" ====> "liuliuac"	// _bio_regex_ogs
 *							------------
 * 		WARNING:	if target string contain the pattern, 
 *				which will fail into an infinite loop.
 *				so here we will check it.
 *		for example: "ababac","aba","liuaba" ===+> ERROR	// _bio_regex_ogs
 *
 * Return:	A pointer to the new string. If unmatch, it is similar to copy string
 *
 */
char* _bio_regex_ogs(const char *s, const char *pattern, const char *target, const int iscase) {
	/* check if target string contain the pattern */
	if (_bio_regex_m(target, pattern, iscase))
		BioDie("target str '%s' contain the pattern '%s' at %d line in %s\n",
					target, pattern, __LINE__, __FILE__);

	char *new, *tmp;
	new = _bio_regex_s(s, pattern, target, iscase);
	while (_bio_regex_m(new, pattern, iscase)) {
		tmp = new;
		new = _bio_regex_s(new, pattern, target, iscase);
		free(tmp);
	}

	return new;
}


//#define BIOREGEX_TEST
#ifdef BIOREGEX_TEST

int main(void) {
	// bio_regex_m
	char *s1 = "hello liu min liu";
	char *pattern = "LIU";
	if(bio_regex_m(s1, pattern))
		printf("\"%s\" match \"%s\"\n", s1, pattern);
	else printf("\"%s\" does not match \"%s\"\n", s1, pattern);

	// bio_regex_mi
	if(bio_regex_mi(s1, pattern))
		printf("\"%s\" match \"%s\" ignore case\n", s1, pattern);

	// bio_regex_mti
	size_t n;
	n = bio_regex_mti(s1, pattern);
	printf("\"%s\" match \"%s\" ignore case %u times\n", s1, pattern, n);

	// bio_regex_msi
	char **str;
	regmatch_t *reg;
	if (bio_regex_msi("LIu min liu","l(iu)", &reg, 2, &str))
		printf("bio_regex_msi match!\n");
	printf("s=%d, e=%d\n", reg[0].rm_so, reg[0].rm_eo);
	printf("s=%d, e=%d\n", reg[1].rm_so, reg[1].rm_eo);
	printf("%s\t%s\n", str[0], str[1]);
	
	// bio_regex_gmsi
	bio_free_2D_array((void **)str, 2); str = NULL;
	free(reg); reg = NULL;
	if (bio_regex_gmsi("LIu min liu","l(iu)", &reg, 2, &str))
		printf("bio_regex_gmsi match!\n");
	printf("s=%d, e=%d\n", reg[0].rm_so, reg[0].rm_eo);
	printf("s=%d, e=%d\n", reg[1].rm_so, reg[1].rm_eo);
	printf("s=%d, e=%d\n", reg[2].rm_so, reg[2].rm_eo);
	printf("s=%d, e=%d\n", reg[3].rm_so, reg[3].rm_eo);
	printf("%s\t%s\t%s\t%s\n", str[0], str[1], str[2], str[3]);

	// bio_regex_s, bio_regex_gs, bio_regex_ogs
	printf("--- bio_regex_s, bio_regex_gs, bio_regex_ogs ---\n");
	printf("s=%s\n",bio_regex_s("ababac", "aba", "hello"));
	printf("gs=%s\n",bio_regex_gs("liuaba", "aba", "hello"));
	printf("gs=%s\n",bio_regex_gs("ababacabaccc", "aba", "hello"));
	printf("###\tgs=%s\n",bio_regex_gs("ababac", "aba", "liua"));	// compare
	printf("###\togs=%s\n",bio_regex_ogs("ababac", "aba", "liua")); // compare
	//printf("%s\n", bio_regex_ogs("abab", "ab", "hello_ab"));	// ERROR

	// Regular_expression
	bio_free_2D_array((void **)str, 2); str = NULL;
	free(reg); reg = NULL;
	size_t len = bio_regex_gmsi("hello liumin, welcome to china","([[:alnum:]]+)", &reg, 2, &str);	// compare
	//size_t len = bio_regex_gmsi("hello liumin, welcome to china","[[:alnum:]]+", &reg, 2, &str);	// compare
	int i;
	if(len){
		for(i=0; i<len*2; i++) {
			printf("%d\t%s\n", i, str[i]);
		}
	}

 }

#endif	// end BIOREGEX_TEST


