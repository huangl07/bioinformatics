#ifndef BIOHASH_H
#define BIOHASH_H


#include "bio.h"
#include "bioerror.h"
#include "biomemory.h"
#include "biostring.h"
#include "uthash.h"




//----------------------------
// HASH: string <==> string
//----------------------------

// define the hash struct for ( string <==> string )
typedef struct {
	char *key;
	char *value;
	UT_hash_handle hh;
}BioStrStrRecord;

// define the BioStrStrHASH class
typedef struct {
	// slots
	BioStrStrRecord *biohash;	// hash body

	// methods
	void (*addrecord)(void *, const char *, const char *);
	void (*changevalue)(void *, const char *, const char *);
	BioStrStrRecord* (*findrecord)(void *, const char *);
	void (*deleterecord)(void *, BioStrStrRecord *);
	void (*showrecords)(void *);
	void (*deleteallrecord)(void *);
	void (*sortbykey)(void *);
	void (*sortbyvalue)(void *);
	size_t (*recordsnum)(void *);
	size_t (*keys)(void *, char ***);
	void (*freerecord)(BioStrStrRecord **);
	void (*freehash)(void **);
}BioStrStrHASH;

// define the constructor function
extern int newBioStrStrHASH(void *p);







//----------------------------
// HASH: string <==> int
//----------------------------

// define the hash struct for ( string <==> int )
typedef struct {
	char *key;
	int value;
	UT_hash_handle hh;
}BioStrIntRecord;

// define the BioStrIntHASH class
typedef struct {
	// slots
	BioStrIntRecord *biohash;	// hash body

	// methods
	void (*addrecord)(void *, const char *, int);
	void (*changevalue)(void *, const char *, int);
	BioStrIntRecord* (*findrecord)(void *, const char *);
	void (*deleterecord)(void *, BioStrIntRecord *);
	void (*showrecords)(void *, FILE *);
	void (*deleteallrecord)(void *);
	void (*sortbykey)(void *);
	void (*sortbyvalue)(void *);
	size_t (*recordsnum)(void *);
	size_t (*keys)(void *, char ***);
	void (*freerecord)(BioStrIntRecord **);
	void (*freehash)(void **);
}BioStrIntHASH;

// define the constructor function
extern int newBioStrIntHASH(void *p);




//----------------------------
// HASH: string <==> (void *)
//----------------------------

// define the hash struct for ( string <==> (void *) )
typedef struct {
	char *key;
	void *value;
	UT_hash_handle hh;
}BioStrPointerRecord;

// define the BioStrPointerHASH class
typedef struct {
	// slots
	BioStrPointerRecord *biohash;	// hash body

	// methods
	void (*addrecord)(void *, const char *, void *);
	BioStrPointerRecord* (*findrecord)(void *, const char *);
	void (*deleterecord)(void *, BioStrPointerRecord *);
	void (*deleteallrecord)(void *);
	void (*sortbykey)(void *);
	size_t (*recordsnum)(void *);
	size_t (*keys)(void *, char ***);
	void (*freerecord)(void *, BioStrPointerRecord **);
	void (*freevaluecontent)(void **);
	void (*freehash)(void **);
}BioStrPointerHASH;

// define the constructor function
extern int newBioStrPointerHASH(void *p, void (*_free_value_content)(void **));


#endif	//end BIOHASH_H


