/* biohash.c
 *
 * HASH for the bioliuc library
 *
 */



#include "biohash.h"


//----------------------------
// HASH: string <==> string
//----------------------------

// define the static methods
// class methods
static void _str_str_add_record(void *p, const char *key, const char *value);
static void _str_str_change_value(void *p, const char *key, const char *value);
static BioStrStrRecord* _str_str_find_record(void *p, const char *key);
static void _str_str_delete_record(void *p, BioStrStrRecord *r);
static void _str_str_show_records(void *p);
static void _str_str_delete_all_records(void *p);
static void _str_str_sort_by_key(void *p);
static void _str_str_sort_by_value(void *p);
static size_t _str_str_records_num(void *p);
static size_t _str_str_keys(void *p, char ***);
// compare function
static int _str_str_key_sort(BioStrStrRecord *a, BioStrStrRecord *b);
static int _str_str_value_sort(BioStrStrRecord *a, BioStrStrRecord *b);
// define free method
static void _str_str_free_record(BioStrStrRecord **a);
static void _str_str_free_hash(void **a);


//----------------------------
// HASH: string <==> int
//----------------------------

// define the static methods
// class methods
static void _str_int_add_record(void *p, const char *key, int value);
static void _str_int_change_value(void *p, const char *key, int value);
static BioStrIntRecord* _str_int_find_record(void *p, const char *key);
static void _str_int_delete_record(void *p, BioStrIntRecord *r);
static void _str_int_show_records(void *p, FILE *fp);
static void _str_int_delete_all_records(void *p);
static void _str_int_sort_by_key(void *p);
static void _str_int_sort_by_value(void *p);
static size_t _str_int_records_num(void *p);
static size_t _str_int_keys(void *p, char ***);
// compare function
static int _str_int_key_sort(BioStrIntRecord *a, BioStrIntRecord *b);
static int _str_int_value_sort(BioStrIntRecord *a, BioStrIntRecord *b);
// define free method
static void _str_int_free_record(BioStrIntRecord **a);
static void _str_int_free_hash(void **a);


//----------------------------
// HASH: string <==> (void *)
//----------------------------

// define the static methods
// class methods
static void _str_pointer_add_record(void *p, const char *key, void *value);
static BioStrPointerRecord* _str_pointer_find_record(void *p, const char *key);
static void _str_pointer_delete_record(void *p, BioStrPointerRecord *r);
static void _str_pointer_delete_all_records(void *p);
static void _str_pointer_sort_by_key(void *p);
static size_t _str_pointer_records_num(void *p);
static size_t _str_pointer_keys(void *p, char ***);
// compare function
static int _str_pointer_key_sort(BioStrPointerRecord *a, BioStrPointerRecord *b);
// define free method
static void _str_pointer_free_record(void *p, BioStrPointerRecord **a);
static void _str_pointer_free_hash(void **a);






//----------------------------
// HASH: string <==> string
//----------------------------

/* Function:	newBioStrStrHASH()
 *
 * Purpose:	the constructor function for BioStrStrHASH
 *
 * Args:	p - a pointer to BioStrStrHASH
 *
 * Return:	TRUE or FALSE <===> 1 or 0
 *
 */
int newBioStrStrHASH(void *p) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;

	// init slots
	h->biohash = NULL;	// important

	// initialization method
	h->addrecord = _str_str_add_record;
	h->changevalue = _str_str_change_value;
	h->findrecord = _str_str_find_record;
	h->deleterecord = _str_str_delete_record;
	h->showrecords = _str_str_show_records;
	h->deleteallrecord = _str_str_delete_all_records;
	h->sortbykey = _str_str_sort_by_key;
	h->sortbyvalue = _str_str_sort_by_value;
	h->recordsnum = _str_str_records_num;
	h->keys = _str_str_keys;
	h->freerecord = _str_str_free_record;
	h->freehash = _str_str_free_hash;

	return TRUE;
}



/* Args for the class method
 *
 * p - a pointer to BioStrStrHASH
 * key - a pointer to a string
 * value - a pointer to a string
 * r - a pointer to one record
 * a - a pointer to one record
 * b - a pointer to another record
 *
 */

/* adding one record */
void _str_str_add_record(void *p, const char *key, const char *value) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;

	// check repeat
	if ( h->findrecord(h, key) != NULL ) BioDie("add record (%s) failed: key repeat", key);

	// add new record
	BioStrStrRecord *r;
	r = (BioStrStrRecord *) BioMalloc(sizeof(BioStrStrRecord));
	r->key = bio_strcpy(key);
	r->value = bio_strcpy(value);
	//HASH_ADD_PTR(h->biohash, key, r);	// Note: 3 [char[]] <=> [HASH_ADD_PTR];
	HASH_ADD_KEYPTR( hh, h->biohash, r->key, strlen(r->key), r);	// Note: 5 [char *] <=> [HASH_ADD_KEYPTR]
}

/* adding change record */
void _str_str_change_value(void *p, const char *key, const char *value) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	BioStrStrRecord *r;
	if( (r = h->findrecord(p, key)) == NULL ) BioDie("record %s does not exists", key);
	h->deleterecord(p, r);
	h->addrecord(p, key, value);
}

/* find one record by the key */
BioStrStrRecord* _str_str_find_record(void *p, const char *key) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	BioStrStrRecord *r;

	HASH_FIND_STR(h->biohash, key, r);	// Note: HASH_FIND_STR
	return r;
}

/* delete one record */
void _str_str_delete_record(void *p, BioStrStrRecord *r) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	HASH_DEL(h->biohash, r);
	h->freerecord(&r);
}

/* show the hash */
void _str_str_show_records(void *p) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	BioStrStrRecord *r;

	for (r=h->biohash; r != NULL; r=r->hh.next) {
		printf("key: %s\tvalue: %s\n", r->key, r->value);
	}
}


/* delete all records from the HASH */
void _str_str_delete_all_records(void *p) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	BioStrStrRecord *r, *tmp;

	HASH_ITER(hh, h->biohash, r, tmp) {
		HASH_DEL(h->biohash, r);
		h->freerecord(&r);
	}
}

/* sort the hash by key  */
void _str_str_sort_by_key(void *p) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	HASH_SORT(h->biohash, _str_str_key_sort);
}

/* sort the hash by value  */
void _str_str_sort_by_value(void *p) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	HASH_SORT(h->biohash, _str_str_value_sort);
}

/* count the number of the records in the hash  */
size_t _str_str_records_num(void *p) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	return HASH_COUNT(h->biohash);
}

/* abstract the keys of the hash */
size_t _str_str_keys(void *p, char ***mstr) {
	BioStrStrHASH *h = (BioStrStrHASH*)p;
	BioStrStrRecord *r;

	/* allocate memory for keys */
	size_t n = HASH_COUNT(h->biohash);
	*mstr = (char **) BioMalloc(sizeof(char *) * n);

	int count = 0;
	for (r=h->biohash; r != NULL; r=r->hh.next) {
		(*mstr)[count] = bio_strcpy(r->key);
		count++;
	}
	return n;
}

/* compare function for key  */
int _str_str_key_sort(BioStrStrRecord *a, BioStrStrRecord *b) {
	return strcmp(a->key, b->key);
}

/* compare function for value  */
int _str_str_value_sort(BioStrStrRecord *a, BioStrStrRecord *b) {
	return strcmp(a->value, b->value);
}

/* free method: record */
void _str_str_free_record(BioStrStrRecord **a){
	free((*a)->key);
	free((*a)->value);
	free(*a);
	*a = NULL;
}

/* free method: class */
void _str_str_free_hash(void **a){
	( (BioStrStrHASH*)(*a) )->deleteallrecord( (BioStrStrHASH*)(*a) );
	free(*a);
	*a = NULL;
}





//----------------------------
// HASH: string <==> int
//----------------------------

/* Function:	newBioStrIntHASH()
 *
 * Purpose:	the constructor function for BioStrIntHASH
 *
 * Args:	p - a pointer to BioStrIntHASH
 *
 * Return:	TRUE or FALSE <===> 1 or 0
 *
 */
int newBioStrIntHASH(void *p) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;

	// init slots
	h->biohash = NULL;	// important

	// initialization method
	h->addrecord = _str_int_add_record;
	h->changevalue = _str_int_change_value;
	h->findrecord = _str_int_find_record;
	h->deleterecord = _str_int_delete_record;
	h->showrecords = _str_int_show_records;
	h->deleteallrecord = _str_int_delete_all_records;
	h->sortbykey = _str_int_sort_by_key;
	h->sortbyvalue = _str_int_sort_by_value;
	h->recordsnum = _str_int_records_num;
	h->keys = _str_int_keys;
	h->freerecord = _str_int_free_record;
	h->freehash = _str_int_free_hash;

	return TRUE;
}



/* Args for the class method
 *
 * p - a pointer to BioStrIntHASH
 * key - a pointer to a string
 * value - a pointer to a string
 * r - a pointer to one record
 * a - a pointer to one record
 * b - a pointer to another record
 *
 */

/* adding one record */
void _str_int_add_record(void *p, const char *key, int value) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	// check repeat
	if ( h->findrecord(h, key) != NULL ) BioDie("add record (%s) failed: key repeat", key);

	// add new key
	BioStrIntRecord *r;
	r = (BioStrIntRecord *) BioMalloc(sizeof(BioStrIntRecord));
	r->key = bio_strcpy(key);
	r->value = value;
	//HASH_ADD_PTR(h->biohash, key, r);	// Note: 3 [char[]] <=> [HASH_ADD_PTR];
	HASH_ADD_KEYPTR( hh, h->biohash, r->key, strlen(r->key), r);	// Note: 5 [char *] <=> [HASH_ADD_KEYPTR]
}

/* adding one record */
void _str_int_change_value(void *p, const char *key, int value) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	BioStrIntRecord *r;
	if( (r = h->findrecord(p, key)) == NULL ) BioDie("record %s does not exists", key);
	h->deleterecord(p, r);
	h->addrecord(p, key, value);
}

/* find one record by the key */
BioStrIntRecord* _str_int_find_record(void *p, const char *key) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	BioStrIntRecord *r;

	HASH_FIND_STR(h->biohash, key, r);	// Note: HASH_FIND_STR
	return r;
}

/* delete one record */
void _str_int_delete_record(void *p, BioStrIntRecord *r) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	HASH_DEL(h->biohash, r);
	h->freerecord(&r);
}

/* show the hash */
void _str_int_show_records(void *p, FILE *fp) {
        BioStrIntHASH *h = (BioStrIntHASH*)p;
        BioStrIntRecord *r;

        int i = 0;
        for (r=h->biohash; r != NULL; r=r->hh.next) {
                i++;
                //printf("key: %s\tvalue: %d\n", r->key, r->value);
                fprintf(fp, "%-16s%-6d\t", r->key, r->value);
                if (i % 3 == 0) fprintf(fp,"\n");
        }
        fprintf(fp,"\n");
}

/* delete all records from the HASH */
void _str_int_delete_all_records(void *p) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	BioStrIntRecord *r, *tmp;

	HASH_ITER(hh, h->biohash, r, tmp) {
		HASH_DEL(h->biohash, r);
		h->freerecord(&r);
	}
}

/* sort the hash by key  */
void _str_int_sort_by_key(void *p) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	HASH_SORT(h->biohash, _str_int_key_sort);
}

/* sort the hash by value  */
void _str_int_sort_by_value(void *p) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	HASH_SORT(h->biohash, _str_int_value_sort);
}

/* count the number of the records in the hash  */
size_t _str_int_records_num(void *p) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	return HASH_COUNT(h->biohash);
}

/* abstract the keys of the hash */
size_t _str_int_keys(void *p, char ***mstr) {
	BioStrIntHASH *h = (BioStrIntHASH*)p;
	BioStrIntRecord *r;

	/* allocate memory for keys */
	size_t n = HASH_COUNT(h->biohash);
	*mstr = (char **) BioMalloc(sizeof(char *) * n);

	int count = 0;
	for (r=h->biohash; r != NULL; r=r->hh.next) {
		(*mstr)[count] = bio_strcpy(r->key);
		count++;
	}
	return n;
}

/* compare function for key  */
int _str_int_key_sort(BioStrIntRecord *a, BioStrIntRecord *b) {
	return strcmp(a->key, b->key);
}

/* compare function for value  */
int _str_int_value_sort(BioStrIntRecord *a, BioStrIntRecord *b) {
	if(a->value == b->value) return 0;
	return( (a->value > b->value) ? 1 : -1 );
}

/* free method: record */
void _str_int_free_record(BioStrIntRecord **a){
	free((*a)->key);
	free((*a));
	*a = NULL;
}

/* free method: class */
void _str_int_free_hash(void **a){
	( (BioStrIntHASH*)(*a) )->deleteallrecord( (BioStrIntHASH*)(*a) );
	free(*a);
	*a = NULL;
}





//----------------------------
// HASH: string <==> pointer
//----------------------------

/* Function:	newBioStrPointerHASH()
 *
 * Purpose:	the constructor function for BioStrPointerHASH
 *
 * Args:	p - a pointer to BioStrPointerHASH
 * 		_free_value_content - a function pointer which is used to free the value content
 *
 * Return:	TRUE or FALSE <===> 1 or 0
 * 
 * NOTE:	the content which pointer by the value will be free
 *
 */
int newBioStrPointerHASH(void *p, void (*_free_value_content)(void **)) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;

	// init slots
	h->biohash = NULL;	// important

	// initialization method
	h->addrecord = _str_pointer_add_record;
	h->findrecord = _str_pointer_find_record;
	h->deleterecord = _str_pointer_delete_record;
	h->deleteallrecord = _str_pointer_delete_all_records;
	h->sortbykey = _str_pointer_sort_by_key;
	h->recordsnum = _str_pointer_records_num;
	h->keys = _str_pointer_keys;
	h->freerecord = _str_pointer_free_record;
	h->freevaluecontent = _free_value_content;
	h->freehash = _str_pointer_free_hash;

	return TRUE;
}



/* Args for the class method
 *
 * p - a pointer to BioStrPointerHASH
 * key - a pointer to a string
 * value - a pointer to void *
 * r - a pointer to one record
 * a - a pointer to one record
 * b - a pointer to another record
 *
 */

/* adding one record */
void _str_pointer_add_record(void *p, const char *key, void *value) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;

	// check repeat
	if ( h->findrecord(h, key) != NULL ) BioDie("add record (%s) failed: key repeat", key);

	// add new record
	BioStrPointerRecord *r;
	r = (BioStrPointerRecord *) BioMalloc(sizeof(BioStrPointerRecord));
	r->key = bio_strcpy(key);
	r->value = value;
	//HASH_ADD_PTR(h->biohash, key, r);	// Note: 3 [char[]] <=> [HASH_ADD_PTR];
	HASH_ADD_KEYPTR( hh, h->biohash, r->key, strlen(r->key), r);	// Note: 5 [char *] <=> [HASH_ADD_KEYPTR]
}

/* find one record by the key */
BioStrPointerRecord* _str_pointer_find_record(void *p, const char *key) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;
	BioStrPointerRecord *r;

	HASH_FIND_STR(h->biohash, key, r);	// Note: HASH_FIND_STR
	return r;
}

/* delete one record */
void _str_pointer_delete_record(void *p, BioStrPointerRecord *r) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;
	HASH_DEL(h->biohash, r);
	h->freerecord(p, &r);
}

/* delete all records from the HASH */
void _str_pointer_delete_all_records(void *p) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;
	BioStrPointerRecord *r, *tmp;

	HASH_ITER(hh, h->biohash, r, tmp) {
		HASH_DEL(h->biohash, r);
		h->freerecord(p, &r);
	}
}

/* sort the hash by key  */
void _str_pointer_sort_by_key(void *p) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;
	HASH_SORT(h->biohash, _str_pointer_key_sort);
}

/* count the number of the records in the hash  */
size_t _str_pointer_records_num(void *p) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;
	return HASH_COUNT(h->biohash);
}

/* abstract the keys of the hash */
size_t _str_pointer_keys(void *p, char ***mstr) {
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;
	BioStrPointerRecord *r;

	/* allocate memory for keys */
	size_t n = HASH_COUNT(h->biohash);
	*mstr = (char **) BioMalloc(sizeof(char *) * n);

	int count = 0;
	for (r=h->biohash; r != NULL; r=r->hh.next) {
		(*mstr)[count] = bio_strcpy(r->key);
		count++;
	}
	return n;
}

/* compare function for key  */
int _str_pointer_key_sort(BioStrPointerRecord *a, BioStrPointerRecord *b) {
	return strcmp(a->key, b->key);
}

/* free method: record */
/* NOTE: here will free element self, meanwhile free the pointer content of value */
void _str_pointer_free_record(void *p, BioStrPointerRecord **a){
	BioStrPointerHASH *h = (BioStrPointerHASH*)p;

	free((*a)->key); 	// free key
	h->freevaluecontent(&((*a)->value)); 	// NOTE NOTE NOTE
	free(*a); 		// free self
	*a = NULL;		// set self to NULL
}

/* free method: class */
void _str_pointer_free_hash(void **a){
	( (BioStrPointerHASH*)(*a) )->deleteallrecord( (BioStrPointerHASH*)(*a) );
	free(*a);
	*a = NULL;
}

/*
static void _myfree(void **p){
	free(*p);
	*p = NULL;
}
*/

//#define BIOHASH_TEST
#ifdef BIOHASH_TEST
int main(void){
	char *a = bio_strcpy("hello owrld");
	BioStrPointerHASH *pp;
	pp = (BioStrPointerHASH *) BioMalloc(sizeof(BioStrPointerHASH));
	newBioStrPointerHASH(pp, _myfree);

	if(pp->findrecord(pp, "liu") == NULL){
		pp->addrecord(pp, "liu", a);
	}

	printf("ok,sum=%d\n", pp->recordsnum(pp));

	pp->deleteallrecord(pp);

	pp->freehash((void**)(&pp));

	free(a); 	// NOTE here is error; a's memory has been freeed by hashfree function

/*
	// BioStrIntHASH test
	//
	BioStrIntHASH *user;
	user = (BioStrIntHASH *) BioMalloc(sizeof(BioStrIntHASH));
	newBioStrIntHASH(user);

	// add a
	char *hello;
	//hello = bio_strcpy("a");
	if(user->findrecord(user,"a")==NULL){
		user->addrecord(user, "a", 3);
	}
	// add a : failed
	if(user->findrecord(user,"a")==NULL){
		user->addrecord(user, "a", 3);
	}
	// add c
	if(user->findrecord(user,"c")==NULL){
		user->addrecord(user, "c", 2);
	}
	// add b
	if(user->findrecord(user,"b")==NULL){
		user->addrecord(user, "b", 1);
	}

	printf("origin user:\n");
	user->showrecords(user);

	printf("sort by key:\n");
	user->sortbykey(user);
	user->showrecords(user);

	printf("sort by value:\n");
	user->sortbyvalue(user);
	user->showrecords(user);

	printf("delete 'b':\n");
	BioStrIntRecord *r;
	r = user->findrecord(user, "b");
	user->deleterecord(user, r);
	user->showrecords(user);

	char **mstr;
	size_t n = user->keys(user, &mstr);
	int i;
	for(i=0; i<n; i++){
		printf("%d\t%s\n", i, mstr[i]);
	}

	printf("keys: %p\n", user->keys);
	printf("biohash: %p\n", user->biohash);
	BioStrIntHASH *userb;
	userb = (BioStrIntHASH *) BioMalloc(sizeof(BioStrIntHASH));
	newBioStrIntHASH(userb);
	printf("keys: %p\n", userb->keys);
	printf("biohash: %p\n", userb->biohash);

	user->freehash( (void **)(&user) );	// free and NULL

	printf("over!\n");
*/
	
	//*/


	// BioStrPointerHASH test
	/*
	BioStrStrHASH *seqs;
	seqs = (BioStrStrHASH *) BioMalloc(sizeof(BioStrStrHASH));
	newBioStrStrHASH(seqs);

	// add a
	if(seqs->findrecord(seqs,"a")==NULL){
		seqs->addrecord(seqs, "a", "TTTTTTT");
	}
	// add a : failed
	if(seqs->findrecord(seqs,"a")==NULL){
		seqs->addrecord(seqs, "a", "TTTTTTT");
	}
	// add c
	if(seqs->findrecord(seqs,"c")==NULL){
		seqs->addrecord(seqs, "c", "AAAAAAA");
	}
	// add b
	if(seqs->findrecord(seqs,"b")==NULL){
		seqs->addrecord(seqs, "b", "GGGGGGG");
	}

	printf("origin seqs:\n");
	seqs->showrecords(seqs);

	printf("sort by key:\n");
	seqs->sortbykey(seqs);
	seqs->showrecords(seqs);

	printf("sort by value:\n");
	seqs->sortbyvalue(seqs);
	seqs->showrecords(seqs);

	printf("delete 'b':\n");
	BioStrStrRecord *r;
	r = seqs->findrecord(seqs, "b");
	seqs->deleterecord(seqs, r);
	seqs->showrecords(seqs);

	char **mstr;
	size_t n = seqs->keys(seqs, &mstr);
	int i;
	for(i=0; i<n; i++){
		printf("%d\t%s\n", i, mstr[i]);
	}
	*/
}
#endif	// end BIOHASH_TEST


