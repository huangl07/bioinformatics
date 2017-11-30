#ifndef BIOERROR_H
#define BIOERROR_H



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


// define
#define BioDie(x, ...)		bio_die(__FILE__,__LINE__, (x), ## __VA_ARGS__)



// extern function
// die and warn
extern void bio_die(char *file, int line, char *format, ...);



#endif	// end BIOERROR_H

