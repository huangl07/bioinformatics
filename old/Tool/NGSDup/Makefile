headers= ngsduplicate.h
objects= ngsduplicate.o

GCC_OPTION= -m64 -O3


all: ngsduplicate.o $(objects)
	gcc -o ngsduplicate $(objects) $(GCC_OPTION) -lz

ngsdumplicate.o: ngsduplication.c $(headers)
	gcc -c ngsduplicate.c -o ngsduplicate.o $(GCC_OPTION) -lz

.PHONY : clean
clean : 
	rm $(objects) ngsduplicate
