###########################################
#Makefile for simple programs
###########################################
INC= -I.
LIB= -lstdc++ 
CC= mpic++
CC_FLAG=-Wall -std=c++11

OBJ=
PRG=pf

$(PRG): *.cpp *.h
		$(CC) $(CC_FLAG) $(INC) $(LIB) -o $@ $@.cpp

run:
	mpirun -np 1 $(PRG)
			
.SUFFIXES: .c .o .cpp
	.cpp.o:
		$(CC) $(CC_FLAG) $(INC) -c $*.cpp -o $*.o

.PHONY : clean
clean:
	@echo "Removing linked and compiled files......"
		rm -f $(OBJ) $(PRG) 
		rm -rf testdata/*
		rm -rf testgraph/*
