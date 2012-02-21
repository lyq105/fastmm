target = helloworld
CC = gcc
CPP = g++
FROTRAN = ifort
#CFLAGS = -O2
CFLAGS = -g
LIBS = 
LINKER = $(CPP)
LINKFLAG = -I.

# source files
srcfiles 	:= $(wildcard *.cpp)

# object files
objects		:= $(patsubst %.cpp, %.o, $(srcfiles))

all: $(target)

$(target): $(objects)
	@echo "linking the program"
	@$(LINKER) $(LINKFLAG) $(objects) -o$(target) 
	@echo "linking is done!"


# compile rules

.cpp.o:
	@$(CPP) $(CFLAGS) -c $< -o $@
.c.o:
	@$(CC) $(CFLAGS) -c $< -o $@ 

clean:
	@rm *.o 
