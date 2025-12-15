# This is a makefile. There are many like it, but this one is mine.
# We should write an if statement to toggle between the setting for XSEDE and MAC




CC = icc

VERSION = $(shell uname)

ifeq ($(VERSION), Darwin)
$(info ${VERSION})


# MAC
CLTOOLS = /Library/Developer/CommandLineTools/usr
CFLAGS = -I/usr/local/include -I${CLTOOLS}/include/
LDFLAGS = -L/usr/local/lib -L${CLTOOLS}/lib/
DEBUG = -Wall -ansi -std=c11
FLAGS = $(CFLAGS) $(LDFLAGS) $(DEBUG)



else ifeq ($(VERSION), Linux)
$(info ${VERSION})



# XSEDE
DEBUG = -Wall -ansi -std=c11
CFLAGS = -I/opt/gsl/2.1/gnu/include 
LDFLAGS = -L/opt/gsl/2.1/gnu/lib
FLAGS = $(CFLAGS) $(LDFLAGS) $(DEBUG)

endif

NAME = a

Hello: Hello.c
	$(CC) $(FLAGS) Hello.c -lgsl -lgslcblas -lfftw3 -lm -O3 -o $(NAME).out


CompRhoGchi: CompRhoGchi.c
	$(CC) $(FLAGS) CompRhoGchi.c -mkl -lgsl -lgslcblas -lm -O3 -o $(NAME).out


CompRhoGchiParams: CompRhoGchiParams.c
	$(CC) $(FLAGS) CompRhoGchiParams.c -mkl -lgsl -lgslcblas -lm -O3 -o $(NAME).out





clean:
	rm *.out .*swp

