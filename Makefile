PROG = dssim
MAINSRC = $(PROG).c

CC = gcc 
CFLAGS = -Wall -O3 
#OPTS = -lm  -lfftw -lgsl -lgslcblas -I/usr/local/include/ -I/usr/local/include/gsl/ -L/usr/local/lib/
OPTS = -lm  -lfftw3 -lgsl -lgslcblas -I/usr/include/ -I/usr/include/gsl/ -L/usr/lib/
#OPTS = -lgslcblas -I/share/apps/include/ -L/share/apps/lib/
LIBS = -lm  -lfftw3 -lgsl 

all: 
	$(CC) $(CFLAGS) $(OPTS) *.c -o $(PROG) $(LIBS) 

clean:
	rm *.o
