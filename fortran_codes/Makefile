SRCS := $(wildcard *.f90)
BINS := $(SRCS:%.f90=%)

all: ${BINS}

% : %.f90
	gfortran -O3 $^ -o $@
