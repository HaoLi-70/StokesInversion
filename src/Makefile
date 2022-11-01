# Makefile 
#-------------------------------------------

CMP = mpicc
CMPS = gnu


ifeq ($(CMPS), gnu)
ifeq ($(debug), y)
	CFLAGS := -c -fallow-parameterless-variadic-functions -fsanitize=address 
else
ifeq ($(hard), y)
	CFLAGS := -c -O3
else
	CFLAGS := -c
endif
endif
endif





all: MCMC_INV

SRC := $(wildcard *.c)
OBJS := $(SRC:.c=.o)
#OBJ := $(patsubst %.c,%.o,$(wildcard *.c))

%.o: %.c
	$(CMP) $(CFLAGS) -o $@ $< 

MCMC_INV: $(OBJS)
	$(CMP) -o MCMC_INV $(OBJS) -L. -lm
	mv MCMC_INV ../


.PHONY:clean

clean:
	find . -maxdepth 2 -name "*.o" -delete 
	find . -maxdepth 2 -name "MCMC_INV" -delete 