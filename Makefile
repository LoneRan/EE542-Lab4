EXECS=ratetest 
MPICC?=mpicc

all: ratetest 

ratetest: ratetest.c
	${MPICC} -o ratetest ratetest.c


clean:
	rm -f ratetest
