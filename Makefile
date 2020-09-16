EXECS=ratetest 
MPICC?=mpicc

all: ratetest Mratetest

ratetest: ratetest.c
	${MPICC} -o ratetest ratetest.c
Mratetest: Mratetest.c
	${MPICC} -o Mratetest Mratetest.c

clean:
	rm -f ratetest Mratetest
