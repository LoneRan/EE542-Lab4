EXECS=ratetest 
MPICC?=mpicc

all: ${EXECS}

ratetest: ratetest.c
	${MPICC} -o ratetest ratetest.c
Mratetest: Mratetest.c
	${MPICC} -o Mratetest Mratetest.c

clean:
	rm -f ${EXECS}
