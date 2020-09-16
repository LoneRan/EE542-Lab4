EXECS=ratetest 
MPICC?=mpicc

all: ${EXECS}

Mratetest: Mratetest.c
	${MPICC} -o Mratetest Mratetest.c

clean:
	rm -f ${EXECS}
