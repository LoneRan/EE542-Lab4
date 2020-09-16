EXECS=ratetest 
MPICC?=mpicc

all: ${EXECS}

ratetest: ratetest.c
	${MPICC} -o ratetest ratetest.c

clean:
	rm -f ${EXECS}