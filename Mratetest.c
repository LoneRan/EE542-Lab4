#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double calTime(struct timeval time1, struct timeval time2)
{
    long elap = (time1.tv_sec - time2.tv_sec) * 1000000 + time1.tv_usec - time2.tv_usec;
    long res = elap / 1000;
    return res;
}
int main(int argc, char **argv)
{
    const int transfer_LIMIT = 100;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    // Find out rank, size
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // We are assuming 2 processes for this task
    if (world_size != 2)
    {
        fprintf(stderr, "World size must be two for %s\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int transfer_count = 0;

    int transfer_num = 1048576;
    int numbers[transfer_num];
    int partner_rank = (world_rank + 1) % 2;

    struct timeval start_time;
    struct timeval end_time;
    int power = 1;
    while(power < (1024/sizeof(MPI_INT))
    {
        long elap = 0;
        transfer_count = 0;

        while (transfer_count < transfer_LIMIT)
        {

            if (transfer_count == 0)
            {
                gettimeofday(&start_time, NULL);
            }
            if (world_rank == transfer_count % 2)
            {
                // Increment the ping pong count before you send it
                transfer_count++;
                MPI_Send(&transfer_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
                for (int j = 0; j < power; j++)
                {
                    MPI_Send(&numbers, transfer_num, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
                }

                // printf("node %d sent %d byte(s) to %d\n",
                //     world_rank, transfer_num, partner_rank);
            }
            else
            {
                MPI_Recv(&transfer_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                for (int j = 0; j < power; j++)
                {
                    MPI_Recv(&numbers, transfer_num, MPI_INT, partner_rank, 0, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                }
                // printf("node %d received %d byte(s) from %d\n",
                //     world_rank, transfer_num, partner_rank);
            }
            if (transfer_count == 100)
            {
                gettimeofday(&end_time, NULL);
            }
            long temp = calTime(end_time, start_time);
            elap += temp;
            power = power * 2;
        }

        printf("It takes %ld ms to transfer %d bytes back and forth %d times\n", elap, transfer_num * power, transfer_LIMIT / 2);
        long throughput = (long)4 * transfer_num * power * transfer_LIMIT * 1000 / elap;
        printf("Throughput = %ld byte/s\n", throughput);
    }

    MPI_Finalize();
}
