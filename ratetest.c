/*
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

long calTime(struct timeval time1, struct timeval time2)
{
    long elap = (time1.tv_sec - time2.tv_sec) * 1000000 + time1.tv_usec - time2.tv_usec;
    long res = elap ;
    return res;
}
int main(int argc, char **argv)
{
    const int transfer_LIMIT = 2000;

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
    size_t mem_size = 1024 * 1024 * 1024;
    int partner_rank = (world_rank + 1) % 2;
    int *numbers = malloc(mem_size);
    printf("size of MPI_INT is %d\n", sizeof(MPI_INT));
    struct timeval start_time;
    struct timeval end_time;
    for (int i = 1; i < 28; i++)
    {
        transfer_count = 0;
        int transfer_num = 1;
        for (int index = 1; index <= i; index++)
        {
            transfer_num = transfer_num * 2;
        }
        long elap = 0;
        //int numbers[transfer_num];
        if (transfer_count == 0)
        {
            gettimeofday(&start_time, NULL);
        }
        while (transfer_count < transfer_LIMIT)
        {

            if (world_rank == transfer_count % 2)
            {
                // Increment the ping pong count before you send it
                transfer_count++;
                MPI_Send(&transfer_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);

                MPI_Send(numbers, transfer_num, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);

                // printf("node %d sent %d byte(s) to %d\n",
                //     world_rank, transfer_num, partner_rank);
            }
            else
            {
                MPI_Recv(&transfer_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

                MPI_Recv(numbers, transfer_num, MPI_INT, partner_rank, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

                // printf("node %d received %d byte(s) from %d\n",
                //     world_rank, transfer_num, partner_rank);
            }
            
        }
         gettimeofday(&end_time, NULL);
        elap = calTime(end_time, start_time);
        printf("It takes %ld us to transfer %d MPI_INTs back and forth %d times\n", elap, transfer_num, transfer_LIMIT / 2);
        if (elap > 0)
        {

            long throughput = sizeof(MPI_INT) * (long)transfer_num * transfer_LIMIT* 1000000 / (elap*1024*1024);
            printf("Throughput = %ld MB/s\n", throughput);
        }
    }
    free(numbers);

    MPI_Finalize();
}
*/
/* Fox's algorithm on Matrix multiplication ***
   This program was developed taking the help from:
   
   Parallel Programming with MPI by Peter S Pacheco
   Chapter 7: Communicators & Topology
*/

#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include"mpi.h"

#define N 14000 /* dimension of the input matrix */

int matrixA[N][N];
int matrixB[N][N];

typedef struct {
	int p; /* number of processors */
	MPI_Comm comm; /* handle to global grid communicator */
	MPI_Comm row_comm; /* row communicator */
	MPI_Comm col_comm; /* column communicator */
	int q; /* dimension of the grid, = sqrt(p) */
	int my_row; /* row position of a processor in a grid */
	int my_col; /* column position of a procesor in a grid */
	int my_rank; /* rank within the grid */
}GridInfo;


void SetupGrid(GridInfo *grid)
{
	int old_rank;
	int dimensions[2];
	int wrap_around[2];
	int coordinates[2];
	int free_coords[2];
	
	/* get the overall information before overlaying cart_grid */

	MPI_Comm_size(MPI_COMM_WORLD,&(grid->p));
	MPI_Comm_rank(MPI_COMM_WORLD,&old_rank);
	
	/* Assumption: p is a perfect square */
	grid->q=(int)sqrt((double)grid->p);
	/* set the dimensions */
	dimensions[0]=dimensions[1]=grid->q;
	
	/* we want a torus on the second dimension, so set it appropriately */

	wrap_around[0]=0;
	wrap_around[1]=1;
	
	MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,wrap_around,1,&(grid->comm));
	/* since we have set reorder to true, this might have changed the ranks */
	MPI_Comm_rank(grid->comm,&(grid->my_rank));
	/* get the cartesian coordinates for the current process */
	MPI_Cart_coords(grid->comm,grid->my_rank,2,coordinates);
	/* set the coordinate values for the current coordinate */
	grid->my_row=coordinates[0];
	grid->my_col=coordinates[1];

        /* create row communicators */
	free_coords[0]=0;
	free_coords[1]=1; /* row is gonna vary */
	MPI_Cart_sub(grid->comm,free_coords,&(grid->row_comm));
	
        /* create column communicators */
	free_coords[0]=1;
	free_coords[1]=0; /* row is gonna vary */
	MPI_Cart_sub(grid->comm,free_coords,&(grid->col_comm));
	
}

/* normal matrix multiplication stuff */

void matmul(int **a, int **b, int **c, int size)
{
	int i,j,k;
       
	int **temp = (int**) malloc(size*sizeof(int*));
	for(i=0;i<size;i++)
		*(temp+i)=(int*) malloc(size*sizeof(int));

	for(i=0;i<size;i++)
	{
			for(j=0;j<size;j++)
			{
				temp[i][j]=0;
				for(k=0;k<size;k++){
					temp[i][j]=temp[i][j]+ (a[i][k] * b[k][j]);
				}
			}
	}
	
	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
			c[i][j]+=temp[i][j];
	
}
void transfer_data_from_buff(int *buff,int **a,int buffsize, int row, int col){
	
  	if(buffsize!=row*col)
	{
		printf("transfer_data_from_buf: buffer size does not match matrix size!\n");
		exit(1);
	}
	int count=0, i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++){
			a[i][j]=buff[count];
			count++;
		}
	}
}

void transfer_data_to_buff(int *buff,int **a,int buffsize, int row, int col){
	
  	if(buffsize!=row*col)
	{
		printf("transfer_data_to_buf: buffer size does not match matrix size!");
		exit(1);
	}
	int count=0, i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++){
			buff[count]=a[i][j];
			count++;
		}
	}
}

void Fox(int n,GridInfo *grid,int **a, int **b, int **c)
{
	int **tempa;
	int *buff; /* buffer for Bcast & send_recv */
	int stage;
	int root;
	int submat_dim; /* = n/q */
	int source;
	int dest;
	int i;
	MPI_Status status;
	
	submat_dim=n/grid->q;
	
	/* Initialize tempa */
	tempa=(int**) malloc(submat_dim*sizeof(int*));
	for(i=0;i<submat_dim;i++)
		*(tempa+i)=(int*) malloc(submat_dim*sizeof(int));
	/* initialize buffer */
	buff=(int*)malloc(submat_dim*submat_dim*sizeof(int));

        /* we are gonna shift the elements of matrix b upwards with the column fixed */
	source = (grid->my_row+1) % grid->q; /* pick the emmediately lower element */
	dest= (grid->my_row+grid->q-1) % grid->q; /* move current element to immediately upper row */
	
	
	for(stage=0;stage<grid->q;stage++)
	{
		root=(grid->my_col+stage)%grid->q;
		if(root==grid->my_col)
		{
			transfer_data_to_buff(buff,a,submat_dim*submat_dim, submat_dim,submat_dim);
			MPI_Bcast(buff,submat_dim*submat_dim,MPI_INT,root,grid->row_comm);
			transfer_data_from_buff(buff,a,submat_dim*submat_dim, submat_dim,submat_dim);
		
			matmul(a,b,c,submat_dim);
		}else
		{
			transfer_data_to_buff(buff,tempa,submat_dim*submat_dim, submat_dim,submat_dim);
			MPI_Bcast(buff,submat_dim*submat_dim,MPI_INT,root,grid->row_comm);
			transfer_data_from_buff(buff,tempa,submat_dim*submat_dim, submat_dim,submat_dim);
			
			matmul(tempa,b,c, submat_dim);
		}
		transfer_data_to_buff(buff,b,submat_dim*submat_dim, submat_dim,submat_dim);
		MPI_Sendrecv_replace(buff,submat_dim*submat_dim,MPI_INT,dest,0,source,0,grid->col_comm,&status);
		transfer_data_from_buff(buff,b,submat_dim*submat_dim, submat_dim,submat_dim);
	}

}

void initialiseAB()
{
	int i,j;
/* *****************************************************************************************
Initialize the input matrix 
Note: This initalization is deterministic & hence is done by every process in the same way 
      I wanted to design a fully distributed program, hence I took this strategy.
      A better strategy could have been to let master alone initialize the matrices & then
      send the slaves their local copy only 
*******************************************************************************************/
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			matrixA[i][j]=rand() % 10 + 1;
			matrixB[i][j]=rand() % 10 + 1;
		}
	}
}


int main(int argc, char *argv[])
{	
	
	int i,j,dim;
	int **localA;
	int **localB;
	int **localC;
	MPI_Init (&argc, &argv);
	
	GridInfo grid;
	/*initialize Grid */

	SetupGrid(&grid);
        /* Initialize matrix A & B */
	initialiseAB();
        /* calculate local matrix dimension */
	dim=N/grid.q;
	/* allocate space for the three matrices */		

	
	localA=(int**) malloc(dim*sizeof(int*));

	localB=(int**) malloc(dim*sizeof(int*));
	
	localC=(int**) malloc(dim*sizeof(int*));
	
	for(i=0;i<dim;i++)
	{
		*(localA+i)=(int*) malloc(dim*sizeof(int));
		*(localB+i)=(int*) malloc(dim*sizeof(int));
		*(localC+i)=(int*) malloc(dim*sizeof(int));
	}


/* Compute local matrices - Ideally the master should do this & pass it onto all the slaves */
/* At the same time initialize localC to all zeros */

	int base_row=grid.my_row*dim;
	int base_col=grid.my_col*dim;

	for(i=base_row;i<base_row+dim;i++)
	{
		for(j=base_col;j<base_col+dim;j++)
		{
		         localA[i-(base_row)][j-(base_col)]=matrixA[i][j];
			 localB[i-(base_row)][j-(base_col)]=matrixB[i][j];
			 localC[i-(base_row)][j-(base_col)]=0;
		}
	}



	Fox(N,&grid,localA, localB, localC);

/* print results */
	printf("rank=%d, row=%d col=%d\n",grid.my_rank,grid.my_row,grid.my_col);
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			//printf("localC[%d][%d]=%d ", i,j,localC[i][j]);
			printf("%d ", localC[i][j]);
		}
		printf("\n");
	}
	MPI_Finalize ();
	exit(0);

}		
