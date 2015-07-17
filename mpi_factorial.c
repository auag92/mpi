#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MASTER 0
#define BEGIN 111
#define WRITE 222
int my_rank;          // rank of process
int numtasks;         // number of porcess
int numworkers;       // no. of worker nodes
int rank;
int source;           // rank of sender
int dest;             // rank of receiver
int tag = 0;          // tag of message
int message;          // storage for messages
MPI_Status status;    // return status for receive
long n;                // Input for factorial calculation
long fact, ans;
long start, end;

long mult(long start, long end){
  int ans;
  if (start != end)
    return start*mult(start-1, end);
  else if(end = 0)
    return 1;
  else
    return end;
}
void get_arg(int *n){
  printf("enter the argument:");
  scanf("%d", n);
  printf("\n");
}
int main(int argc, char *argv[]){

  /*Start up MPI*/
  MPI_Init(&argc, &argv);
  /*Find out the number of processes*/
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  numworkers  = numworkers - 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if(my_rank == 0){
    get_arg(&n);
    int q, rem;
    q   =   n/numworkers;
    rem =   n%numworkers;
    end    =    0;
    start  =    q;
    for (rank=1; rank <= (numworkers); rank++) {
     dest   =    rank;
     end    =    (long)(q*(my_rank-1));
     start  =    (long)(q*my_rank - 1);
     if(  rem >= my_rank ){
       start++;
     }
     MPI_Send(&start,    1,    MPI_LONG,    dest,   BEGIN,  MPI_COMM_WORLD);
     MPI_Send(&end,      1,    MPI_LONG,    dest,   BEGIN,  MPI_COMM_WORLD);
   }
   fact = 1;
   for (rank=1; rank <= numworkers; rank++) {
     source = rank;
     MPI_Recv(&ans,      1,    MPI_LONG,    source, WRITE,  MPI_COMM_WORLD, &status);
     fact = fact*ans;
   }
   printf("The factorial is:%ld",fact);
  }
  else{
    source = MASTER;
    MPI_Recv(&start,      1,   MPI_LONG,    source,  BEGIN,  MPI_COMM_WORLD, &status);
    MPI_Recv(&end,        1,   MPI_LONG,    source,  BEGIN,  MPI_COMM_WORLD, &status);
    ans = mult(start, end);
    MPI_Send(&ans,        1,   MPI_LONG,    MASTER,   BEGIN,  MPI_COMM_WORLD);
  }

  /*Shutdown MPI*/
  MPI_Finalize();
}
