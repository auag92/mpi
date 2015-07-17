#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MASTER 0

int mult(int start, int end);
void get_arg(int *n);
int main(int argc, char *argv[]){
  int my_rank;          // rank of process
  int p;                // number of porcess
  int source;           // rank of sender
  int dest;             // rank of receiver
  int tag = 0;          // tag of message
  int message;          // storage for messages
  MPI_Status status;    // return status for receive

  int n;                // Input for factorial calculation
  get_arg(&n);

  /*Start up MPI*/
  MPI_Init(&argc, &argv);
  /*Find out the number of processes*/
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  if(my_rank != 0){
    /*Create Message*/
    dest = 0;
    int q = n/p;
    int rem = n - q*p;
    int start,end;
    if (my_rank<p){
      end = q*(my_rank-1) + 1;
      start = q*my_rank;
    }
    else{
      end = q*(my_rank-1) + 1;
      start = n;
    }
    int message = mult(start, end);
    MPI_Send(message, strlen(message)+1, MPI_CHAR,dest,tag, MPI_COMM_WORLD);
  }
  else{
    for(source=1; source < p; source++){
      MPI_Recv(message,100,MPI_CHAR,source,tag,MPI_COMM_WORLD,&status);
      int fact = fact*message;
      printf("factorial is %d\n",fact);
    }
  }
  /*Shutdown MPI*/
  MPI_Finalize();
}

int mult(int start, int end){
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
