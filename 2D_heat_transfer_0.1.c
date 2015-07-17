#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "Mpiinfo.h"
#include <math.h>
#define MESH 100
//---------Initial Profiles--------------
#define centre
#ifdef centre
  #define RADIUS2 100
#endif
//---------functions---------------------
void laplacian(double *f, double *lap, int M);
void heat_source(double *phi, int M);
void boundary(double *phi, int M);
void global_boundary();
int main(int argc, char *argv[]){
  int my_rank;          // rank of process
  int p;                // number of porcess
  int source;           // rank of sender
  int dest;             // rank of receiver
  int tag = 0;          // tag of message
  char message[100];    // storage for messages
  MPI_Status status;    // return status for receive

  /*Start up MPI*/
  MPI_Init(&argc, &argv);
  /*Find out the number of processes*/
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  if(my_rank != 0){
    /*Create Message*/
    sprintf(message,"Greetings from process %d",my_rank);
    dest = 0;
    /*use strlen+1 so that '\0' gets transmitted*/
    MPI_Send(message, strlen(message)+1, MPI_CHAR,dest,tag, MPI_COMM_WORLD);
  }
  else{
    for(source=1; source < p; source++){
      MPI_Recv(message,100,MPI_CHAR,source,tag,MPI_COMM_WORLD,&status);
      printf("%s\n",message);
    }
  }
  /*Shutdown MPI*/
  MPI_Finalize();
}


void laplacian(double *f, double *lap, int M) {
  long i,j,z;

  for (i=1; i< M -1; i++)
  {
    for (j=1; j< M -1; j++)
    {
      z= i*M + j;
      lap[z] = (f[z-1] + f[z+1] + f[z+M] + f[z-M] -4.0*f[z])*inv_deltax2;
    }
  }
}
void heat_source(double *phi, int M){
  long i,j,z;
  double r;
#ifdef Centre
  for ( i = 0; i < M; i++)
  {
    for ( j=0; j < M; j++)
    {
      r= (i-M*0.5)*(i-M*0.5) + (j-M*0.5)*(j-M*0.5);
      z= i*M + j;
      if(r < RADIUS2){
	      phi_old[z] = 1.0;
      }
      else{
	      phi_old[z] = 0.0;
      }
    }
  }
#endif
}
void global_boundary();
