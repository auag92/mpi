#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MESHX 100
#define MASTER 0
#define NONE 0

#define BEGIN 999
#define LTAG 777
#define RTAG 666
#define WRITE 555

#define ntimesteps 10000
#define saveT 500

#define D 1.0
#define deltax 1
#define deltat 1

int numtasks, taskid;
int numworkers;

double *comp;
long start, end, rows, t;
int avgrow, extra, dest, rank, offset;
int left_node, right_node;
int source, msgtype;
int dest;
MPI_Status status;

void mpiexchange(int taskid){
  if((taskid%2)==0){
    if(taskid != numworkers){
      MPI_Send(&comp[end],        1, MPI_DOUBLE, right_node, LTAG,    MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&comp[end+1],      1, MPI_DOUBLE, source,     msgtype, MPI_COMM_WORLD, &status);
    }
    MPI_Send(&comp[start],        1, MPI_DOUBLE, left_node,  RTAG,    MPI_COMM_WORLD);
    msgtype = LTAG;
    MPI_Recv(&comp[0],            1, MPI_DOUBLE, source,     msgtype, MPI_COMM_WORLD, &status);
  }
  else{
    if(taskid != 1){
      source  = left_node;
      msgtype = LTAG;
      MPI_Recv(&comp[0],          1, MPI_DOUBLE, source,    msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&comp[start],      1, MPI_DOUBLE, left_node, RTAG,     MPI_COMM_WORLD);
    }
    if(taskid != numworkers){
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&comp[end+1],      1, MPI_DOUBLE, source,    msgtype,   MPI_COMM_WORLD, &status);
      MPI_Send(&comp[start],      1, MPI_DOUBLE, left_node, RTAG,      MPI_COMM_WORLD);
    }
  }
}
void sendtomaster(){
  MPI_Send(&offset,     1, MPI_INT, dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,       1, MPI_INT, dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,  1, MPI_INT, dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node, 1, MPI_INT, dest, WRITE, MPI_COMM_WORLD);
  if(taskid == 1){
    MPI_Send(&comp[0],  rows, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
  }else{
    MPI_Send(&comp[1],  rows, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
  }
}
void initialize(){
  long i;
  for(i=0; i<MESHX; i++){
    if(i<MESHX/2)
      comp[i]=1.0;
    else
      comp[i]=0.0;
  }
}
void receivefrmworker(){
  int rank;
  for(rank = 1; rank <= numworkers; rank++){
    source  =rank;
    MPI_Recv(&offset,     1, MPI_INT,          source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,       1, MPI_INT,          source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,  1, MPI_INT,          source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node, 1, MPI_INT,          source, WRITE, MPI_COMM_WORLD, &status);
    MPI_Recv(&comp[offset], rows, MPI_DOUBLE,  source, WRITE, MPI_COMM_WORLD, &status);
  }
}
void solverloop(){
  double flux_back, flux_front;
  long i;
  flux_back = (comp[1]-comp[0])/deltax;
  for(i=start; i<=end; i++){
    flux_front     = (comp[i+1]-comp[i])/deltax;
    comp[i]        = comp[i]+deltat*D*(flux_front - flux_back)/deltax;
    flux_back      = flux_front;
  }
}

void writetofile(long t){
  char name[1000];
  long i;
  sprintf(name,"composition_%ld.dat",t);
  FILE *fp = fopen(name,"w");
  for(i = 0; i < MESHX; i++){
    fprintf(fp,"%ld %le\n",i*deltax, comp[i]);
  }
  fclose(fp);
}
void apply_boundary_conditions(){
  if(taskid == 0){
    comp[0] = comp[1];
  }
  else if(taskid == numworkers){
    comp[end + 1] = comp[end];
  }
}

int main( int argc, char *argv[]){
  initialize();
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  numworkers = numtasks - 1;
  if(taskid == MASTER){

    comp = (double *)malloc(MESHX*sizeof(double));

    avgrow = MESHX/numworkers;
    extra  = MESHX%numworkers;
    offset = 0;
    for (rank = 1; rank <=numworkers; rank++){
      rows = (rank <= extra) ? avgrow+1 : avgrow;
      if (rank == 1)
      left_node  =  NONE;
      if (left_node == (numworkers)){
        right_node = NONE;
      }

      MPI_Send(&offset,     1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&rows,       1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&left_node,  1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&right_node, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);

      MPI_Send(&comp[offset],  rows, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
      offset = offset + rows;
    }
    for(t=1; t<ntimesteps; t++){
      if (t%saveT == 0){
        receivefrmworker();
        writetofile(t);
      }
    }
    free(comp);
  }
  else {
    source    = MASTER;
    msgtype   = BEGIN;
    MPI_Recv(&offset,     1, MPI_INT, source,msgtype,MPI_COMM_WORLD,&status);
    MPI_Recv(&rows,       1, MPI_INT, source,msgtype,MPI_COMM_WORLD,&status);
    MPI_Recv(&left_node,  1, MPI_INT, source,msgtype,MPI_COMM_WORLD,&status);
    MPI_Recv(&right_node, 1, MPI_INT, source,msgtype,MPI_COMM_WORLD,&status);

    start = 1;
    if (taskid == 1){
      end = rows - 1;
      comp = (double *)malloc((rows+1)*sizeof(double));
      MPI_Recv(&comp[0], rows, MPI_DOUBLE, source,msgtype,MPI_COMM_WORLD,&status);
    }
    else if(taskid == numworkers){
      end = rows - 1;
      comp = (double *)malloc((rows+1)*sizeof(double));
      MPI_Recv(&comp[1], rows, MPI_DOUBLE, source,msgtype,MPI_COMM_WORLD,&status);
    }
    else{
      end = rows;
      comp = (double *)malloc((rows+2)*sizeof(double));
      MPI_Recv(&comp[1], rows, MPI_DOUBLE, source,msgtype,MPI_COMM_WORLD,&status);
    }
    for(t=1; t<ntimesteps; t++){
      mpiexchange(taskid);
      solverloop(comp);
      apply_boundary_conditions(taskid);
      if (t%saveT == 0){
        sendtomaster();
      }
    }
    free(comp);
  }
  //printf("hello Apaar, No. of tasks = %d, My name = %d\n",numtasks, taskid);

  MPI_Finalize();
}
