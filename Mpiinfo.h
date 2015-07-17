void Mpiinfo(long taskid) {
  int rank;
  if (taskid == MASTER) {
    averow = (numx[0])/numworkers;
    extra  = (numx[0])%numworkers;
    offset = 0;
    for (rank=1; rank<=numworkers; rank++) {
	rows = (rank <= extra) ? averow+1 : averow;
	/* Tell each worker who its neighbors are, since they must exchange */
	/* data with each other. */  
	if (rank == 1)
  #ifdef PERIODIC
	  left_node = numworkers;
  #else
	  left_node = NONE;
  #endif
	else
	  left_node = rank - 1;
	      
	if (rank == numworkers)
  #ifdef PERIODIC
	  right_node = 1;
  #else
	  right_node = NONE;
  #endif
	else
	  right_node = rank + 1;
      /*  Now send startup information to each worker  */
      dest = rank;
      MPI_Send(&offset,                                      1, MPI_INT,      dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&rows,                                        1, MPI_INT,      dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&left_node,                                   1, MPI_INT,      dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(&right_node,                                  1, MPI_INT,      dest, BEGIN, MPI_COMM_WORLD);
      MPI_Send(gridinfo_levels[0]+offset*numy[0], rows*numy[0], MPI_gridinfo, dest, BEGIN, MPI_COMM_WORLD);
      printf("Sent to task %d: rows= %d offset= %d ",dest,rows,offset);
      printf("left= %d right= %d, rows=%d\n",left_node,right_node,rows);
      offset_worker[rank] = offset;
      rows_worker[rank]   = rows;
      offset              = offset + rows;
    }
  } else {
    source  = MASTER;
    msgtype = BEGIN;
    MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, offset=%d\n",taskid,offset);
    MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, taskid=%ld, rows*(layer_size)=%ld\n",taskid,rows*(numy[0]));
    MPI_Recv(&left_node,  1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, left=%d\n",left_node);
    MPI_Recv(&right_node, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
    printf("recieved, right=%d\n",right_node);
    gridinfo_levels   = (struct variables** )malloc((1)*sizeof(**gridinfo_levels));
    if ((offset==0) || (offset+rows == numx[0])) {
      if (offset == 0) {
	start = 2;
	end = rows-1;
	gridinfo_levels[0] = (struct variables *)malloc((((rows)+2)*numy[0])*(sizeof(*gridinfo_levels[0])));
	MPI_Recv(gridinfo_levels[0],rows*numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      } else {
	start = 2;
	end = rows-1;
	gridinfo_levels[0] = (struct variables*)malloc((((rows)+2)*numy[0])*(sizeof(*gridinfo_levels[0])));
	MPI_Recv(gridinfo_levels[0]+2*numy[0],rows*numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      }
    } else {
      start = 2;
      end  = rows+1;
      gridinfo_levels[0] = (struct variables*)malloc(((rows+4)*numy[0])*(sizeof(*gridinfo_levels[0])));
      MPI_Recv(gridinfo_levels[0]+2*numy[0],rows*numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
//     writetofile_struct_worker(gridinfo_levels,0,start,end,taskid);
  }
}
void receivefrmmaster(long taskid) {
  source  = MASTER;
  msgtype = BEGIN;
  if ((offset==0) || (offset+rows == numx[0])) {
    if (offset == 0) {
      MPI_Recv(gridinfo_levels[0],rows*numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    } else {
      MPI_Recv(gridinfo_levels[0]+2*numy[0],rows*numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
  } else {
    MPI_Recv(gridinfo_levels[0]+2*numy[0],rows*numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
  }
//   writetofile_struct_worker(gridinfo_levels,1,start,end,taskid);
}
void sendtoworker() {
  int rank;
  for (rank=1; rank<=numworkers; rank++) {
    dest = rank;
    MPI_Send(gridinfo_levels[0]+offset_worker[rank]*numy[0], rows_worker[rank]*numy[0], MPI_gridinfo, dest, BEGIN, MPI_COMM_WORLD);
  }
}
void sendtomaster() {
  MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
  MPI_Send(&rows,   1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
  if(offset==0) {
    MPI_Send(gridinfo_levels[0],           rows*numy[0], MPI_gridinfo, MASTER, DONE, MPI_COMM_WORLD);
  } else {
    MPI_Send(gridinfo_levels[0]+2*numy[0], rows*numy[0], MPI_gridinfo, MASTER, DONE, MPI_COMM_WORLD);
  }
}
void receivefrmworker() {
  int rank;
  for (rank=1; rank <=numworkers; rank++) {
    source  = rank;
    msgtype = DONE;
    MPI_Recv(&offset,                                         1, MPI_INT,      source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,                                           1, MPI_INT,      source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Recv(&gridinfo_levels[0][offset*numy[0]],  rows*numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
  }
}
void mpiexchange(long taskid) {
  if (taskid %2) {
    if (taskid != 1) {
      MPI_Send(&gridinfo_levels[0][2*numy[0]], numy[0], MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
      MPI_Send(&gridinfo_levels[0][3*numy[0]], numy[0], MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
    
      source = left_node;
      msgtype = LTAG;
      MPI_Recv(&gridinfo_levels[0][1*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&gridinfo_levels[0][0*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
    if (taskid != numworkers) {
      MPI_Send(&gridinfo_levels[0][end*numy[0]], numy[0], MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
      MPI_Send(&gridinfo_levels[0][(end-1)*numy[0]], numy[0], MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
      
      source = right_node;
      msgtype = RTAG;
      MPI_Recv(&gridinfo_levels[0][(end+1)*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&gridinfo_levels[0][(end+2)*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    }
  } else {
    if (taskid != numworkers) {
      source = right_node;
      msgtype = RTAG;
      MPI_Recv(&gridinfo_levels[0][(end+1)*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&gridinfo_levels[0][(end+2)*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
	
      MPI_Send(&gridinfo_levels[0][end*numy[0]], numy[0], MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
      MPI_Send(&gridinfo_levels[0][(end-1)*numy[0]], numy[0], MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
    }
    if(taskid !=1) {
      source = left_node;
      msgtype = LTAG;
      MPI_Recv(&gridinfo_levels[0][1*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&gridinfo_levels[0][0*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
      
      MPI_Send(&gridinfo_levels[0][2*numy[0]], numy[0], MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
      MPI_Send(&gridinfo_levels[0][3*numy[0]], numy[0], MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
    }
  }
}
void mpiboundary(long taskid) {
#ifdef PERIODIC
  if(taskid==1) {
    MPI_Send(&gridinfo_levels[0][2*numy[0]], numy[0], MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
    source = left_node;
    msgtype = LTAG;
    MPI_Recv(&gridinfo_levels[0][1*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Send(&gridinfo_levels[0][3*numy[0]], numy[0], MPI_gridinfo, left_node, RTAG, MPI_COMM_WORLD);
    source = left_node;
    msgtype = LTAG;
    MPI_Recv(&gridinfo_levels[0][0*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
  }
  if(taskid==numworkers) {
    source = right_node;
    msgtype = RTAG;
    MPI_Recv(&gridinfo_levels[0][(end+1)*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Send(&gridinfo_levels[0][end*numy[0]], numy[0], MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
    
    source = right_node;
    msgtype = RTAG;
    MPI_Recv(&gridinfo_levels[0][(end+2)*numy[0]], numy[0], MPI_gridinfo, source, msgtype, MPI_COMM_WORLD, &status);
    MPI_Send(&gridinfo_levels[0][(end-1)*numy[0]], numy[0], MPI_gridinfo, right_node, LTAG, MPI_COMM_WORLD);
  }
#endif
}