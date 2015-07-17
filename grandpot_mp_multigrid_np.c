#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "header_files.h"
#include <math.h>
int main(int argc, char *argv[]) {

  initialize_grid();

  initialize_gradlayer();

  //Finished initialization of gradient layer
  init_propertymatrices(T);
  initialize_variables();

  //Finished Initializing
  printf("Lowest grid size: X= %d Time steps= %d\n",MESH_X,ntimesteps);
  printf("Initializing grid and writing initial.dat file...\n");

#ifdef MPI
  //Initializing MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  //Creating new data type
  Build_derived_type(&gridinfo_instance, &MPI_gridinfo);
  numworkers = numtasks-1;

  long shift_position=0;
  long to = atol(argv[1]);

  if (atol(argv[1]) !=0) {
#ifdef SHIFT
      long position,file_iter;
      long time;
      FILE *fp;
      fp = fopen("shift.dat","r");

      for(file_iter=0; file_iter <= atol(argv[1])/saveT; file_iter++) {
	fscanf(fp,"%ld %ld\n",&time, &position);
      }
      fclose(fp);
      shift_position = position;
#endif
   }
  if (taskid == MASTER) {
    gridinfo_levels = (struct variables** )malloc((REFINE_LEVELSY)*sizeof(**gridinfo_levels));

    for (levels=0; levels < (REFINE_LEVELSY); levels++) {
      gridinfo_levels[levels] = (struct variables* )malloc((grid_points[levels])*sizeof(*gridinfo_levels[levels]));
    }
    rows_worker   = (long *)malloc(numtasks*sizeof(*rows_worker));
    offset_worker = (long *)malloc(numtasks*sizeof(*offset_worker));

    if (to==0) {
      initdomain_struct(gridinfo_levels);
//       rdfrmfile_struct(gridinfo_levels,to);
    } else {
      rdfrmfile_struct(gridinfo_levels,to);
    }
#ifndef ISOTHERMAL
#ifdef TEMPGRADY
  temperature_gradientY.base_temp       = BASETEMP;
  temperature_gradientY.DeltaT          = DELTAT;
  temperature_gradientY.Distance        = DISTANCE;
  temperature_gradientY.gradient_OFFSET = OFFSET;
  temperature_gradientY.velocity        = VELOCITY;
  temperature_gradientY.GRADIENT        = temperature_gradientY.DeltaT/temperature_gradientY.Distance;
#ifdef SHIFT
  if (atol(argv[1]) !=0) {
    temperature_gradientY.gradient_OFFSET = (temperature_gradientY.gradient_OFFSET) + floor((temperature_gradientY.velocity)*(atol(argv[1])*deltat)) - shift_position*deltay;
  }
#else
  if (atol(argv[1]) !=0) {
    temperature_gradientY.gradient_OFFSET = (temperature_gradientY.gradient_OFFSET) + floor((temperature_gradientY.velocity)*(atol(argv[1])*deltat));
  }
#endif
#endif
#endif
//Starting of calculation scheme
//   distribute_levels();
    levels_exchange();
    apply_boundaryconditions();
    writetofile_struct(gridinfo_levels,to);
  }
  printf("taskid =%d\n", taskid);
  Mpiinfo(taskid);
#else
  if (to==0) {
    initdomain_struct(gridinfo_levels);
//     rdfrmfile_struct(gridinfo_levels,to);
  } else {
    rdfrmfile_struct(gridinfo_levels,to);
  }
  if (atol(argv[1]) !=0) {
#ifdef SHIFT
    long position,file_iter;
    long time;
    FILE *fp;
    fp = fopen("shift.dat","r");

    for(file_iter=0; file_iter <= atol(argv[1])/saveT; file_iter++) {
      fscanf(fp,"%ld %ld\n",&time, &position);
    }
    fclose(fp);
    shift_position = position;
#endif
  }
#ifndef ISOTHERMAL
#ifdef TEMPGRADY
  temperature_gradientY.base_temp       = BASETEMP;
  temperature_gradientY.DeltaT          = DELTAT;
  temperature_gradientY.Distance        = DISTANCE;
  temperature_gradientY.gradient_OFFSET = OFFSET;
  temperature_gradientY.velocity        = VELOCITY;
  temperature_gradientY.GRADIENT        = temperature_gradientY.DeltaT/temperature_gradientY.Distance;
#ifdef SHIFT
  if (atol(argv[1]) !=0) {
    temperature_gradientY.gradient_OFFSET = (temperature_gradientY.gradient_OFFSET) + floor((temperature_gradientY.velocity)*(atol(argv[1])*deltat)) - shift_position*deltay;
  }
#else
  if (atol(argv[1]) !=0) {
    temperature_gradientY.gradient_OFFSET = (temperature_gradientY.gradient_OFFSET) + floor((temperature_gradientY.velocity)*(atol(argv[1])*deltat));
  }
#endif
#endif
#endif
  levels_exchange();
  apply_boundaryconditions();
  writetofile_struct(gridinfo_levels,to);
#endif

#ifdef MPI
 if (taskid == MASTER) {
    FILE *fp;
    //Smoothing iterations
#ifdef TEMPGRADY
    apply_temperature_gradientY(gridinfo_levels, offset, rows, shift_OFFSET, t);
#endif
    for(t=1;t<=nsmooth;t++) {
      for (levels=0; levels < REFINE_LEVELSY; levels++) {
        if (levels ==0) {
          sendtoworker();
          receivefrmworker();
        } else {
          smooth(levels, 2, numx[levels]-2);
        }
      }
    }
    for(t=1;t<=ntimesteps;t++) {
#ifdef TEMPGRADY
      apply_temperature_gradientY(gridinfo_levels, offset, rows, shift_OFFSET, t);
#endif
      printf("Iteration = %d\n", t, taskid);
      for (levels=0; levels < REFINE_LEVELSY; levels++) {
	if (levels ==0) {
	  sendtoworker();
	  receivefrmworker();
	} else {
	  solverloop(levels, 2, numx[levels]-2);
	}
      }
      levels_exchange();
      apply_boundaryconditions();
      if(t%saveT==0) {
	/* Write final output, call X graph and finalize MPI */
	printf("Writing final.dat file and generating graph...\n");
  //       assemble_layers();
	writetofile_struct(gridinfo_levels,to+t);
	fp=fopen("shift.dat","a");
	fprintf(fp,"%ld %ld\n",t+to, shift_OFFSET+shift_position);
	fclose(fp);
      }
//Apply post-condition
#ifdef SHIFT
      for (rank=1; rank <=numworkers; rank++) {
	source  = rank;
	msgtype = SHIFT_SIGNAL;

	MPI_Recv(&SHIFT_Y,       1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
	msgtype = SHIFT_POS;
	MPI_Recv(&INTERFACE_POS, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);

	if (SHIFT_Y ==1) {
	  shift_ON=1;
	}
	if (INTERFACE_POS > MAX_INTERFACE_POS) {
	  MAX_INTERFACE_POS = INTERFACE_POS;
	}
      }
      //Calculate the global shift position here!!...................................................
      if (shift_ON ==1) {
	printf("shift in the negative y direction=%ld\n",shift_gridpoints[0]);
	apply_shiftY(gridinfo_levels, MAX_INTERFACE_POS, 0, 0);
	levels_exchange();
	apply_boundaryconditions();
  //       shift_OFFSET += (INTERFACE_POS_GLOBAL-shiftj);
	shift_OFFSET += shift_gridpoints[0];
      }
      shift_ON = 0;
#endif
    }
    Free_variables();
    MPI_Type_free(&MPI_gridinfo);
    MPI_Finalize();
  }

  if (taskid != MASTER) {
    printf("start=%d,end=%d\n", start, end);
    //Smoothing iterations
#ifndef ISOTHERMAL
#ifdef TEMPGRADY
  temperature_gradientY.base_temp       = BASETEMP;
  temperature_gradientY.DeltaT          = DELTAT;
  temperature_gradientY.Distance        = DISTANCE;
  temperature_gradientY.gradient_OFFSET = OFFSET;
  temperature_gradientY.velocity        = VELOCITY;
  temperature_gradientY.GRADIENT        = DELTAT/DISTANCE;
#endif
#endif
    for(t=1;t<=nsmooth;t++) {
      receivefrmmaster(taskid);
      mpiexchange(taskid);
      smooth(0, start, end);
      mpiboundary(taskid);
      sendtomaster();
    }
    for(t=1;t<=ntimesteps;t++) {
      receivefrmmaster(taskid);
      mpiexchange(taskid);
      solverloop(0, start, end);
      mpiboundary(taskid);
      sendtomaster();
#ifdef SHIFT
      if (MAX_INTERFACE_POS > shift_MAX) {
	shift_ON = 1;
      }
      MPI_Send(&shift_ON, 1, MPI_INT, MASTER, SHIFT_SIGNAL, MPI_COMM_WORLD);
      MPI_Send(&MAX_INTERFACE_POS, 1, MPI_INT, MASTER, SHIFT_POS, MPI_COMM_WORLD);
      shift_ON = 0;
#endif
    }
    free(gridinfo_levels[0]);
    for(levels=0; levels < REFINE_LEVELSY; levels++) {
      for (layer=0; layer < NUM_BUFFER_LAYERS; layer++) {
	free(gradient1_levels[levels][layer]);
      }
      free(gradient1_levels[levels]);
    }
    free(gradient1_levels);
    MPI_Type_free(&MPI_gridinfo);
    MPI_Finalize();
  }
  return(0);
#else
  FILE *fp;
  for(t=1;t<=ntimesteps;t++) {
    printf("Iteration = %d\n", t);
    for (levels=0; levels < REFINE_LEVELSY; levels++) {
      solverloop(levels, 2, numx[levels]-2);
    }
    levels_exchange();
    apply_boundaryconditions();
    if(t%saveT==0) {
      /* Write final output, call X graph and finalize MPI */
      printf("Writing final.dat file and generating graph...\n");
//       assemble_layers();
      writetofile_struct(gridinfo_levels,to+t);
      fp=fopen("shift.dat","a");
      fprintf(fp,"%ld %ld\n",t+to, shift_OFFSET+shift_position);
      fclose(fp);
    }
    //Apply post-condition
#ifdef SHIFT
    if (shift_ON ==1) {
      printf("shift in the negative y direction=%ld\n",shift_gridpoints[0]);
      apply_shiftY(gridinfo_levels, MAX_INTERFACE_POS, 0, 0);
      levels_exchange();
      apply_boundaryconditions();
//       shift_OFFSET += (INTERFACE_POS_GLOBAL-shiftj);
      shift_OFFSET += shift_gridpoints[0];
    }
    shift_ON = 0;
#endif
  }
  Free_variables();
  return(0);
#endif
}
// void writetofile_struct_worker(struct variables** gridinfo,int t, int start, int end,int taskid) {
//    long x, y, gidy;
//    long b;
//    char name[1000];
//    FILE *fp;
//    for (b=0; b < NUMPHASES; b++) {
//     sprintf(name,"%s_%d_worker_%d.dat",Phases[b], t, taskid);
//     fp=fopen(name,"w");
//     for(x=start-2; x <=(end+2); x++) {
//       for(y=0; y < numy[0]; y++) {
// 	gidy = x*numy[0] + y;
// 	fprintf(fp, "%le %le %le \n", x*geometry[0].DeltaX, y*geometry[0].DeltaY, gridinfo[0][gidy].phia[b]);
//       }
//       fprintf(fp,"\n");
//     }
//     fclose(fp);
//   }
// }
