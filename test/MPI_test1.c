// required MPI include file
#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
   int  numtasks, rank, len, rc;
   char hostname[MPI_MAX_PROCESSOR_NAME];

   // initialize MPI
   MPI_Init(&argc,&argv);

   // get number of tasks
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

   // get my rank
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   // this one is obvious
   MPI_Get_processor_name(hostname, &len);
   printf ("Number of tasks= %d My rank= %d Running on %s Cmd args %s %s Length %d\n", numtasks,rank,hostname,argv[0], argv[1], len);
   int a = 0;
   int b = 0;
   for(int i = 0; i < 1000; ++i) {
      a = i*(2+i)-i/(i+2);
      b = a*i/(i+1)*9;
   }

        // do some work with message passing


   // done with MPI
   MPI_Finalize();
}
