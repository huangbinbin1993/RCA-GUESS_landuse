#ifdef MPI_SRC
#include"mpi.h"
#endif
#include<stdlib.h>
#include<iostream>
#include"froutines.h"
using namespace std;

void terminate_program( int code )
{
#ifdef MPI_SRC
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );
  cout << "terminate_program called from proc " << myid  << "  code = "<< code<<endl;
  MPI_Abort( MPI_COMM_WORLD, code );
  MPI_Finalize();
#endif
  exit( code );
}



