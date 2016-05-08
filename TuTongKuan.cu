include <stdio.h>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <deque>
#include <vector>
#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include <sstream>

#include <cuda_runtime.h>
#include <cstring>

#include "caffe/caffe.hpp"
#include "caffe/util/io.hpp"
#include "caffe/blob.hpp"

#include <tfs.h>
#include <ImageReader.h>
#include <CSVReader.h>

using namespace std;
using namespace caffe;

#define BUFSIZE 256
#define TAG 0
void MPI_init(int argc, char** argv, int &device_id, int& myid, int& numprocs)
{
  int devCount;
  char idstr[256];
  char idstr2[256];
  char buff[BUFSIZE];
  int i;
  int rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  /**
   * for Mvapich2
   */
//      rank = atoi(getenv("MV2_COMM_WORLD_RANK"));
  /**
   * for OpenMPI
   */
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  rank = atoi(getenv("OMPI_COMM_WORLD_LOCAL_RANK"));
  myid = rank;

  cudaGetDeviceCount(&devCount);
  device_id = myid % devCount;
  cudaSetDevice(device_id);
  printf("rank=%d, devCount=%d, device_id=%d\n", rank, devCount, device_id);
  
  PI_Status stat;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  myid = rank;
  MPI_Get_processor_name(processor_name, &namelen);

  if (myid == 0)
  {
    printf("  We have %d processors\n", numprocs);
    printf("  Spawning from %s \n", processor_name);
    printf("  CUDA MPI\n");
    printf("\n");
    for (i = 1; i < numprocs; i++)
    {
      buff[0] = 'I';
      MPI_Send(buff, BUFSIZE, MPI_CHAR, i, TAG, MPI_COMM_WORLD);
    }

    //cudaGetDeviceCount(&devCount);
    //device_id = myid % devCount;
    buff[1] = '\0';
    idstr[0] = '\0';
    if (devCount == 0)
    {
      sprintf(idstr, "- %-11s %5d %4d NONE", processor_name, rank, devCount);
    }
    else
    {
      if (devCount >= 1)
      {
        sprintf(idstr, "+ %-11s %5d %4d", processor_name, rank, devCount);
        idstr2[0] = '\0';
        //        for (int i = 0; i < devCount; ++i)
        {

          cudaDeviceProp devProp;
          cudaGetDeviceProperties(&devProp, device_id);
          sprintf(idstr2, " %s (%d) ", devProp.name, device_id);
          strncat(idstr, idstr2, BUFSIZE);
        }
      }
      else
      {
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, i);
        sprintf(idstr, "%-11s %5d %4d %s", processor_name, rank, devCount, devProp.name);
      }
    }
    strncat(buff, idstr, BUFSIZE);
    
    printf("  Probing nodes...\n");
    printf("     Node        Psid  CUDA Cards (devID)\n");
    printf("     ----------- ----- ---- ----------\n");

    printf("%s\n", buff);

    for (i = 1; i < numprocs; i++)
    {
      MPI_Recv(buff, BUFSIZE, MPI_CHAR, i, TAG, MPI_COMM_WORLD, &stat);
      printf("%s\n", buff);
    }
    printf("\n");
    //    MPI_Finalize();

  }
  else
  {
    MPI_Recv(buff, BUFSIZE, MPI_CHAR, 0, TAG, MPI_COMM_WORLD, &stat);
    MPI_Get_processor_name(processor_name, &namelen);
    buff[1] = '\0';
    idstr[0] = '\0';
    if (devCount == 0)
    {
      sprintf(idstr, "- %-11s %5d %4d NONE", processor_name, rank, devCount);
    }
    else
    {
      if (devCount >= 1)
      {
        sprintf(idstr, "+ %-11s %5d %4d", processor_name, rank, devCount);
        idstr2[0] = '\0';

        //        for (int i = 0; i < devCount; ++i)
        {
          cudaDeviceProp devProp;
          cudaGetDeviceProperties(&devProp, device_id);
          sprintf(idstr2, " %s (%d) ", devProp.name, device_id);
          strncat(idstr, idstr2, BUFSIZE);
        }
      }
      else
      {
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, device_id);
        sprintf(idstr, "%-11s %5d %4d %s", processor_name, rank, devCount, devProp.name);
      }
    }
    strncat(buff, idstr, BUFSIZE);
    MPI_Send(buff, BUFSIZE, MPI_CHAR, 0, TAG, MPI_COMM_WORLD);
  }

}

int main(int argc, char **argv)
{
  if (argc != 3)
  {
    LOG(ERROR) << "TuTongKuan TableName1 TableName2";
    return 1;
  }
  int device_id;
  int rank_id;
  int np;
  MPI_init(argc, argv, device_id, rank_id, np);
  char buf[128];
  FILE *pp;
  char cmd[] = "dship download ";

  if ((pp = popen(cmd, "r")) == NULL)
  {
    printf("popen() error!/n");
    exit(1);
  }

  while (fgets(buf, sizeof buf, pp))
  {
    printf("%s", buf);
  }
  pclose(pp);

  MPI_Finalize();
  return 0;
}