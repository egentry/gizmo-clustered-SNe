#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char** argv) {
    // // // This all needs to be run from the scripts folder, to get the python correct
    // // // argv[1] should be $end_file, relative to the scripts folder
    // // // argv[2] should be $run_dir, relative to the scripts folder


    if(argc < 3)
    {
        return 0;
    }

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank==0)
    {

        int status = 0;
        char command[256] = "python prepare_for_restart.py ";
        strcat(command, argv[1]);
        strcat(command, " ");
        strcat(command, argv[2]);
        strcat(command, " ");
        printf("argv[3] = %s \n", argv[3]);
        if(argc > 3)
        {
            if (strcmp(argv[3],"--from-checkpoint") == 0)
            {
                strcat(command, argv[3]);
            }
        }

        printf("calling: '%s' \n", command);
        status = system(command);

    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Finalize the MPI environment.
    MPI_Finalize();
}