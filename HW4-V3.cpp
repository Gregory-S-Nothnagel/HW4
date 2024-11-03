/*

Gregory Nothnagel
gsnothnagel@crimson.ua.edu 
CS 581
Homework #4

To Compile, run the following script:

	#!/bin/bash
	source /apps/profiles/modules_asax.sh.dyn
	module load openmpi

	mpicxx -O3 -march=native -flto -funroll-loops -ffast-math -Wall -o HW4-V3 HW4-V3.cpp
 
To Run, use the following command:

	mpirun -n <NUM_PROCESSES> ./HW4-V3 <BOARD_SIZE> <MAX_ITERATIONS> <BOARD_OUTPUT_PATH>
	
For more details, see "Compilation" section of README.md in the github repo

* Then you can inspect the board output file to confirm that everything looks as expected

*/

#include <mpi.h>
#include <iostream>

using namespace std;

// Function to print the grid on process 0
void printGrid(FILE* f, char* grid, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
			if (grid[i * N + j] == 1) fprintf(f, "█");
			else fprintf(f, "-");
        }
        fprintf(f, "\n");
    }
}

bool changeDetected;
// Function to apply Conway's Game of Life rules to a given section of the grid
void applyConwayRules(char* subgrid, char* nextSubgrid, int M, int N) {
    for (int i = 1; i <= M; ++i) {
        for (int j = 0; j < N; ++j) {
            int neighbors = subgrid[(i - 1) * N + j] +
                            (j > 0) * subgrid[i * N + (j - 1 + N) % N] +
                            (j < N - 1) * subgrid[i * N + (j + 1) % N] +
                            subgrid[(i + 1) * N + j] +
                            (j > 0) * subgrid[(i - 1) * N + (j - 1 + N) % N] +
                            (j < N - 1) * subgrid[(i - 1) * N + (j + 1) % N] +
                            (j > 0) * subgrid[(i + 1) * N + (j - 1 + N) % N] +
                            (j < N - 1) * subgrid[(i + 1) * N + (j + 1) % N];

            // Conway's Game of Life rules
            nextSubgrid[i * N + j] = (subgrid[i * N + j] && (neighbors == 2 || neighbors == 3)) ||
                                     (!subgrid[i * N + j] && neighbors == 3);
			
			changeDetected = changeDetected || (nextSubgrid[i * N + j] != subgrid[i * N + j]);
			
			//nextSubgrid[i * N + j] = neighbors;
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	
	// Check if the correct number of arguments is provided
    if (argc < 4) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <board_size> <max_iterations> <output_directory>" << endl;
        }
        MPI_Finalize();
        return 1; // Exit with error
    }

	int MAX_ITERATIONS = std::atoi(argv[2]); 	// Maximum number of iterations
    const char* outputDir = argv[3];      		// Output directory
    int N = std::atoi(argv[1]);     			// Board size
    int M = N / num_procs;  					// Number of rows each process will handle

    char* fullGrid = nullptr;

    if (rank == 0) {
        // Process 0 initializes the full grid
        fullGrid = new char[N * N];
        for (int i = 0; i < N * N; i++) {
            fullGrid[i] = rand() % 2;
        }
    }
	
	double startTime = MPI_Wtime();


    // Each process allocates memory for (M+2)xN subgrid (M rows + 2 buffer rows)
    char* subgrid = new char[(M + 2) * N];
	char* nextSubgrid = new char[(M + 2) * N];
	for (int i = 0; i < (M + 2) * N; i++) {
		subgrid[i] = false;
		nextSubgrid[i] = false;
	}

    // Scatter rows to each process, sending M rows to each process
    MPI_Scatter(fullGrid, M * N, MPI_CHAR, subgrid + N, M * N, MPI_CHAR, 0, MPI_COMM_WORLD);

	bool globalChangeDetected = true;
	int iteration;
	for (iteration = 0; iteration < MAX_ITERATIONS && globalChangeDetected; iteration++) {
		changeDetected = false;
		globalChangeDetected = false;

		MPI_Request sendRequest[2];
		MPI_Request recvRequest[2];
		int req_count = 0;

		// Initialize non-blocking communication for the top and bottom rows
		if (rank > 0) {
			MPI_Irecv(subgrid, N, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, &recvRequest[0]);       // Receive from above
			MPI_Isend(subgrid + N, N, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, &sendRequest[0]);   // Send to above
			MPI_Wait(&recvRequest[0], MPI_STATUS_IGNORE);  // Wait for top border to arrive
			MPI_Wait(&sendRequest[0], MPI_STATUS_IGNORE);  // Complete send to above
		}
		if (rank < num_procs - 1) {
			MPI_Irecv(subgrid + (M + 1) * N, N, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD, &recvRequest[1]);  // Receive from below
			MPI_Isend(subgrid + M * N, N, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD, &sendRequest[1]);        // Send to below
			MPI_Wait(&recvRequest[1], MPI_STATUS_IGNORE);  // Wait for bottom border to arrive
			MPI_Wait(&sendRequest[1], MPI_STATUS_IGNORE);  // Complete send to below
		}
		
		// Apply Conway’s rules and swap grids and reduce to detect changes
		applyConwayRules(subgrid, nextSubgrid, M, N);
		swap(subgrid, nextSubgrid);

		// Reduce to detect changes
		MPI_Allreduce(&changeDetected, &globalChangeDetected, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
	}



	// Gather all rows from each process back into fullGrid
	MPI_Gather(subgrid + N, M * N, MPI_CHAR, fullGrid, M * N, MPI_CHAR, 0, MPI_COMM_WORLD);

	double endTime = MPI_Wtime();

    if (rank == 0) {
        // Process 0 prints the combined grid
		FILE* f = fopen(outputDir, "w");
		cout << "last iteration: " << iteration << ", time: " << endTime - startTime << endl;
        printGrid(f, fullGrid, N);
		fclose(f);
        delete[] fullGrid;
    }

    MPI_Finalize();
    return 0;
}
