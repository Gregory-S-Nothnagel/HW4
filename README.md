# MPI Implementation of Conway's Game of Life

This repository contains an implementation of Conway's Game of Life using C++ and the Message Passing Interface (MPI). There are three versions of this program:

- **HW4-V1.cpp**: Uses `MPI_Sendrecv`, blocking, no improvements.
- **HW4-V2.cpp**: Uses `MPI_Sendrecv`, blocking, with improvements (e.g., using pointers instead of vectors for faster swapping and condensing `MPI_Reduce` and `MPI_Scatter` into `MPI_Allreduce`).
- **HW4-V3.cpp**: Uses point-to-point communication, non-blocking, with improvements (e.g., using pointers instead of vectors for faster swapping and condensing `MPI_Reduce` and `MPI_Scatter` into `MPI_Allreduce`).

## Performance Testing

To measure performance, three rounds of tests were conducted on all three HW4 versions using `mpicxx`. Each round of each version produces its own output file, totaling nine output files in this repository. The final average times for each thread count use the times from each of the three rounds for that version.

### Test Setup

- **Board Size**: 5000 x 5000
- **Max Iterations**: 5000
- **Thread Counts**: 1, 2, 3, 4, 8, 10, 20  
  _(16 is not included due to issues with distributing uneven workloads, as mentioned in the report.)_

### Example Output

The following lines are from the output file `HW4-V1-ROUND_1.o203640`:

```
last iteration: 5000, time: 739.254
last iteration: 5000, time: 369.96
last iteration: 5000, time: 185.237
last iteration: 5000, time: 92.9542
last iteration: 5000, time: 74.4761
last iteration: 5000, time: 37.7634
```

Each output file lists times in order of threads, so:
- **1 thread**: `739.254`
- **2 threads**: `369.96`
- **4 threads**: `185.237`
- **8 threads**: `92.9542`
- **10 threads**: `74.4761`
- **20 threads**: `37.7634`

The numbers used in the report table are averages of the times from three separate output files, rounded to two decimal digits. The board size (5000) is not indicated in the output files but can be confirmed by checking the job submission script files.

### Output Files

Output files follow the naming format `HW4-VX-ROUND_Y.oZ`:
- **X**: Version of the assignment (HW4-V1.cpp / HW4-V2.cpp / HW4-V3.cpp)
- **Y**: Round number
- **Z**: Unique identifier

## Scripts

There are three job submission script files provided:
- **script_1.sh** for **HW4-V1.cpp**
- **script_2.sh** for **HW4-V2.cpp**
- **script_3.sh** for **HW4-V3.cpp**

The only difference between these scripts is the source file used for the job. All other features are identical. The script for each source file is run 3 times to get the 3 rounds for that source file.

## Compilation

To compile `HW4-V2.cpp` (the MPI_Sendrecv blocking, improved version), use the following script (contents of `script_1.sh`):

```bash
#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load openmpi

mpicxx -O3 -march=native -flto -funroll-loops -ffast-math -Wall -o HW4-V2 HW4-V2.cpp
```

### Running the Program

To run the compiled program, use the following command:

```bash
mpirun -n <NUM_PROCESSES> ./HW4-V2 <BOARD_SIZE> <MAX_ITERATIONS> <BOARD_OUTPUT_PATH>
```

#### Arguments

- **`<NUM_PROCESSES>`**: Number of processes to split the board iteration work into (e.g., `5`).
- **`<BOARD_SIZE>`**: Size of the board (e.g., `50` for a 50x50 board).  
  **Note**: `<BOARD_SIZE>` **must be a multiple of** `<NUM_PROCESSES>` **for proper functionality.** (see report, section "Program Design:", line beginning "*All versions (V1, V2, V3)...")
- **`<MAX_ITERATIONS>`**: Maximum number of iterations (e.g., `100`).
- **`<BOARD_OUTPUT_PATH>`**: Path to save the final board state.

#### Example Command

To test a small board (10x10) with a maximum of 1000 iterations:

```bash
mpirun -n 5 ./HW4-V2 10 1000 /scratch/$USER/test.txt
```
