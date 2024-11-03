#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load openmpi

mpicxx -O3 -march=native -flto -funroll-loops -ffast-math -Wall -o HW4-V1 HW4-V1.cpp
mpirun -n 1 ./HW4-V1 5000 5000 /scratch/$USER/HW4-V1-board_1-thread.txt
mpirun -n 2 ./HW4-V1 5000 5000 /scratch/$USER/HW4-V1-board_2-thread.txt
mpirun -n 4 ./HW4-V1 5000 5000 /scratch/$USER/HW4-V1-board_4-thread.txt
mpirun -n 8 ./HW4-V1 5000 5000 /scratch/$USER/HW4-V1-board_8-thread.txt
mpirun -n 10 ./HW4-V1 5000 5000 /scratch/$USER/HW4-V1-board_10-thread.txt
mpirun -n 20 ./HW4-V1 5000 5000 /scratch/$USER/HW4-V1-board_20-thread.txt

diff /scratch/$USER/HW4-V1-board_1-thread.txt /scratch/$USER/HW4-V1-board_2-thread.txt
diff /scratch/$USER/HW4-V1-board_1-thread.txt /scratch/$USER/HW4-V1-board_4-thread.txt
diff /scratch/$USER/HW4-V1-board_1-thread.txt /scratch/$USER/HW4-V1-board_8-thread.txt
diff /scratch/$USER/HW4-V1-board_1-thread.txt /scratch/$USER/HW4-V1-board_10-thread.txt
diff /scratch/$USER/HW4-V1-board_1-thread.txt /scratch/$USER/HW4-V1-board_20-thread.txt