HIPMCL_EXE=/hipmcl/bin/hipmcl
IN_FILE=/data/network.out
OUT_FILE=/data/network.out.hipmcl

export OMP_NUM_THREADS=3
mpirun -np 1 $HIPMCL_EXE -M $IN_FILE -I 2 -per-process-mem 12 -o $OUT_FILE
