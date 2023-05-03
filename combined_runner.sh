# The function for running everything together.
# The main reason for creating this function is that running each function
# takes tons of time, and to evaluate import values can be costly in human time.
# It runs a series of script for generating performance data for the whole process.
# Help from https://learnxinyminutes.com/docs/bash/
#!/usr/bin/bash
# set up global variables here
export I_MPI_PIN_DOMAIN=omp


make
K=$1
maxIter=$2
dataset="./data/custom_abalone.data"
echo "Start the runner process..."
# export OMP_NUM_THREADS=1
# TODO: update this array in PSC
echo "Using Mass Abalone Data (54000 abalones)"
echo "For combined..."
echo "For combined..."
echo "For combined..."
echo "Start KNN Combined..."
echo "Start KNN Combined..."
echo "Start KNN Combined..."
# export OMP_NUM_THREADS=$item
array=(2 4 8)
array1=(4 8 16)
for item in "${array[@]}"; do
    export OMP_NUM_THREADS=$item
    echo "current number of threads per MPI process : $item"
    for item1 in "${array1[@]}"; do
	echo "current number of MPI processes: $item1"
	mpirun --bind-to none -np $item1 ./knn_combined_v2 $dataset $K
	echo ""

    done
done

K=$1
max_iter=$2
echo "Start K-means Combined..."
echo "Start K-means Combined..."
echo "Start K-means Combined..."
#  mpirun --np 8 ./knn_mpi ./data/abalone.data 20
array=(2 4 8)
array1=(4 8 16)
for item in "${array[@]}"; do
    export OMP_NUM_THREADS=$item
    echo "current number of threads per MPI process : $item"
    for item1 in "${array1[@]}"; do
	echo "current number of MPI processes: $item1"
	mpirun --bind-to none -np $item1 ./k_means_combined_v2 $dataset $K $max_iter
	echo ""
    done
done
