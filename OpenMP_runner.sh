# The function for running everything together.
# The main reason for creating this function is that running each function
# takes tons of time, and to evaluate import values can be costly in human time.
# It runs a series of script for generating performance data for the whole process.
# Help from https://learnxinyminutes.com/docs/bash/
#!/usr/bin/bash
# set up global variables here
make
K=10
maxIter=100
dataset="./data/abalone.data"
echo "Start the runner process..."
# export OMP_NUM_THREADS=1
# TODO: update this array in PSC

echo "For OpenMP..."
echo "For OpenMP..."
echo "For OpenMP..."
echo "Start KNN..."
echo "Start KNN..."
echo "Start KNN..."
# export OMP_NUM_THREADS=$item
array=(1 2 4 8)
for item in "${array[@]}"; do
    export OMP_NUM_THREADS=$item
    echo "current number of threads: $item"
    ./knn $dataset $K
    echo ""
done

K=5
echo "Start K-means..."
echo "Start K-means..."
echo "Start K-means..."
#  mpirun --np 8 ./knn_mpi ./data/abalone.data 20
array=(1 2 4 8)
for item in "${array[@]}"; do
    export OMP_NUM_THREADS=$item
    echo "current number of threads: $item"
    ./k_means $dataset $K $maxIter
    echo ""
done
