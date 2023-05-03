# 15618-Final-Project

website & description: https://tianqi-wu.github.io/course-projects/index.html
Dataset: https://archive.ics.uci.edu/ml/datasets/abalone
This is data of abalone from UCI.


How to run it: To run it with automated script, type the following in your console/terminal:

Please note that the two scripts, ./OpenMP_runner.sh and ./MPI_runner.sh, will need the users themselves to specify
by modifying the script themselves.

```
make
./OpenMP_runner.sh
./MPI_runner.sh
./combined_runner.sh <K-number> <max-iteration-number>
```

You can also run each of the executable with the following commands:

```
./knn <dataset> <K>
./k_means <dataset> <K> <maxIter>
mpirun --np <thread-num> ./knn_mpi <dataset> <K>
mpirun --np <thread-num> ./k_means_mpi <dataset> <K> <maxIter>
```

## Structure of the Repository

The repository consists of the final versions of our implementation.

.cpp, .h files: main files with our implementations

validator folder: notebook for validating correctness of algorithms

generator.py file: script for generating abalones. You need to change variables to change values.

data folder: our vanilla and testing data are included. You can use abalones.data or mass_abalones.data.

