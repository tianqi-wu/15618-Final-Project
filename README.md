# 15618-Final-Project

website & description: https://tianqi-wu.github.io/course-projects/index.html
Dataset: https://archive.ics.uci.edu/ml/datasets/abalone
This is data of abalone from UCI.


How to run it: type the following in your console/terminal:



```
make
```


You can also try running 

```
./OpenMP_runner.sh 
./MPI_runner.sh 
```

for inspecting simple results. 

## Advanced

```
The requests are generated via command line: we have provided usage messages. For KNN/K-means in OpenMP, please use formats such as ./knn <filename> <K-num> to specify the data. For K-means, please use mpirun --np <num_thread> ./k_means <filename> <K-num> <maxIter-num>. You can also look at the script files ending in .sh files for more information.
```
