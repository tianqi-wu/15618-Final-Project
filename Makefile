all: knn k_means knn_mpi k_means_mpi

knn: knn.cpp abalone.cpp abalone.h
	g++ -fopenmp -o knn knn.cpp -I.

k_means: k_means.cpp abalone.cpp abalone.h
	g++ -fopenmp -o k_means k_means.cpp -I.

knn_mpi: knn_mpi_main.cpp abalone.cpp abalone.h
	mpic++ knn_mpi_main.cpp -o knn_mpi

k_means_mpi: k_means_mpi_main.cpp abalone.cpp abalone.h
	mpic++ k_means_mpi_main.cpp -o k_means_mpi