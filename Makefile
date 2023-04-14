all: knn k_means

knn: knn.cpp abalone.cpp abalone.h
	g++ -fopenmp -o knn knn.cpp -I.

k_means: k_means.cpp abalone.cpp abalone.h
	g++ -fopenmp -o k_means k_means.cpp -I.
