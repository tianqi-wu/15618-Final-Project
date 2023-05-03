PROGNAME = knn k_means knn_mpi k_means_mpi k_means_combined knn_combined knn_combined_v2 k_means_combined_v2


all: $(PROGNAME)

knn: knn.cpp abalone.cpp abalone.h
	g++ -fopenmp -o knn knn.cpp -I.

k_means: k_means.cpp abalone.cpp abalone.h
	g++ -fopenmp -o k_means k_means.cpp -I.

knn_mpi: knn_mpi_main.cpp abalone.cpp abalone.h
	mpic++ knn_mpi_main.cpp -o knn_mpi

k_means_mpi: k_means_mpi_main.cpp abalone.cpp abalone.h
	mpic++ k_means_mpi_main.cpp -o k_means_mpi

knn_combined: knn_combined.cpp abalone.cpp abalone.h
	mpic++ -fopenmp -o knn_combined knn_combined.cpp

k_means_combined: k_means_combined.cpp abalone.cpp abalone.h
	mpic++ -fopenmp -o k_means_combined k_means_combined.cpp
knn_combined_v2: knn_combined_v2.cpp abalone.cpp abalone.h
	mpic++ -fopenmp -o knn_combined_v2 knn_combined_v2.cpp

k_means_combined_v2: k_means_combined_v2.cpp abalone.cpp abalone.h
	mpic++ -fopenmp -o k_means_combined_v2 k_means_combined_v2.cpp

clean:
	rm $(PROGNAME)
