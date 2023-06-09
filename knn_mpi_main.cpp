/**
 * KNN based on UCI abalone data.
 * we have specified the input correctluy so
 * https://github.com/ispc/ispc/tree/v1.15.0/examples/simple
 * mpirun -np 8 knn_mpi
 */

using namespace std;
#include "abalone.h"
#include "mpi.h"
#include "common.h"
#include "timing.h"
#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <queue>

#define MASTER 0

typedef pair<double, int> abaloneKeyValue;
/**
 * Sequential version of KNN. Uses PriorityQueue to help.
 * An abstraction of what a single thread will do.
 * Will extend to OpenMP friendly version in the future.
 */
double KNN_sequential(vector<Abalone> data, int K, Abalone someAbalone)
{
    priority_queue<abaloneKeyValue, vector<abaloneKeyValue>, greater<abaloneKeyValue>> pq;
    // populate the priority queue

    for (int i = 0; i < data.size(); i++)
    {
        double differenceValue = calculateDistanceEuclidean(data[i], someAbalone, true);
        int ringNumber = data[i].rings;
        pq.push(make_pair(differenceValue, ringNumber));
        
    }
    double sum = 0;
    for (int i = 0; i < K; i++)
    {
        pair<double, int> top = pq.top();
        sum += top.second;
        pq.pop();
    }
    return sum / K;
}

/**
 * Separates data in a 70%/30% fashion.
 * This is guaranteed to be deterministic.
 */

void pesudo_training_test_parse(vector<Abalone> &training, vector<Abalone> &testing)
{
    vector<Abalone> updatedTrainingData;
    for (int i = 0; i < training.size(); i++)
    {
        if (i % 8 == 0 || i % 6 == 0)
        {
            testing.push_back(training[i]);
        }
        else
        {
            updatedTrainingData.push_back(training[i]);
        }
    }
    training = updatedTrainingData;
}

/**
 * Driver function.
 * Used to run the thing correctly
 */
int main(int argc, char *argv[])
{
    int pid;   // current PID
    int nproc; // process rank
    int offset;
    int source;
    int dest;
    MPI_Status status;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    // Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // assign tags to the values

    int tag4 = 1;
    int tag3 = 2;
    int tag2 = 3;
    int tag1 = 4;

  string location = "";
  int K = 20;
    if(argc <= 1) {
      //printf("Warning: no location or K specified: will run the default version.\n");
      location = "./data/abalone.data";
      K = 20;
    }else if(argc == 3){
      location = argv[1];
      K = atoi(argv[2]);
      if(K == 0) {
         printf("usage:mpirun --np <num_thread> %s <filename> <K-num>\n", argv[0]);
         return 2;
      }
    }else {
      printf("usage: mpirun --np <num_thread> %s <filename> <K-num>\n", argv[0]);
      return 2;
    }
    // https://stackoverflow.com/questions/37532631/read-class-objects-from-file-c
    ifstream fin;
    fin.open(location);
    if (!fin)
    {
        cerr << "Error in opening the file!" << endl;
        return 1; // if this is main
    }

    // reading information into the correct place.
    vector<Abalone> abalones;
    vector<float> abaloneResult; 
    vector<float> tempAbaloneResultStorage;

    vector<Abalone> training;
    vector<Abalone> testing;

    if (pid == 0)
    {
        Abalone temp;
        while (fin >> temp.sex >> temp.length >> temp.diameter >>
               temp.height >> temp.whole_height >> temp.shucked_weight >> temp.viscera_weight >> temp.shell_weight >> temp.rings)
        {
            abalones.push_back(temp);
        }
        training = abalones;
        pesudo_training_test_parse(training, testing);
        abaloneResult.resize(testing.size());
    }


    Timer totalSimulationTimer;
    int trainingSize = training.size();
    int testingSize = testing.size();

    int chunksize = (testingSize / nproc);

    int leftover = (testingSize % nproc);

    MPI_Bcast(&trainingSize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&testingSize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    if (pid == 0)
    {
        offset = chunksize + leftover;

        for (dest = 1; dest < nproc; dest++)
        {
            MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            offset = offset + chunksize;
        }

        int chunksize = (testingSize / nproc);
        int leftover = (testingSize % nproc);
        tempAbaloneResultStorage.resize(chunksize + leftover);
    }
    else
    {
        MPI_Recv(&offset, 1, MPI_INT, MASTER, tag1, MPI_COMM_WORLD, &status);
        int chunksize = (testingSize / nproc);
        int leftover = (testingSize % nproc);
        tempAbaloneResultStorage.resize(chunksize);
    }

    training.resize(trainingSize);
    testing.resize(testingSize);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(training.data(), sizeof(Abalone) * trainingSize, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(testing.data(), sizeof(Abalone) * testingSize, MPI_CHAR, MASTER, MPI_COMM_WORLD);

    if (pid == MASTER)
    {

        offset = 0;
        // START of update. range: (offset, chunksize+leftover, taskid);
        std::vector<Abalone> subAbalones =
            {testing.begin() + offset, testing.begin() + offset + chunksize + leftover};
        
        

        for (int i = 0; i < subAbalones.size(); i++)
        {
            tempAbaloneResultStorage[i] = KNN_sequential(training, K, subAbalones[i]);
        }

        // overwrite Storage to perform update for the master itself
        std::copy(std::begin(tempAbaloneResultStorage), std::end(tempAbaloneResultStorage), std::begin(abaloneResult) + offset);

        // Wait to receive results from each task

        for (int i = 1; i < nproc; i++)
        {
            source = i;
            offset = leftover + chunksize * i;
            MPI_Recv(&abaloneResult[offset], chunksize * sizeof(float), MPI_CHAR, source, tag2, MPI_COMM_WORLD, &status);
        }
    }
    /* end of master section */

    /***** Non-master tasks only *****/

    if (pid > MASTER)
    {
        /* Receive my portion of array from the master task */
        source = MASTER;
        chunksize = (testingSize / nproc);

        // START of update. range: (offset, chunksize+leftover, taskid);

        std::vector<Abalone> subAbalones =
            {testing.begin() + offset, testing.begin() + offset + chunksize};

        // KNN!!!!!
        for (int i = 0; i < subAbalones.size(); i++)
        {
            tempAbaloneResultStorage[i] = KNN_sequential(training, K, subAbalones[i]);
        }

        // END of update

        dest = MASTER;

        MPI_Send(tempAbaloneResultStorage.data(), sizeof(float) * tempAbaloneResultStorage.size(), MPI_CHAR, MASTER, tag2, MPI_COMM_WORLD);
    }

   

    MPI_Barrier(MPI_COMM_WORLD);
    double totalSimulationTime = totalSimulationTimer.elapsed();

    if (pid == 0)
    {
        printf("total simulation time: %.6fs\n", totalSimulationTime);
        saveToFile("output.txt", abaloneResult);
    }

    MPI_Finalize();
    return 0;
}
