/**
 * KNN based on UCI abalone data.
 * we have specified the input correctluy so
 */

using namespace std;
#include "abalone.h"
#include "mpi.h"
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

/* Pesudo-random initialization.
 * For testing purposes. We want this algorithm to be deterministic.
 */

vector<Abalone> intializeRandomPesudo(vector<Abalone> data, int K)
{
    vector<Abalone> newVector;
    int slice = data.size() / K;
    for (int i = 0; i < K; i++)
    {
        int rand_num = slice * i;
        // copy construct. They are no longer what they were any more.
        Abalone newAbalone = Abalone(data[rand_num]);
        newVector.push_back(newAbalone);
        // std::cout << newAbalone.show() << std::endl;
    }

    return newVector;
}

/**
 * Random initialization. So-called one of the worst in K-means.
 * Implemented this only to make the algorithm work. Will extend when I have time.
 * We don't need to get the id. It's not necessary at all.
 */
vector<Abalone> intializeRandom(vector<Abalone> data, int K)
{
    vector<Abalone> newVector;

    for (int i = 0; i < K; i++)
    {
        int rand_num = rand() % data.size();
        // copy construct. They are no longer what they were any more.
        newVector.push_back(Abalone(data[rand_num]));
    }

    return newVector;
}

/**
 * KNN furthest point herustic.
*/

vector<Abalone> intializeFurthestPointHerustic(vector<Abalone> data, int K) {
    vector<Abalone> newVector;
    int lastAbaloneCluster = 0;
    for(int i = 0; i < K; i++) {
        int currFarthestPointIndex = 0;
        if(i != 0) {
            double maxSum = 0;
            for(int j = 0; j < data.size(); j++) {
                double currentSum = 0;
                for(int k = 0; k < i; k++) {
                    currentSum += calculateDistanceEuclidean(data[j], newVector[k], false);
                }
                if(maxSum < currentSum) {
                    maxSum = currentSum;
                    currFarthestPointIndex = j;
                }                
            }
        }
        // copy construct. They are no longer what they were any more.
        Abalone newAbalone = Abalone(data[currFarthestPointIndex]);
        newVector.push_back(newAbalone);
        //std::cout << newAbalone.show() << std::endl; 
    }

    return newVector;
}




/**
 * Sequential version of K-means. Uses PriorityQueue to help.
 * An abstraction of what a single thread will do.
 * Will extend to OpenMP friendly version in the future.
 */

/**
 * Driver function.
 * Used to run the thing correctly
 */
int main(int argc, char **argv)
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
    if(argc <= 1) {
      location = "./data/abalone.data";
    }else if(argc == 2){
      location = argv[1];
    }else {
      printf("usage: %s <filename>", argv[0]);
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

    // pass it in the K-means hyperparameter function
    int K = 5;
    int maxIter = 100;

    // reading information into the correct place.

    vector<Abalone> abalones;
    vector<int> clusterAssignment; // the final result for cluster assignment
    // resize so it can be assigned without pushing and popping.
    vector<Abalone> data;
    if (pid == 0)
    {
        Abalone temp;
        while (fin >> temp.sex >> temp.length >> temp.diameter >>
               temp.height >> temp.whole_height >> temp.shucked_weight >> temp.viscera_weight >> temp.shell_weight >> temp.rings)
        {
            abalones.push_back(temp);
        }
        data = abalones;
        clusterAssignment.resize(data.size());
    }

    std::vector<int> tempClusterAssignmentStorage;

    int abalonesArraySize = abalones.size();

    int chunksize = (abalonesArraySize / nproc);

    int leftover = (abalonesArraySize % nproc);

    MPI_Bcast(&abalonesArraySize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    if (pid == 0)
    {
        offset = chunksize + leftover;

        for (dest = 1; dest < nproc; dest++)
        {
            MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            offset = offset + chunksize;
        }

        int chunksize = (abalonesArraySize / nproc);
        int leftover = (abalonesArraySize % nproc);
        tempClusterAssignmentStorage.resize(chunksize + leftover);
    }
    else
    {
        MPI_Recv(&offset, 1, MPI_INT, MASTER, tag1, MPI_COMM_WORLD, &status);
        int chunksize = (abalonesArraySize / nproc);
        int leftover = (abalonesArraySize % nproc);
        tempClusterAssignmentStorage.resize(chunksize);
    }
    data.resize(abalonesArraySize);

    MPI_Barrier(MPI_COMM_WORLD);
    Timer totalSimulationTimer;
    vector<Abalone> clusterCenter;
    int clusterCenterSize;
    if (pid == 0)
    {
        clusterCenter = intializeRandomPesudo(data, K);
        clusterCenterSize = clusterCenter.size();
    }

    MPI_Bcast(&clusterCenterSize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    clusterCenter.resize(clusterCenterSize);
    MPI_Bcast(data.data(), sizeof(Abalone) * abalonesArraySize, MPI_CHAR, MASTER, MPI_COMM_WORLD);

    for (int iter = 0; iter < maxIter; iter++)
    {
        MPI_Bcast(clusterCenter.data(), sizeof(Abalone) * clusterCenterSize, MPI_CHAR, MASTER, MPI_COMM_WORLD);
        // pick each cluster assignment to minimize distance
        if (pid == MASTER)
        {
            /* ANDY: do we need to initialize the array */

            offset = 0;
            // START of update. range: (offset, chunksize+leftover, taskid);
            std::vector<Abalone> subAbalones =
                {data.begin() + offset, data.begin() + offset + chunksize + leftover};
            // std::vector<Particle> tempParticleStorage;

            // K-means place
            for (int i = 0; i < subAbalones.size(); i++)
            {
                // which cluster should the abalone belong to?
                int clusterBelong = 0;
                double minDistance = 1000000;
                for (int j = 0; j < clusterCenter.size(); j++)
                {
                    double distance = calculateDistanceEuclidean(subAbalones[i], clusterCenter[j], false);
                    // Found smaller stuff: we should update the distance.
                    if (distance < minDistance)
                    {
                        clusterBelong = j;
                        minDistance = distance;
                    }
                }
                tempClusterAssignmentStorage[i] = clusterBelong;
            }

            // overwrite tempParticleStorage to perform update for the master itself
            std::copy(std::begin(tempClusterAssignmentStorage), std::end(tempClusterAssignmentStorage), std::begin(clusterAssignment) + offset);
            // END of update

            // Wait to receive results from each task

            for (int i = 1; i < nproc; i++)
            {
                source = i;

                offset = leftover + chunksize * i;
                MPI_Recv(&clusterAssignment[offset], chunksize * sizeof(int), MPI_CHAR, source, tag2, MPI_COMM_WORLD, &status);
            }
        }
        /* end of master section */

        /***** Non-master tasks only *****/

        if (pid > MASTER)
        {
            /* Receive my portion of array from the master task */
            source = MASTER;
            // particles.resize(particleArraySize);
            //  MPI_Recv(particles.data(), particleArraySize * sizeof(Particle), MPI_CHAR, source, tag4, MPI_COMM_WORLD, &status);
            chunksize = (abalonesArraySize / nproc);

            // START of update. range: (offset, chunksize+leftover, taskid);

            std::vector<Abalone> subAbalones =
                {data.begin() + offset, data.begin() + offset + chunksize + leftover};

            for (int i = 0; i < subAbalones.size(); i++)
            {
                // which cluster should the abalone belong to?
                int clusterBelong = 0;
                double minDistance = 1000000;
                for (int j = 0; j < clusterCenter.size(); j++)
                {
                    double distance = calculateDistanceEuclidean(subAbalones[i], clusterCenter[j], false);
                    // Found smaller stuff: we should update the distance.
                    if (distance < minDistance)
                    {
                        clusterBelong = j;
                        minDistance = distance;
                    }
                }
                tempClusterAssignmentStorage[i] = clusterBelong;
            }

            // END of update

            dest = MASTER;

            MPI_Send(tempClusterAssignmentStorage.data(), sizeof(int) * tempClusterAssignmentStorage.size(), MPI_CHAR, MASTER, tag2, MPI_COMM_WORLD);
        }
        /*
        for (int i = 0; i < data.size(); i++)
        {
            // which cluster should the abalone belong to?
            int clusterBelong = 0;
            double minDistance = 1000000;
            for (int j = 0; j < clusterCenter.size(); j++)
            {
                double distance = calculateDistanceEuclidean(data[i], clusterCenter[j], false);
                // Found smaller stuff: we should update the distance.
                if (distance < minDistance)
                {
                    clusterBelong = j;
                    minDistance = distance;
                }
            }
            clusterAssignment[i] = clusterBelong;*/

        // update each cluster to a new value.
        // let master gather and reassign this stuff
        if (pid == MASTER)
        {
            // work on the ith cluster, and check if jth assignment belongs to it.
            for (int i = 0; i < clusterCenter.size(); i++)
            {
                Abalone clusterCenterAbalone = Abalone('M', 0, 0, 0, 0, 0, 0, 0, 0);
                int clusterBelongingCount = 0;
                for (int j = 0; j < clusterAssignment.size(); j++)
                {
                    if (clusterAssignment[j] == i)
                    {
                        abaloneAddition(clusterCenterAbalone, data[j]);
                        clusterBelongingCount++;
                    }
                }
                abaloneAverage(clusterCenterAbalone, clusterBelongingCount);
                clusterCenter[i] = clusterCenterAbalone;
                // printf("%d\n", clusterBelongingCount);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double totalSimulationTime = totalSimulationTimer.elapsed();

    if (pid == 0)
    {
        for (int i = 0; i < clusterCenter.size(); i++)
        {
            std::cout << clusterCenter[i].show() << std::endl;
        }
        printf("total simulation time: %.6fs\n", totalSimulationTime);
        // saveToFile(options.outputFile, particles);
    }

    MPI_Finalize();
}