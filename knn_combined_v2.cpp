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


float KNN_parallel(vector<Abalone> data,int K, Abalone someAbalone) {
  //priority_queue<abaloneKeyValue, vector<abaloneKeyValue>, greater<abaloneKeyValue>> pq;
    // populate the priority queue

  int sliced_size = data.size() / 8;
  std::vector<std::pair<double,int>> vector_for_each[8];

  #pragma omp parallel for
  for (int i = 0; i < 8; i++){
    int start_idx = sliced_size * i;
    int end_idx = i== 7 ? data.size(): sliced_size * (i+1);

    #pragma omp simd
    for (int j = start_idx; j < end_idx; j++){
      double differenceValue = calculateDistanceEuclidean(data[j], someAbalone, true);
      int ringNumber = data[j].rings;
      vector_for_each[i].push_back(make_pair(differenceValue, ringNumber));
    }
    
  }
  
  std::vector<std::pair<double,int>> joined_vector;
  for (int i = 0; i < 8; i++){
  
    joined_vector.insert(joined_vector.end(), vector_for_each[i].begin(), vector_for_each[i].end());
  }


  
  //std::sort(joined_vector.begin(),joined_vector.end());
 
  std::sort(joined_vector.begin(),joined_vector.end());
  //joined_vector.reverse();
  //std::reverse(joined_vector.begin(),joined_vector.end());

  
  double sum = 0;
  for (int i = 0; i < K; i++){
    sum += joined_vector[i].second;
  }
  //assert(joined_vector.size() > K);
  return sum / K;
}



int main(int argc, char *argv[])
{
  int pid;
  int nproc;
  int offset;
  int source;
  int dest;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int tag4 = 1;
  int tag3 = 2;
  int tag2 = 3;
  int tag1 = 4;
  
  int K = 20;
  
  string location = "./data/abalone.data";
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
  vector<float> abaloneResult; // newParticles
  vector<float> tempAbaloneResultStorage;
  
  vector<Abalone> training;
  vector<Abalone> testing;

  //initialization
  if (pid == 0){
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

  int trainingSize = training.size();
  int testingSize = testing.size();

  int chunksize = (testingSize / nproc);
  int leftover = (testingSize % nproc);

  //tempAbaloneResultStorage.resize(chunksize + leftover);

  Timer totalSimulationTimer;

  MPI_Bcast(&trainingSize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&testingSize, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

  int offset_arr[nproc];
  
  if (pid == 0)
    {
      offset = chunksize + leftover;
      
      for (dest = 1; dest < nproc; dest++)
        {
	  offset_arr[dest] = offset;
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
  //MPI_Bcast(testing.data(), sizeof(Abalone) * testingSize, MPI_CHAR, MASTER, MPI_COMM_WORLD);

  if (pid == MASTER){
    for(int i = 1; i < nproc; i++){
      MPI_Send(testing[offset_arr[i]]


    
    offset = 0;


    
    std::vector<Abalone> subAbalones =
      {testing.begin() + offset, testing.begin() + offset + chunksize + leftover};

    for (int i = 0; i < subAbalones.size(); i++){
      tempAbaloneResultStorage[i] = KNN_parallel(training, K, subAbalones[i]);
    }
    
    // overwrite tempParticleStorage to perform update for the master itself
    std::copy(std::begin(tempAbaloneResultStorage), std::end(tempAbaloneResultStorage), std::begin(abaloneResult) + offset);
    
        // Wait to receive results from each task
    
    for (int i = 1; i < nproc; i++)
      {
	source = i;
	offset = leftover + chunksize * i;
	MPI_Recv(&abaloneResult[offset], chunksize * sizeof(float), MPI_CHAR, source, tag2, MPI_COMM_WORLD, &status);
      }
  }

  if (pid > MASTER)
    {
        /* Receive my portion of array from the master task */
        source = MASTER;
        // test.resize(testingSize);
        // MPI_Recv(particles.data(), particleArraySize * sizeof(Particle), MPI_CHAR, source, tag4, MPI_COMM_WORLD, &status);
        chunksize = (testingSize / nproc);

        // START of update. range: (offset, chunksize+leftover, taskid);

        std::vector<Abalone> subAbalones =
            {testing.begin() + offset, testing.begin() + offset + chunksize};

        // KNN!!!!!
        for (int i = 0; i < subAbalones.size(); i++)
        {
            tempAbaloneResultStorage[i] = KNN_parallel(training, K, subAbalones[i]);
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
