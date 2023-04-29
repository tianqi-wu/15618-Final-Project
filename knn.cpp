/**
 * KNN based on UCI abalone data.
 * we have specified the input correctluy so
 */


using namespace std;
#include "abalone.h"
#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <queue>
#include <omp.h>
#include <cmath>
#include "timing.h"
#include <assert.h>
#include <algorithm>
typedef pair<double, int> abaloneKeyValue;
/**
 * Sequential version of KNN. Uses PriorityQueue to help.
 * An abstraction of what a single thread will do.
 * Will extend to OpenMP friendly version in the future.
 */
double KNN_sequential(vector<Abalone> data, int K, Abalone someAbalone) {
    priority_queue<abaloneKeyValue, vector<abaloneKeyValue>, greater<abaloneKeyValue>> pq;
    // populate the priority queue
    
    for(int i = 0; i < data.size(); i++) {
        double differenceValue = calculateDistanceEuclidean(data[i], someAbalone, true);
        int ringNumber = data[i].rings;
        pq.push(make_pair(differenceValue, ringNumber));
        //printf("%f %d\n", differenceValue, ringNumber);
    }
    double sum = 0;
    for(int i = 0; i < K; i++) {
        pair<double, int> top = pq.top();
        sum += top.second;
        pq.pop();
    }
    return sum / K;
}


double KNN_sequential_simd(vector<Abalone> data, int K, Abalone someAbalone) {
    priority_queue<abaloneKeyValue, vector<abaloneKeyValue>, greater<abaloneKeyValue>> pq;
    // populate the priority queue

    #pragma omp simd
    for(int i = 0; i < data.size(); i++) {
        double differenceValue = calculateDistanceEuclideanSimd(data[i], someAbalone, true);
        int ringNumber = data[i].rings;
        pq.push(make_pair(differenceValue, ringNumber));
    }
    double sum = 0;
    for(int i = 0; i < K; i++) {
        pair<double, int> top = pq.top();
        sum += top.second;
        pq.pop();
    }
    return sum / K;
}

double KNN_parallel(vector<Abalone> data, int K, Abalone someAbalone) {
  //priority_queue<abaloneKeyValue, vector<abaloneKeyValue>, greater<abaloneKeyValue>> pq;
    // populate the priority queue

  int sliced_size = data.size() / 128;
  std::vector<std::pair<double,int>> vector_for_each[128];

  #pragma omp parallel for
  for (int i = 0; i < 128; i++){
    int start_idx = sliced_size * i;
    int end_idx = i== 127 ? data.size(): sliced_size * (i+1);

    for (int j = start_idx; j < end_idx; j++){
     
      
      double differenceValue = calculateDistanceEuclidean(data[j], someAbalone, true);
      int ringNumber = data[j].rings;
      vector_for_each[i].push_back(make_pair(differenceValue, ringNumber));
    }
    
  }
  
  std::vector<std::pair<double,int>> joined_vector;
  for (int i = 0; i < 128; i++){
    joined_vector.insert(joined_vector.end(), vector_for_each[i].begin(), vector_for_each[i].end());
  }
  
  std::sort(joined_vector.begin(),joined_vector.end());
  
  double sum = 0;
  for (int i = 0; i < K; i++){
    sum += joined_vector[i].second;
  }
  
  
  /*double sum = 0;
    
    for(int i = 0; i < K; i++) {
        pair<double, int> top = pq.top();
        sum += top.second;
        pq.pop();
	}*/
  return sum / K;
}


/**
 * Separates data in a 70%/30% fashion.
 * This is guaranteed to be deterministic. 
 */

void pesudo_training_test_parse(vector<Abalone> &training, vector<Abalone> &testing)  {
    vector<Abalone> updatedTrainingData;
    for(int i = 0; i < training.size(); i++) {
        if(i % 8 == 0 || i % 6 == 0) {
            testing.push_back(training[i]);
        }else {
            updatedTrainingData.push_back(training[i]);
        }
    }
    training = updatedTrainingData;
}


/**
 * Driver function.
 * Used to run the thing correctly
 */
int main(int argc, const char **argv)
{
    string location = "./data/mass_abalone.data";
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
    Abalone temp;
    while (fin >> temp.sex >> temp.length >> temp.diameter >> 
        temp.height >> temp.whole_height >> temp.shucked_weight
        >> temp.viscera_weight >> temp.shell_weight >> temp.rings)
    {
        abalones.push_back(temp);
    }

    // pass it in the KNN function
    int K = 20;
    
    // tentative testing data:
    //M 0.71 0.555 0.195 1.9485 0.9455 0.3765 0.495 12


    // generate testing data from the training data by random selection.
    // we follow a non-strict 7/3 rule: we extract data from it.

    
    // this is the part of code that has to be timed

    vector<Abalone> training = abalones;
    vector<Abalone> testing;

    pesudo_training_test_parse(training, testing);
    printf("Training Size: %lu\n", training.size());
    printf("Testing Size: %lu\n", testing.size());
    double myTime = 0;

    vector<double> sequential_result;
    sequential_result.resize(testing.size());

    vector<double> parallel_result;
    parallel_result.resize(testing.size());

    vector<double> over_training_result;
    over_training_result.resize(testing.size());
    
    /*Timer seqTimer;
    
    for(int i =0; i < testing.size(); i++) {
        Abalone currAbalone = testing[i];
        sequential_result[i] = KNN_sequential(training, K, currAbalone);
    }

    double seqTime = seqTimer.elapsed();
    std::cout << "seq runtime" << seqTime << std::endl;
    */

    
    Timer parallelTimer;

    #pragma omp parallel for schedule(static)
    for(int i =0; i < testing.size(); i++) {
        Abalone currAbalone = testing[i];
        //parallel_result[i] = KNN_parallel(training, K, currAbalone);
	parallel_result[i] = KNN_sequential(training, K, currAbalone);
    }

    double parallelTime = parallelTimer.elapsed();
    //std::cout<< "parallel output" << parallel_output << std::endl;
    std::cout<< "parallel runtime " << parallelTime << std::endl;


    Timer overTrainingTimer;
    for(int i = 0; i < testing.size(); i++){
      Abalone currAbalone = testing[i];
      over_training_result[i] = KNN_parallel(training, K, currAbalone);
    }
    double overTrainingTime = overTrainingTimer.elapsed();
    std::cout << "parallel over training runtime" << overTrainingTime << std::endl;

    

    /*bool testsPassed = true;
    for(int i = 0; i < testing.size(); i++) {
        if(sequential_result[i] != parallel_result[i]) {
            printf("Error at sample %d, should be %f, got %f\n", i, sequential_result[i], parallel_result[i]);
            testsPassed = false;
            break;
        }
    }
    if(testsPassed) {
        printf("All tests passed! Congratulations!\n");
	}*/
    
}
