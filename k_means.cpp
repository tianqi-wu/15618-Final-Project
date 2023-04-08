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

typedef pair<double, int> abaloneKeyValue;

/* Pesudo-random initialization.
 * For testing purposes. We want this algorithm to be deterministic.
 */

vector<Abalone> intializeRandomPesudo(vector<Abalone> data, int K) {
    vector<Abalone> newVector;
    int slice = data.size() / K;
    for(int i = 0; i < K; i++) {
        int rand_num = slice * i;
        // copy construct. They are no longer what they were any more.
        Abalone newAbalone = Abalone(data[rand_num]);
        newVector.push_back(newAbalone);
        std::cout << newAbalone.show() << std::endl; 
    }

    return newVector;
}

/**
 * Random initialization. So-called one of the worst in K-means.
 * Implemented this only to make the algorithm work. Will extend when I have time.
 * We don't need to get the id. It's not necessary at all.
*/
vector<Abalone> intializeRandom(vector<Abalone> data, int K) {
    vector<Abalone> newVector;
    
    for(int i = 0; i < K; i++) {
        int rand_num = rand() % data.size();
        // copy construct. They are no longer what they were any more.
        newVector.push_back(Abalone(data[rand_num]));
    }

    return newVector;
}

/**
 * Sequential version of K-means. Uses PriorityQueue to help.
 * An abstraction of what a single thread will do.
 * Will extend to OpenMP friendly version in the future.
 */

vector<Abalone> K_means_sequential(vector<Abalone> data, int K, int maxIter) {
    priority_queue<abaloneKeyValue, vector<abaloneKeyValue>, greater<abaloneKeyValue>> pq;
    // populate the priority queue
    vector<Abalone> clusterCenter = intializeRandomPesudo(data, K);
    vector<int> clusterAssignment; // each give its own id
    // resize so it can be assigned without pushing and popping.
    clusterAssignment.resize(data.size());

    for(int iter = 0; iter < maxIter; iter++) {
        // pick each cluster assignment to minimize distance
        for(int i = 0; i < data.size(); i++) {
            // which cluster should the abalone belong to?
            int clusterBelong = 0;
            double minDistance = 1000000;
            for(int j = 0; j < clusterCenter.size(); j++) {
                double distance =  calculateDistanceEuclidean(data[i], clusterCenter[j], false);
                // Found smaller stuff: we should update the distance.
                if(distance < minDistance) {
                    clusterBelong = j;
                    minDistance = distance;
                }
            }
            clusterAssignment[i] = clusterBelong;
        }

        // update each cluster to a new value.
        
        // work on the ith cluster, and check if jth assignment belongs to it.
        for(int i = 0; i < clusterCenter.size(); i++) {
            Abalone clusterCenterAbalone = Abalone('M',0,0,0,0,0,0,0,0);
            int clusterBelongingCount = 0;
            for(int j = 0; j < clusterAssignment.size(); j++) {
                if(clusterAssignment[j] == i) {
                    abaloneAddition(clusterCenterAbalone, data[j]);
                    clusterBelongingCount++;
                }
            }
            abaloneAverage(clusterCenterAbalone, clusterBelongingCount);
            clusterCenter[i] = clusterCenterAbalone;
            //printf("%d\n", clusterBelongingCount);
        }
    }
   return clusterCenter;
}


/**
 * Driver function.
 * Used to run the thing correctly
 */
int main(int argc, const char **argv)
{
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
    Abalone temp;
    while (fin >> temp.sex >> temp.length >> temp.diameter >> 
        temp.height >> temp.whole_height >> temp.shucked_weight
        >> temp.viscera_weight >> temp.shell_weight >> temp.rings)
    {
        abalones.push_back(temp);
    }

    // pass it in the K-means hyperparameter function
    int K = 5;
    int maxIter = 1;
        
    vector<Abalone> result = K_means_sequential(abalones, K, maxIter);
    for(int i = 0; i < result.size(); i++) {
        std::cout << result[i].show() << std::endl; 
    }
}
