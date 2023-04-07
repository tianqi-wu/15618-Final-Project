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
/**
 * Sequential version of KNN. Uses PriorityQueue to help.
 * An abstraction of what a single thread will do.
 * Will extend to OpenMP friendly version in the future.
 */
double KNN_sequential(vector<Abalone> data, int K, Abalone someAbalone) {
    priority_queue<abaloneKeyValue, vector<abaloneKeyValue>, greater<abaloneKeyValue>> pq;
    // populate the priority queue
    
    for(int i = 0; i < data.size(); i++) {
        double differenceValue = calculateDistanceEuclidean(data[i], someAbalone);
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

    // pass it in the KNN function
    int K = 20;
    
    // tentative testing data:
    //M 0.71 0.555 0.195 1.9485 0.9455 0.3765 0.495 12
    
    Abalone randAbalone = Abalone('M', 0.71, 0.555,0.195,1.9485,0.9455,0.3765,0.495, 12);
    printf("%f\n", KNN_sequential(abalones, K, randAbalone));
}
