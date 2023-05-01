#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <queue>
#include <sstream>
#include <omp.h>

class Abalone_Simd
{
public:
  char sex;
  double attr_arr[8];
     
  Abalone_Simd(){} // Constructor without parameters
	
	// copy constructor
	// https://www.mygreatlearning.com/blog/copy-constructor-in-cpp/
	Abalone_Simd(const Abalone_Simd &obj) {
	  this->sex = obj.sex;
	  for (int i = 0; i < 8; i++){
	    this->attr_arr[i] = obj.attr_arr[i];
	  }
	}

	// default constructor
        Abalone_Simd(char sex,
		double length,
		double diameter,
		double height,
		double whole_height,
		double shucked_weight,
		double viscera_weight,
		double shell_weight,
		double rings) { // Constructor with parameters
	  this->sex = sex;
	  this->attr_arr[0] = length;
	  this->attr_arr[1] = diameter;
	  this->attr_arr[2] = height;
	  this->attr_arr[3] = whole_height;
	  this->attr_arr[4] = shucked_weight;
	  this->attr_arr[5] = viscera_weight;
	  this->attr_arr[6] = shell_weight;
	  this->attr_arr[7] = rings;
	}


// https://stackoverflow.com/questions/30904499/printing-out-a-class-objects-variables-from-a-stack

  /*string show() {
   	stringstream ss;
   	        
		this->attr_arr[0]
		<< "," << this->attr_arr[1]
		<< "," << this->attr_arr[2]
		<< "," << this->attr_arr[3]
		<< "," << this->attr_arr[4]
		<< "," << this->attr_arr[5]
		<< "," << this->attr_arr[6]
		<< "," << this->attr_arr[7]
		<< std::endl;
   		return (ss.str());
		}*/
};

// https://www.geeksforgeeks.org/priority-queue-of-pairs-in-c-ordered-by-first/
/**
 * Just a simple Euclidean calculation function.
 */
double euclidean(double p1, double p2) {
	return (p1 - p2) * (p1 - p2);
}

/* Should have 8 in total. Ring is the y data that we should be learning from.
 * For categorical, we used jaccard. For numerical, euclidian.
 * given the goals of KNN/K-means is different, I added this bool value.
 */
double calculateDistanceEuclidean(Abalone_Simd abalone1, Abalone_Simd abalone2, bool isKNN){
  double sum = 0.0;
  double *arr1 = abalone1.attr_arr;
  double *arr2 = abalone2.attr_arr;
  
  #pragma omp simd reduction(+:sum)
  for(int i = 0; i < 8;i++){
    sum += arr1[i] + arr2[i];
  }
    
    
  return sum;
}

/**
 * Creating a monster-size abalone by adding up their dimensions.
 */
void abaloneAddition(Abalone_Simd &abalone1, Abalone_Simd &abalone2) {
  #pragma omp simd
  for(int i = 0; i < 8; i++){
    abalone1.attr_arr[i] += abalone2.attr_arr[i];
  }
}

/**
 * Creating a normal-size abalone from a monster-size by adding up their dimensions.
 */
void abaloneAverage(Abalone_Simd &abalone1, int clusterBelongingCount) {
  if(clusterBelongingCount == 0) {
    return;
  }
  #pragma omp simd
  for(int i = 0; i < 8; i++){
    abalone1.attr_arr[i] /= clusterBelongingCount;
  }
  
}
