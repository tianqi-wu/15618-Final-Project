/*
	Name		Data Type	Meas.	Description
	----		---------	-----	-----------
	Sex		nominal			M, F, and I (infant)
	Length		continuous	mm	Longest shell measurement
	Diameter	continuous	mm	perpendicular to length
	Height		continuous	mm	with meat in shell
	Whole weight	continuous	grams	whole abalone
	Shucked weight	continuous	grams	weight of meat
	Viscera weight	continuous	grams	gut weight (after bleeding)
	Shell weight	continuous	grams	after being dried
	Rings		integer			+1.5 gives the age in years
*/
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

class Abalone
{
public:
	char sex;
	double length;
	double diameter;
	double height;
	double whole_height;
	double shucked_weight;
	double viscera_weight;
	double shell_weight;
	double rings;
	Abalone(){} // Constructor without parameters
	
	// copy constructor
	// https://www.mygreatlearning.com/blog/copy-constructor-in-cpp/
	Abalone(const Abalone &obj) {
   		this->sex = obj.sex;
		this->sex = obj.sex;
		this->length = obj.length;
		this->diameter = obj.diameter;
		this->height = obj.height;
		this->whole_height = obj.whole_height;
		this->shucked_weight = obj.shucked_weight;
		this->viscera_weight = obj.viscera_weight;
		this->shell_weight = obj.shell_weight;
		this->rings = obj.rings;
	}

	// default constructor
	Abalone(char sex,
			double length,
			double diameter,
			double height,
			double whole_height,
			double shucked_weight,
			double viscera_weight,
			double shell_weight,
			double rings) { // Constructor with parameters
		this->sex = sex;
		this->length = length;
		this->diameter = diameter;
		this->height = height;
		this->whole_height = whole_height;
		this->shucked_weight = shucked_weight;
		this->viscera_weight = viscera_weight;
		this->shell_weight = shell_weight;
		this->rings = rings;
	}


// https://stackoverflow.com/questions/30904499/printing-out-a-class-objects-variables-from-a-stack

	string show() {
   	stringstream ss;
   		ss << this->sex 
		<< "," << this->length
		<< "," << this->diameter
		<< "," << this->height 
		<< "," << this->whole_height
		<< "," << this->shucked_weight
		<< "," << this->viscera_weight
		<< "," << this->shell_weight
		<< "," << this->rings
		<< std::endl;
   		return (ss.str());
	}
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
double calculateDistanceEuclidean(Abalone abalone1, Abalone abalone2, bool isKNN) {
	int sexDiff = 0;//abalone1.sex != abalone2.sex ? 1 : 0;
    double lengthDiff = euclidean(abalone1.length, abalone2.length);
	double diameterDiff = euclidean(abalone1.diameter, abalone2.diameter);
	double heightDiff = euclidean(abalone1.height, abalone2.height);
	double wholeHeightDiff = euclidean(abalone1.whole_height, abalone2.whole_height);
	double shuckedWeightDiff = euclidean(abalone1.shucked_weight, abalone2.shucked_weight);
	double visceraWeightDiff = euclidean(abalone1.viscera_weight, abalone2.viscera_weight);
	double shellWeightDiff = euclidean(abalone1.shell_weight, abalone2.shell_weight);
	double ringDiff = 0;
	if(!isKNN) {
		ringDiff = euclidean(abalone1.rings, abalone2.rings);
	}
	return sexDiff + lengthDiff + diameterDiff + heightDiff + 
		wholeHeightDiff + shuckedWeightDiff + visceraWeightDiff + shellWeightDiff + ringDiff;
}

double calculateDistanceEuclideanSimd(Abalone abalone1, Abalone abalone2, bool isKNN){
  double sum = 0.0;
  
#pragma omp simd reduction(+:sum)
  for(int i = 0; i < 8;i++){
    switch(i){
      case 0:
	sum += (abalone1.length - abalone2.length)* (abalone1.length - abalone2.length);
	break;
      case 1:
	sum += (abalone1.diameter - abalone2.diameter) * (abalone1.diameter - abalone2.diameter);
	break;
      case 2:
        sum += (abalone1.height - abalone2.height) * (abalone1.height - abalone2.height);
	break;
      case 3:
	sum += (abalone1.whole_height - abalone2.whole_height) * (abalone1.whole_height - abalone2.whole_height);
	break;
      case 4:
	sum += (abalone1.shucked_weight - abalone2.shucked_weight) * (abalone1.shucked_weight - abalone2.shucked_weight);
	break;
      case 5:
	sum += (abalone1.viscera_weight - abalone2.viscera_weight) * (abalone1.viscera_weight - abalone2.viscera_weight);
	break;
      case 6:
	sum += (abalone1.shell_weight - abalone2.shell_weight) * (abalone1.shell_weight - abalone2.shell_weight);
	break;
      case 7:
	if(!isKNN){
	  sum += (abalone1.rings - abalone2.rings) * (abalone1.rings - abalone2.rings);
	}
	break;
    
    
    }
    
    
  }
  return sum;
}

/**
 * Creating a monster-size abalone by adding up their dimensions.
 */
void abaloneAddition(Abalone &abalone1, Abalone &abalone2) {
	abalone1.length += abalone2.length;
	abalone1.diameter += abalone2.diameter;
	abalone1.height += abalone2.height;
	abalone1.whole_height += abalone2.whole_height;
	abalone1.shucked_weight += abalone2.shucked_weight;
	abalone1.viscera_weight += abalone2.viscera_weight;
	abalone1.shell_weight += abalone2.shell_weight;
	abalone1.rings += abalone2.rings;
}

/**
 * Creating a normal-size abalone from a monster-size by adding up their dimensions.
 */
void abaloneAverage(Abalone &abalone1, int clusterBelongingCount) {
	if(clusterBelongingCount == 0) {
		return;
	}
	abalone1.length /= clusterBelongingCount;
	abalone1.diameter /= clusterBelongingCount;
	abalone1.height /= clusterBelongingCount;
	abalone1.whole_height /= clusterBelongingCount;
	abalone1.shucked_weight /= clusterBelongingCount;
	abalone1.viscera_weight /= clusterBelongingCount;
	abalone1.shell_weight /= clusterBelongingCount;
	abalone1.rings /= clusterBelongingCount;
}
