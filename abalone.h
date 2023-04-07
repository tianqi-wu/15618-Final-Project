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
	int rings;
	Abalone(){} // Constructor without parameters
	
	Abalone(char sex,
			double length,
			double diameter,
			double height,
			double whole_height,
			double shucked_weight,
			double viscera_weight,
			double shell_weight,
			int rings)
	{ // Constructor with parameters
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
 */
double calculateDistanceEuclidean(Abalone abalone1, Abalone abalone2) {
	int sexDiff = abalone1.sex != abalone2.sex ? 1 : 0;
    double lengthDiff = euclidean(abalone1.length, abalone2.length);
	double diameterDiff = euclidean(abalone1.diameter, abalone2.diameter);
	double heightDiff = euclidean(abalone1.height, abalone2.height);
	double wholeHeightDiff = euclidean(abalone1.whole_height, abalone2.whole_height);
	double shuckedWeightDiff = euclidean(abalone1.shucked_weight, abalone2.shucked_weight);
	double visceraWeightDiff = euclidean(abalone1.viscera_weight, abalone2.viscera_weight);
	double shellWeightDiff = euclidean(abalone1.shell_weight, abalone2.shell_weight);
	return sexDiff + lengthDiff + diameterDiff + heightDiff + 
		wholeHeightDiff + shuckedWeightDiff + visceraWeightDiff + shellWeightDiff;
}