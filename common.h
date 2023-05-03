#ifndef COMMON_H_
#define COMMON_H_

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

inline void saveToFile(std::string fileName,
                       const std::vector<float> &abaloneResult) {
  std::ofstream f(fileName);
  assert((bool)f && "Cannot open output file");

  f << std::setprecision(9);
  for (float p : abaloneResult) {
    f << p << std::endl;
  }
  assert((bool)f && "Failed to write to output file");
}

inline void saveToFileKMeans(std::string fileName,
                       const std::vector<int> &kMeansResult) {
  std::ofstream f(fileName);
  assert((bool)f && "Cannot open output file");

  f << std::setprecision(9);
  for (int p : kMeansResult) {
    f << p << std::endl;
  }
  assert((bool)f && "Failed to write to output file");
}

#endif