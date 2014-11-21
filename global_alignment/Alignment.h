#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "GlobalAlignment.h"

#include <set>
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>

using namespace std;
using namespace spaceGlobalAlignment;

const string ErrorInfo =
    "Please input correct parameters! For example:\n\
		>galign.exe -m2 -s1 -i2 seq1.txt seq2.txt\n\
		scoring function with +2 for match, -1 for mismatch and -2 for indels\
		to align two DNA sequences in seq1.txt and seq2.txt.\n";

int checkFilePath(const string & strPath) {
  ifstream fin(strPath.c_str());
  int r = 0;
  if (!fin.good())
    r = 0;
  else
    r = 1;
  fin.close();
  return r;
}
#endif /* ALIGNMENT_H_ */
