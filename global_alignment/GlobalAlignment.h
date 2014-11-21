/*
 * GlobalAlignment.cpp
 * MATH 578A Homework1
 * Global Alignment Algorithm
 *
 * Author: Haifeng Chen
 * Contact: haifengc at usc dot edu
 * Compiler: MinGW GCC 4.7.2
 * Created on: Feb 21, 2013
 */

#ifndef GLOBALALIGNMENT_H_
#define GLOBALALIGNMENT_H_

#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <Windows.h>

using namespace std;

namespace spaceGlobalAlignment {

#define I (i - 1)
#define J (j - 1)
#define P (p - 1)
#define Q (q - 1)
#define DIAG ('a')
#define UP ('b')
#define LEFT ('c')
#define mismatchColor FOREGROUND_RED
#define IndelColor FOREGROUND_GREEN
#define matchColor FOREGROUND_BLUE
#define background FOREGROUND_INTENSITY
#define SetConCol SetConsoleTextAttribute
#define setBGcolor SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), background)

enum MatchLabel {
  MATCH = 0,
  MISMATCH = 1,
  INDEL = 2
};

class CGlobalAlignment {
 public:
  CGlobalAlignment(const string & filePath1, const string & filePath2,
                   const vector<int> & weigth);
  ~CGlobalAlignment();
  void runGlobalAlignment();
 private:
  string filePathU;
  string filePathV;
  string U;
  string V;
  vector<int> w;
  string rU;
  string rV;
  int alignScore;
  double lfIdentity;

  vector<vector<int> > s;
  vector<vector<char> > l;

  size_t n;
  size_t m;

  HANDLE Handle;
  clock_t TimeStart, TimeEnd;
  double memory;

  ostringstream ossGlLog;

  void setWindowWL();
  void resultDisplay();
  void globalAlignAlgorithm();
  void stringReverse(string & str);
  void outputLog(const string & strOut);
  void outputfastaFormat(const string & str);
  int readString(const string & strPath, string & str);
  MatchLabel charMatch(const char & a, const char & b);
  void outputResultString(const string & str, const int & start,
                          const int & end, const vector<int> & color);pair<int, char> max(const int & s1, const int & s2, const int & s3);
};

}

#endif /* GLOBALALIGNMENT_H_ */
