/*
 * MATH 578A Homework2
 * Banded Global Alignment Algorithm
 * BandedGlobalAlignment.cpp
 *
 * Author: Haifeng Chen
 * Contact: haifengc at usc dot edu

 * Compiler:
 * (1) Ubuntu 12.10 32-bit, g++ (Ubuntu/Linaro 4.7.2-2ubuntu1) 4.7.2
 * (2) Windows 7 32-bit, MinGW GCC 4.7.2

 * Created on: Mar 3, 2013 - 2pm
 */

#ifndef BANDEDGLOBALALIGNMENT_H_
#define BANDEDGLOBALALIGNMENT_H_

#include <ctime>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>

using namespace std;

namespace spaceBandedGlobalAlignment
{

#define I (i - 1)
#define J (j - 1)
#define P (p - 1)
#define Q (q - 1)
#define DIAG ('a')
#define UP ('b')
#define LEFT ('c')
#define slCOL(row, col) ((col) - L[(row)])
#define inf (std::numeric_limits<int>::max())

enum MatchLabel
{
	MATCH = 0, MISMATCH = 1, INDEL = 2
};

class CBandedGlobalAlignment
{
public:
	CBandedGlobalAlignment(const string & filePath1, const string & filePath2,
			const vector<int> & weigth);
	~CBandedGlobalAlignment();
	void runBandedGlobalAlignment();
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

	vector<int> L;
	vector<int> R;

	int n;
	int m;
	int d;

	clock_t TimeStart, TimeEnd;
	double memory;

	ostringstream ossGlLog;

	void outPutsANDl();
	void resultDisplay();
	void setLR(const int & k);
	void stringReverse(string & str);
	void outputLog(const string & strOut);
	void outputfastaFormat(const string & str);
	bool bandedGlobalAlignAlgorithm(const int & k);
	int readString(const string & strPath, string & str);
	MatchLabel charMatch(const char & a, const char & b);
	void outputResultString(const string & str, const int & start,
			const int & end);
	pair<int, char> max(const int & s1, const int & s2, const int & s3);
};

}

#endif /* BANDEDGLOBALALIGNMENT_H_ */
