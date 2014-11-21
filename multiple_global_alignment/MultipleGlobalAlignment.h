/*
 * MATH 578A Homework3
 * Multiple Global Alignment Algorithm
 * MultipleGlobalAlignment.cpp
 *
 * Author: Haifeng Chen
 * Contact: haifengc at usc dot edu

 * Compiler:
 * (1) Ubuntu 12.10 32-bit, g++ (Ubuntu/Linaro 4.7.2-2ubuntu1) 4.7.2
 * (2) Windows 7 32-bit, MinGW GCC 4.7.2

 * Created on: Mar 9, 2013
 */

#ifndef MULTIPLEGLOBALALIGNMENT_H_
#define MULTIPLEGLOBALALIGNMENT_H_

#include <set>
#include <map>
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

namespace spaceMultipleGlobalAlignment
{

#define I (i - 1)
#define J (j - 1)
#define P (p - 1)
#define Q (q - 1)
#define DIAG ('a')
#define UP ('b')
#define LEFT ('c')
#define lf_inf (std::numeric_limits<double>::max())

struct multisetcmp
{
	bool operator()(const pair<double, int> & a,
			const pair<double, int> & b) const
	{
		return a.first > b.first;
	}
};

enum MatchLabel
{
	MATCH = 0, MISMATCH = 1, INDEL = 2, TWOindel = 3
};

enum Nucleotide
{
	A = 0, C = 1, G = 2, T = 3, d = 4
};

const char alphabet[5] = { 'A', 'C', 'G', 'T', '-' };

class CMultipleGlobalAlignment
{
public:
	CMultipleGlobalAlignment(const string & file, const vector<int> & weigth);
	~CMultipleGlobalAlignment();
	void runMultipleGlobalAlignment();
private:
	string filePath;
	vector<string> seqs;
	vector<int> w;

	vector<vector<double> > s;
	//vector<vector<double> > epq; //the profile value gotten from Pro1(i) and Pro2(j)
	vector<vector<char> > l;

	/* Dist stores the distance from i to j, which is n×n matrix.
	 *
	 * We should use a data structure, which can insert, delete and
	 * find element in log(n) time. Red-black tree insert, delete
	 * and find emelmetn in log(n) time.
	 *
	 * In STL, set, map, multiset, multimap are implemented by red-black tree.
	 * Here, we don't need to change the value of Distance, we just need to
	 * delete or insert, so we use multiset.
	 * reference http://nlp.stanford.edu/IR-book/html/htmledition/time-complexity-of-hac-1.html
	 */
	vector<multiset<pair<double, int>, multisetcmp> > setDist;
	vector<vector<double> > Dist;
	vector<int> Indicator; // tag[i] marks whether sequence i has been clustered.
	vector<vector<string> > cluster; // the sequence id in each cluster
	vector<vector<int> > clusterID;
	vector<vector<vector<double> > > clusterProfile;

	map<int, string> mapres;

	clock_t TimeStart, TimeEnd;
	double memory;

	void outPutCluster();
	void resultDisplay();
	void multipleGlobalAlign();
	void stringReverse(string & str);
	int readString(const string & strPath);
	void outputfastaFormat(const string & str);
	MatchLabel charMatch(const char & a, const char & b);
	void outputfastaFormat(ofstream & fout, const string & str);
	pair<double, char> max(const double & s1, const double & s2,
			const double & s3);
	void outputResultString(const string & str, const int & start,
			const int & end);
	void setMulitiStrAlignProfile(const vector<string> & str,
			vector<vector<double> > & lfPro);
	void setMulitiStrAlignProfile(const string & str,
			vector<vector<double> > & lfPro);
	double computeProfileIJ(const vector<vector<double> > & lfPro,
			const int & si);
	double computeProfileIJ(const vector<vector<double> > & lfPro1,
			const int & si, const vector<vector<double> > & lfPro2,
			const int & sj);
	double globalAlignScore(const vector<string> & str1,
			const vector<string> & str2, const vector<vector<double> > & lfPro1,
			const vector<vector<double> > & lfPro2);
	double globalAlignScore(const string & U, const string & V);
	void globalAlignPath(vector<string> & str1, vector<string> & str2,
			const vector<vector<double> > & lfPro1,
			const vector<vector<double> > & lfPro2);

};

}

#endif /* MULTIPLEGLOBALALIGNMENT_H_ */
