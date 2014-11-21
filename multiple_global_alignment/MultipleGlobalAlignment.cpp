#include "MultipleGlobalAlignment.h"

using namespace spaceMultipleGlobalAlignment;

CMultipleGlobalAlignment::CMultipleGlobalAlignment(const string & file,
		const vector<int> & weight)
{
	filePath = file;
	w = weight;
	memory = 0.0;

	TimeStart = 0.0;
	TimeEnd = 0.0;
}

CMultipleGlobalAlignment::~CMultipleGlobalAlignment()
{

}

int CMultipleGlobalAlignment::readString(const string & strPath)
{
	ifstream fin(strPath.c_str());
	if (!fin.good()) {
		printf("Cannot open the file: %s.\n", strPath.c_str());
		fin.close();
		return 0;
	}
	bool vaild = true;
	string strTmp, strIn;
	strIn.clear();
	seqs.clear();
	while (!fin.eof()) {
		getline(fin, strTmp, '\n');
		if (strTmp.size() == 0)
			continue;
		for (int i = 0; i < (int) strTmp.size(); i++) {
			if (strTmp[0] == '>') {
				if (strIn.size() != 0) {
					seqs.push_back(strIn);
					strIn.clear();
				}
				break;
			}

			if (strTmp[i] != 'A' && strTmp[i] != 'C' && strTmp[i] != 'G'
					&& strTmp[i] != 'T' && strTmp[i] != 'N' && strTmp[i] != 'a'
					&& strTmp[i] != 'c' && strTmp[i] != 'g' && strTmp[i] != 't'
					&& strTmp[i] != 'n') {
				vaild = false;
			} else {
				strIn += strTmp[i];
			}
		}
	}
	if (strIn.size() != 0) {
		seqs.push_back(strIn);
		strIn.clear();
	}
	if (vaild == false) {
		printf(
				"The DNA file %s contains characters which is not a, c, g, t or n. \
				These characters have been deleted.\n",
				strPath.c_str());
	}
	fin.close();
	return 1;
}

inline void CMultipleGlobalAlignment::stringReverse(string & str)
{
	int n = str.size();
	char c;
	for (int i = 0; i < n / 2; i++) {
		c = str[i];
		str[i] = str[n - i - 1];
		str[n - i - 1] = c;
	}
}

inline MatchLabel CMultipleGlobalAlignment::charMatch(const char & a,
		const char & b)
{
	if (a == '-' && b == '-')
		return TWOindel;
	else if (a == '-' || b == '-')
		return INDEL;
	else if (a == b)
		return MATCH;
	else
		return MISMATCH;
}

inline pair<double, char> CMultipleGlobalAlignment::max(const double & s1,
		const double & s2, const double & s3)
{
	/*if two of them are equal, then there are more than one optimal path*/
	if (s1 >= s2) {
		if (s1 >= s3)
			return pair<double, char>(s1, DIAG);
		else
			return pair<double, char>(s3, LEFT);
	} else {
		if (s2 >= s3)
			return pair<double, char>(s2, UP);
		else
			return pair<double, char>(s3, LEFT);
	}
}

void CMultipleGlobalAlignment::outputfastaFormat(const string & str)
{
	for (size_t t = 0; t < str.size(); t++) {
		if (t % 80 == 0 && t != 0) {
			printf("\n");
		}
		printf("%c", str[t]);
	}
	printf("\n");
}

void CMultipleGlobalAlignment::outputfastaFormat(ofstream & fout,
		const string & str)
{
	for (size_t t = 0; t < str.size(); t++) {
		if (t % 80 == 0 && t != 0) {
			fout << endl;
		}
		fout << str[t];
	}
	fout << endl;
}

void CMultipleGlobalAlignment::resultDisplay()
{
	string strfile = filePath;
	size_t pos = strfile.find(".txt");
	strfile = strfile.substr(0, pos);

	char chr[100];
	sprintf(chr, "%s_res.txt", strfile.c_str());

	ofstream fout(chr);
	size_t length = 0;
	for (size_t i = 0; i < cluster.size(); i++) {
		if (Indicator[i] && cluster[i].size() != 0) {
			if (cluster[i].size() > 0)
				length = cluster[i][0].size();
			for (size_t j = 0; j < cluster[i].size(); j++) {
				mapres[clusterID[i][j]] = cluster[i][j];
			}
		}
	}

	fout << "Progressive multiple global alignment result: "
			<< "(alignment string length = " << length
			<< ", number of string = " << mapres.size() << ")" << endl;

	for (map<int, string>::iterator it = mapres.begin(); it != mapres.end();
			it++) {
		outputfastaFormat(fout, it->second);
	}
	fout << endl;
	fout << "Running time: " << (double) (TimeEnd - TimeStart) / CLOCKS_PER_SEC
			<< "s" << endl;
	fout.close();
}

void CMultipleGlobalAlignment::setMulitiStrAlignProfile(
		const vector<string> & str, vector<vector<double> > & lfPro)
{
	/* Time  and Space complexity
	 * Time O(5n)
	 * Space O(5n)
	 * */
	if (str.size() == 0)
		return;
	lfPro.resize(5);
	size_t n = str[0].size();
	vector<double> lfv(n);
	for (size_t t = 0; t < lfPro.size(); t++) {
		lfPro[t] = lfv;
	}
	lfv.clear();
	size_t r = str.size();
	for (size_t i = 0; i < n; i++) {
		int cnt[5] = { 0 };
		for (size_t j = 0; j < r; j++) {
			if (str[j][i] == 'a' || str[j][i] == 'A')
				cnt[A]++;
			else if (str[j][i] == 'c' || str[j][i] == 'C')
				cnt[C]++;
			else if (str[j][i] == 'g' || str[j][i] == 'G')
				cnt[G]++;
			else if (str[j][i] == 't' || str[j][i] == 'T')
				cnt[T]++;
			else if (str[j][i] == '-')
				cnt[d]++;
		}
		for (int j = 0; j < 5; j++) {
			if (cnt[j] == 0)
				lfPro[j][i] = 0.0;
			else
				lfPro[j][i] = (double) cnt[j] / (double) str.size();
		}
	}
}

void CMultipleGlobalAlignment::setMulitiStrAlignProfile(const string & str,
		vector<vector<double> > & lfPro)
{
	vector<string> sstr;
	sstr.push_back(str);
	setMulitiStrAlignProfile(sstr, lfPro);
}

inline double CMultipleGlobalAlignment::computeProfileIJ(
		const vector<vector<double> > & lfPro1, const int & si,
		const vector<vector<double> > & lfPro2, const int & sj)
{
	/* Time  and Space complexity
	 * Time O(25)
	 * no extra Space
	 * */
	double sum = 0.0;

	double s = 0.0;
	for (size_t i = 0; i < 4; i++) {
		s += lfPro1[i][si] * lfPro2[i][sj];
	}
	sum += s * w[MATCH];

	s = 0.0;
	for (size_t i = 0; i < 4; i++) {
		s += lfPro1[i][si]
				* (lfPro2[0][sj] + lfPro2[1][sj] + lfPro2[2][sj] + lfPro2[3][sj]
						- lfPro2[i][sj]);
	}
	sum += s * w[MISMATCH];

	s = 0.0;
	s += lfPro1[4][si] * (1 - lfPro2[4][sj]);
	s += lfPro2[4][sj] * (1 - lfPro1[4][si]);
	sum += s * w[INDEL];

	return sum;
}

inline double CMultipleGlobalAlignment::computeProfileIJ(
		const vector<vector<double> > & lfPro, const int & si)
{
	double sum = 0.0;
	sum = lfPro[0][si] + lfPro[1][si] + lfPro[2][si] + lfPro[3][si];
	return sum * w[INDEL];
}

double CMultipleGlobalAlignment::globalAlignScore(const string & U,
		const string & V)
{
	/* Time  and Space complexity
	 * Time O(mn)
	 * Space O(mn)
	 * here we can use linear space pairwise alignment
	 * actually, we can use pairwiseGlobalAlignScore(const vector<string> & str1...)
	 * to calculate the score of two string, but this function will save time.
	 * */
	size_t n = U.size(), m = V.size();

	s.resize(n + 1);

	for (size_t i = 0; i <= n; i++) {
		s[i].resize(m + 1, 0);
	}

	s[0][0] = 0;
	for (size_t j = 1; j <= m; j++) {
		s[0][j] = s[0][j - 1] + w[INDEL];
	}

	for (size_t i = 1; i <= n; i++) {
		s[i][0] = s[i - 1][0] + w[INDEL];
		for (size_t j = 1; j <= m; j++) {
			double s1, s2, s3;
			s1 = s[i - 1][j - 1] + w[charMatch(U[I], V[J])];
			s2 = s[i - 1][j] + w[INDEL];
			s3 = s[i][j - 1] + w[INDEL];
			pair<double, char> charMatchResult = max(s1, s2, s3);
			s[i][j] = charMatchResult.first;
		}
	}

	return s[n][m];
}

double CMultipleGlobalAlignment::globalAlignScore(const vector<string> & str1,
		const vector<string> & str2, const vector<vector<double> > & lfPro1,
		const vector<vector<double> > & lfPro2)
{
	/* Time  and Space complexity
	 * Time O(25mn)
	 * Space O(mn)
	 * here we can use linear space pairwise alignment
	 * */
	if (str1.size() == 1 && str2.size() == 1) {
		string U = str1[0];
		string V = str2[0];
		return globalAlignScore(U, V);
	}

	size_t n = str1[0].size(), m = str2[0].size();

	s.resize(n + 1);

	for (size_t i = 0; i <= n; i++) {
		s[i].resize(m + 1, 0);
	}

	s[0][0] = 0;
	for (size_t j = 1; j <= m; j++) {
		s[0][j] = s[0][j - 1] + computeProfileIJ(lfPro2, J);
	}

	for (size_t i = 1; i <= n; i++) {
		s[i][0] = s[i - 1][0] + computeProfileIJ(lfPro1, I);
		for (size_t j = 1; j <= m; j++) {
			double s1, s2, s3;
			s1 = s[i - 1][j - 1] + computeProfileIJ(lfPro1, I, lfPro2, J);
			s2 = s[i - 1][j] + computeProfileIJ(lfPro1, I);
			s3 = s[i][j - 1] + computeProfileIJ(lfPro2, J);
			pair<double, char> charMatchResult = max(s1, s2, s3);
			s[i][j] = charMatchResult.first;
		}
	}

	return s[n][m];
}

void CMultipleGlobalAlignment::globalAlignPath(vector<string> & str1,
		vector<string> & str2, const vector<vector<double> > & lfPro1,
		const vector<vector<double> > & lfPro2)
{
	/*pairwiseGlobalAlignScore pairwiseGlobalAlignPath have much common code,
	 * but I use two different function, and don't use one to call another, because
	 * for pairwiseGlobalAlignScore, we just need to calculate score, we do not
	 * need to store the array l. we can use linear space pairwise alignment.
	 * this will save some space;
	 */
	/* Time  and Space complexity
	 * Time O(mn)
	 * Space O(2mn)
	 * */
	size_t n = str1[0].size(), m = str2[0].size();

	s.resize(n + 1);
	l.resize(n + 1);
	for (size_t i = 0; i <= n; i++) {
		s[i].resize(m + 1, 0);
		l[i].resize(m + 1, LEFT);
	}

	s[0][0] = 0;
	for (size_t j = 1; j <= m; j++) {
		s[0][j] = s[0][j - 1] + computeProfileIJ(lfPro2, J);
		l[0][j] = LEFT;
	}
	for (size_t i = 1; i <= n; i++) {
		s[i][0] = s[i - 1][0] + computeProfileIJ(lfPro1, I);
		l[i][0] = UP;
		for (size_t j = 1; j <= m; j++) {
			double s1, s2, s3;
			s1 = s[i - 1][j - 1] + computeProfileIJ(lfPro1, I, lfPro2, J);
			s2 = s[i - 1][j] + computeProfileIJ(lfPro1, I);
			s3 = s[i][j - 1] + computeProfileIJ(lfPro2, J);
			pair<double, char> charMatchResult = max(s1, s2, s3);
			s[i][j] = charMatchResult.first;
			l[i][j] = charMatchResult.second;
		}
	}
	int p = n, q = m;
	vector<string> rstr1(str1.size());
	vector<string> rstr2(str2.size());
	while (p >= 0 && q >= 0 && (p + q != 0)) { // trace back from s[n][m] to s[0][0]
		if (l[p][q] == DIAG) {
			for (size_t t = 0; t < str1.size(); t++) {
				rstr1[t].push_back(str1[t][P]);
			}
			for (size_t t = 0; t < str2.size(); t++) {
				rstr2[t].push_back(str2[t][Q]);
			}
			p = p - 1;
			q = q - 1;
		} else if (l[p][q] == UP) {
			for (size_t t = 0; t < str1.size(); t++) {
				rstr1[t].push_back(str1[t][P]);
			}
			for (size_t t = 0; t < str2.size(); t++) {
				rstr2[t].push_back('-');
			}
			p = p - 1;
		} else if (l[p][q] == LEFT) {
			for (size_t t = 0; t < str1.size(); t++) {
				rstr1[t].push_back('-');
			}
			for (size_t t = 0; t < str2.size(); t++) {
				rstr2[t].push_back(str2[t][Q]);
			}
			q = q - 1;
		}
	}
	for (size_t t = 0; t < str1.size(); t++) {
		stringReverse(rstr1[t]);
	}
	for (size_t t = 0; t < str2.size(); t++) {
		stringReverse(rstr2[t]);
	}
	str1 = rstr1;
	str2 = rstr2;
}

void CMultipleGlobalAlignment::multipleGlobalAlign()
{
	/* Time  and Space complexity
	 * Time O(r^2*(mn + 2log(r)) + r(r + mn + r * (4log(r) + mn))
	 * Space O(mn)
	 * */
	size_t r = seqs.size();
	setDist.resize(r);
	Dist.resize(r);
	for (size_t t = 0; t < Dist.size(); t++) {
		Dist[t].resize(r);
	}
	Indicator.resize(r, 1);
	cluster.resize(r);
	clusterID.resize(r);
	clusterProfile.resize(r);
	clusterProfile.resize(r);
	vector<vector<double> > lfPro;
	for (size_t i = 0; i < r; i++) {
		Indicator[i] = 1;
		cluster[i].push_back(seqs[i]);
		clusterID[i].push_back(i);
		setMulitiStrAlignProfile(seqs[i], lfPro);
		clusterProfile[i] = lfPro;
	}
	seqs.clear();
	for (size_t i = 0; i < r; i++) {
		for (size_t j = 0; j < i; j++) {
			double score = globalAlignScore(cluster[i], cluster[j],
					clusterProfile[i], clusterProfile[j]);
			Dist[i][j] = score;
			Dist[j][i] = score;
			setDist[i].insert(pair<double, int>(score, j));
			setDist[j].insert(pair<double, int>(score, i));
		}

	}
	for (size_t t = 1; t <= r - 1; t++) {
		double maxdis = -lf_inf;
		size_t maxIndex = 0;
		for (size_t i = 0; i < r; i++) {
			if (Indicator[i] && setDist[i].begin()->first > maxdis) {
				maxdis = setDist[i].begin()->first;
				maxIndex = i;
			}
		}
		size_t k1 = maxIndex;
		size_t k2 = setDist[k1].begin()->second;
		/*merge k2 to k1*/
		setDist[k1].clear();
		setDist[k2].clear();
		/*
		 * merge the strings in cluster1 and cluster2,
		 * and after this, they have the same length
		 * */
		globalAlignPath(cluster[k1], cluster[k2], clusterProfile[k1],
				clusterProfile[k2]);
		for (size_t j = 0; j < cluster[k2].size(); j++) {
			cluster[k1].push_back(cluster[k2][j]);
			clusterID[k1].push_back(clusterID[k2][j]);
		}
		Indicator[k2] = 0;
		cluster[k2].clear();
		clusterID[k2].clear();
		clusterProfile[k2].clear();
		setMulitiStrAlignProfile(cluster[k1], lfPro);
		clusterProfile[k1] = lfPro;
		for (size_t j = 0; j < r; j++) {
			if (Indicator[j] && j != k1) {
				setDist[j].erase(pair<double, int>(Dist[j][k1], k1));
				setDist[j].erase(pair<double, int>(Dist[j][k2], k2));

				Dist[j][k1] = globalAlignScore(cluster[k1], cluster[j],
						clusterProfile[k1], clusterProfile[j]);
				Dist[k1][j] = Dist[j][k1];

				setDist[j].insert(pair<double, int>(Dist[j][k1], k1));
				setDist[k1].insert(pair<double, int>(Dist[k1][j], j));
			}
		}
	}
}

void CMultipleGlobalAlignment::runMultipleGlobalAlignment()
{
	time_t rawtime;
	time(&rawtime);
	if (!readString(filePath))
		return;

	TimeStart = clock();
	multipleGlobalAlign();
	TimeEnd = clock();
	resultDisplay();
}
