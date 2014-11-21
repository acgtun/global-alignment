#include "BandedGlobalAlignment.h"

using namespace spaceBandedGlobalAlignment;

CBandedGlobalAlignment::CBandedGlobalAlignment(const string & filePath1,
                                               const string & filePath2,
                                               const vector<int> & weight) {
  filePathU = filePath1;
  filePathV = filePath2;
  w = weight;
  rU.clear();
  rV.clear();
  alignScore = 0;
  lfIdentity = 0.0;

  TimeStart = 0.0;
  TimeEnd = 0.0;
  memory = 0.0;

  m = 0;
  n = 0;
  d = 0;
}

CBandedGlobalAlignment::~CBandedGlobalAlignment() {
  filePathU.clear();
  filePathV.clear();

  U.clear();
  V.clear();

  w.clear();

  rU.clear();
  rV.clear();

  for (int i = 0; i <= n; i++) {
    s[i].clear();
    l[i].clear();
  }
  s.clear();
  l.clear();
  L.clear();
  R.clear();

  ossGlLog.clear();
}

void CBandedGlobalAlignment::outputLog(const string & strOut) {
  ifstream fin("bandedGlobalAlignment.log");

  ostringstream strTemp;

  if (fin.good()) {
    string strVal;
    while (!fin.eof()) {
      getline(fin, strVal, '\n');
      strTemp << strVal << endl;
    }
  }

  fin.close();

  ofstream fout("bandedGlobalAlignment.log");
  fout << strOut;
  fout << strTemp.str();
  fout.close();
}

void CBandedGlobalAlignment::outPutsANDl() {
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= m; j++) {
      if (j < L[i] || j > R[i])
        printf("* ");
      else
        printf("%d ", s[i][slCOL(i, j)]);
    }
    printf("\n");
  }

  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= m; j++) {
      if (j < L[i] || j > R[i])
        printf("* ");
      else
        printf("%c ", l[i][slCOL(i, j)]);
    }
    printf("\n");
  }
}

int CBandedGlobalAlignment::readString(const string & strPath, string & str) {
  ifstream fin(strPath.c_str());
  if (!fin.good()) {
    printf("Cannot open the file: %s.\n", strPath.c_str());
    ossGlLog << "Cannot open the file: %s." << endl;
    fin.close();
    return 0;
  }
  bool vaild = true;
  string strTmp;
  str.clear();
  while (!fin.eof()) {
    getline(fin, strTmp, '\n');
    for (int i = 0; i < (int) strTmp.size(); i++) {
      if (strTmp[0] == '>')
        continue;
      if (strTmp[i] != 'A' && strTmp[i] != 'C' && strTmp[i] != 'G'
          && strTmp[i] != 'T' && strTmp[i] != 'N' && strTmp[i] != 'a'
          && strTmp[i] != 'c' && strTmp[i] != 'g' && strTmp[i] != 't'
          && strTmp[i] != 'n')

        vaild = false;
      else
        str += strTmp[i];
    }
  }
  if (vaild == false) {
    printf(
        "The DNA file %s contains characters which is not a, c, g, t or n. \
				These characters have been deleted.\n",
        strPath.c_str());

    ossGlLog << "The DNA file " << strPath
        << " contains characters which is not a, c, g, t or n. \
				These characters have been deleted."
        << endl;
  }
  fin.close();
  return 1;
}

void CBandedGlobalAlignment::stringReverse(string & str) {
  int n = str.size();
  char c;
  for (int i = 0; i < n / 2; i++) {
    c = str[i];
    str[i] = str[n - i - 1];
    str[n - i - 1] = c;
  }
}

MatchLabel CBandedGlobalAlignment::charMatch(const char & a, const char & b) {
  if (a == '-' && b == '-')
    return INDEL;
  else if (a == b)
    return MATCH;
  else
    return MISMATCH;
}

pair<int, char> CBandedGlobalAlignment::max(const int & s1, const int & s2,
                                            const int & s3) {
  /*if two of them are equal, then there are more than one optimal path*/
  if (s1 >= s2) {
    if (s1 >= s3)
      return pair<int, char>(s1, DIAG);
    else
      return pair<int, char>(s3, LEFT);
  } else {
    if (s2 >= s3)
      return pair<int, char>(s2, UP);
    else
      return pair<int, char>(s3, LEFT);
  }
}

void CBandedGlobalAlignment::outputfastaFormat(const string & str) {
  for (int t = 0; t < (int) str.size(); t++) {
    if (t % 80 == 0 && t != 0) {
      printf("\n");
      ossGlLog << endl;
    }
    printf("%c", str[t]);
    ossGlLog << str[t];
  }
  printf("\n");
  ossGlLog << endl;
}

void CBandedGlobalAlignment::outputResultString(const string & str,
                                                const int & start,
                                                const int & end) {
  for (int i = start; i <= end; i++) {
    printf("%c", str[i]);
    ossGlLog << str[i];
  }
  printf("\n");
  ossGlLog << endl;
}

void CBandedGlobalAlignment::resultDisplay() {
  printf("The two sequences are:\n");
  ossGlLog << "The two sequences are:\n" << endl;

  outputfastaFormat (U);
  printf("\n");
  outputfastaFormat (V);
  printf("\n");
  ossGlLog << endl;

  printf("The alignment of the two sequences is:\n");
  ossGlLog << "The alignment of the two sequences is:" << endl;
  string midline;
  int start = 0;
  int i = 0;
  for (i = 0; i < (int) rU.size(); i++) {
    if (i % 80 == 0 && i != 0) {
      outputResultString(rU, start, i - 1);
      outputResultString(midline, start, i - 1);
      outputResultString(rV, start, i - 1);
      printf("\n");
      ossGlLog << endl;
      start = i;
    }
    if (rU[i] == '-' || rV[i] == '-') {
      midline.push_back(' ');
    } else if (rU[i] == rV[i]) {
      midline.push_back('|');
    } else {
      midline.push_back(' ');
    }
  }
  outputResultString(rU, start, i - 1);
  outputResultString(midline, start, i - 1);
  outputResultString(rV, start, i - 1);
  printf("\n");
  ossGlLog << endl;

  printf("Alignment score:  %6d\n", alignScore);
  ossGlLog << "Alignment score: " << alignScore << endl;

  printf("%% of Identity: %.2lf%%\n", (double) lfIdentity * 100);
  ossGlLog << "% of Identity: " << (double) lfIdentity * 100 << endl;

  printf("Running time: %.2lfs\n",
         (double) (TimeEnd - TimeStart) / CLOCKS_PER_SEC);
  ossGlLog << "Running time: "
      << (double) (TimeEnd - TimeStart) / CLOCKS_PER_SEC << endl;

  printf("Memory used: %.2lf kilobytes\n", memory / 1024);
  char chr[100];
  sprintf(chr, "%.10lf", memory / 1024);
  ossGlLog << "Memory used: " << chr << " kilobytes" << endl;
  ossGlLog << "The length of the two sequences: n = " << n << " m = " << m
      << endl;
}

void CBandedGlobalAlignment::setLR(const int & k) {
  if (m >= n) {
    for (int i = 0; i <= n; i++) {
      L[i] = i - k / 2 > 0 ? i - k / 2 : 0;
      R[i] = i + d + k / 2 < m ? i + d + k / 2 : m;
    }
  } else {
    for (int i = 0; i <= n; i++) {
      L[i] = i - d - k / 2 > 0 ? i - d - k / 2 : 0;
      R[i] = i + k / 2 < m ? i + k / 2 : m;
    }
  }
}

bool CBandedGlobalAlignment::bandedGlobalAlignAlgorithm(const int & k) {
  setLR(k);
  for (int i = 0; i <= n; i++) {
    s[i].clear();
    l[i].clear();
    for (int j = 0; j < R[i] - L[i] + 1; j++) {
      s[i].push_back(0);
      l[i].push_back(DIAG);
    }
  }

  s[0][0] = 0;
  for (int j = L[0]; j <= R[0]; j++) {
    s[0][slCOL(0, j)] = j * w[INDEL];
    l[0][slCOL(0, j)] = LEFT;
  }

  for (int i = 1; i <= n; i++) {
    s[i][0] = i * w[INDEL];
    l[i][0] = UP;
    for (int j = L[i]; j <= R[i]; j++) {
      int s1 = -inf, s2 = -inf, s3 = -inf;
      int s1Col = slCOL(i - 1, j) - 1;
      int s2Col = slCOL(i - 1, j);
      int s3Col = slCOL(i, j) - 1;

      if (j - 1 >= L[i - 1] && j - 1 <= R[i - 1]) {
        s1 = s[i - 1][s1Col] + w[charMatch(U[I], V[J])];
      }
      if (j >= L[i - 1] && j <= R[i - 1]) {
        s2 = s[i - 1][s2Col] + w[INDEL];
      }
      if (j - 1 >= L[i] && j - 1 <= R[i]) {
        s3 = s[i][s3Col] + w[INDEL];
      }

      pair<int, char> charMatchResult = max(s1, s2, s3);
      s[i][slCOL(i, j)] = charMatchResult.first;
      l[i][slCOL(i, j)] = charMatchResult.second;
    }
  }

  int nDiff = 0;
  int p = n, q = m;
  rU.clear();
  rV.clear();
  while (p >= 0 && q >= 0 && (p + q != 0)) {  // trace back from s[n][m] to s[0][0]
    if (l[p][slCOL(p, q)] == DIAG) {
      rU.push_back(U[P]);
      rV.push_back(V[Q]);
      if (U[P] != V[Q])
        nDiff++;
      p = p - 1;
      q = q - 1;
    } else if (l[p][slCOL(p, q)] == UP) {
      rU.push_back(U[P]);
      rV.push_back('-');
      p = p - 1;
      nDiff++;
    } else if (l[p][slCOL(p, q)] == LEFT) {
      rU.push_back('-');
      rV.push_back(V[Q]);
      q = q - 1;
      nDiff++;
    }

    if (nDiff > d + k) {
      return false;
    }

  }

  stringReverse (rU);
  stringReverse (rV);
  alignScore = s[n][slCOL(n, m)];

  int cnt = 0;
  for (int t = 0; t < (int) rU.size(); t++) {
    if (rU[t] == rV[t])
      cnt++;
  }
  lfIdentity = (double) cnt / rV.size();
  memory = sizeof(rU[0]) * rU.size() + sizeof(rV[0]) * rV.size()
      + sizeof(U[0]) * U.size() + sizeof(V[0]) * V.size();
  for (int i = 0; i <= n; i++) {
    memory += sizeof(l[0][0]) + sizeof(s[0][0]) * ((double) R[i] - L[i] + 1);
  }
  return true;
}

void CBandedGlobalAlignment::runBandedGlobalAlignment() {
  ossGlLog << "\n------------------------------------------------------------"
      << endl;
  time_t rawtime;
  time(&rawtime);
  ossGlLog << asctime(localtime(&rawtime)) << endl;
  ossGlLog << filePathU << endl;
  ossGlLog << filePathV << endl;
  if (!readString(filePathU, U) || !readString(filePathV, V))
    return;
  n = U.size();
  m = V.size();
  d = m >= n ? m - n : n - m;
  if (n == 0 && m == 0) {
    printf("The two strings are empty. Please check the DNA fasta files.\n");
    ossGlLog << "The two strings are empty. Please check the DNA fasta files."
        << endl;
    return;
  }

  L.resize(n + 1);
  R.resize(n + 1);

  s.resize(n + 1);
  l.resize(n + 1);

  TimeStart = clock();
  for (int k = 1;; k *= 2) {
    if (bandedGlobalAlignAlgorithm(k))
      break;
  }
  TimeEnd = clock();
  resultDisplay();
  outputLog(ossGlLog.str());
}
