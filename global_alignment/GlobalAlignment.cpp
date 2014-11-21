#include "GlobalAlignment.h"

using namespace spaceGlobalAlignment;

CGlobalAlignment::CGlobalAlignment(const string & filePath1,
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
  Handle = GetStdHandle(STD_OUTPUT_HANDLE);
}

CGlobalAlignment::~CGlobalAlignment() {
  filePathU.clear();
  filePathV.clear();

  U.clear();
  V.clear();

  w.clear();

  rU.clear();
  rV.clear();

  for (size_t i = 0; i <= n; i++) {
    s[i].clear();
    l[i].clear();
  }
  s.clear();
  l.clear();

  ossGlLog.clear();
}

void CGlobalAlignment::outputLog(const string & strOut) {
  ifstream fin("globalAlignment.log");

  ostringstream strTemp;

  if (fin.good()) {
    string strVal;
    while (!fin.eof()) {
      getline(fin, strVal, '\n');
      strTemp << strVal << endl;
    }
  }

  fin.close();

  ofstream fout("globalAlignment.log");
  fout << strOut;
  fout << strTemp.str();
  fout.close();
}

int CGlobalAlignment::readString(const string & strPath, string & str) {
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
    for (size_t i = 0; i < strTmp.size(); i++) {
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

void CGlobalAlignment::stringReverse(string & str) {
  int n = str.size();
  char c;
  for (int i = 0; i < n / 2; i++) {
    c = str[i];
    str[i] = str[n - i - 1];
    str[n - i - 1] = c;
  }
}

MatchLabel CGlobalAlignment::charMatch(const char & a, const char & b) {
  if (a == '-' && b == '-')
    return INDEL;
  else if (a == b)
    return MATCH;
  else
    return MISMATCH;
}

pair<int, char> CGlobalAlignment::max(const int & s1, const int & s2,
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

void CGlobalAlignment::outputfastaFormat(const string & str) {
  for (size_t t = 0; t < str.size(); t++) {
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

void CGlobalAlignment::setWindowWL() {
  _COORD coord;
  coord.X = 82;
  coord.Y = 800;

  _SMALL_RECT Rect;
  Rect.Top = 0;
  Rect.Left = 0;
  Rect.Bottom = coord.X - 1;
  Rect.Right = coord.Y - 1;

  SetConsoleScreenBufferSize(Handle, coord);  // Set Buffer Size
  SetConsoleWindowInfo(Handle, TRUE, &Rect);  // Set Window Size
}

void CGlobalAlignment::outputResultString(const string & str, const int & start,
                                          const int & end,
                                          const vector<int> & color) {
  for (int i = start; i <= end; i++) {
    SetConCol(Handle, color[i]);
    printf("%c", str[i]);
    ossGlLog << str[i];
  }
  printf("\n");
  ossGlLog << endl;
}

void CGlobalAlignment::resultDisplay() {
  setWindowWL();
  SetConCol(Handle, FOREGROUND_GREEN);
  printf("The two sequences are:\n");
  ossGlLog << "The two sequences are:" << endl;
  SetConCol(Handle, FOREGROUND_BLUE);
  outputfastaFormat (U);
  outputfastaFormat (V);
  printf("\n");
  ossGlLog << endl;

  SetConCol(Handle, FOREGROUND_GREEN);
  printf("The alignment of the two sequences is:\n");
  ossGlLog << "The alignment of the two sequences is:" << endl;
  string midline;
  vector<int> color;
  int start = 0;
  size_t i = 0;
  for (i = 0; i < rU.size(); i++) {
    if (i % 80 == 0 && i != 0) {
      outputResultString(rU, start, i - 1, color);
      outputResultString(midline, start, i - 1, color);
      outputResultString(rV, start, i - 1, color);
      printf("\n");
      ossGlLog << endl;
      start = i;
    }
    if (rU[i] == '-' || rV[i] == '-') {
      color.push_back(IndelColor);
      midline.push_back(' ');
    } else if (rU[i] == rV[i]) {
      color.push_back(matchColor);
      midline.push_back('|');
    } else {
      color.push_back(mismatchColor);
      midline.push_back(' ');
    }
  }
  outputResultString(rU, start, i - 1, color);
  outputResultString(midline, start, i - 1, color);
  outputResultString(rV, start, i - 1, color);
  printf("\n");
  ossGlLog << endl;

  SetConCol(Handle, FOREGROUND_GREEN);
  printf("Alignment score:  %6d\n", alignScore);
  ossGlLog << "Alignment score: " << alignScore << endl;
  printf("%% of Identity: %.2lf%%\n", (double) lfIdentity * 100);
  ossGlLog << "% of Identity: " << (double) lfIdentity * 100 << endl;
  printf("Running time: %.2lfs\n",
         (double) (TimeEnd - TimeStart) / CLOCKS_PER_SEC);
  ossGlLog << "Running time: "
      << (double) (TimeEnd - TimeStart) / CLOCKS_PER_SEC << endl;

  printf("Memory used: %.2lf kilobytes\n", memory / 1024);

  if (memory > pow(10.0, 9))
    ossGlLog << "Memory used: " << memory / (1024 * 1024 * 1024) << " gigabytes"
        << endl;
  else if (memory > pow(10.0, 6))
    ossGlLog << "Memory used: " << memory / (1024 * 1024) << " megabytes"
        << endl;
  else
    ossGlLog << "Memory used: " << memory / 1024 << " kilobytes" << endl;

  ossGlLog << "Memory used: " << memory / 1024 << " kilobytes" << endl;
  ossGlLog << "The length of the sequence is " << U.size() << endl;
  setBGcolor;
}

void CGlobalAlignment::globalAlignAlgorithm() {
  TimeStart = clock();
  vector<int> rint(m + 1, 0);
  vector<char> rchar(m + 1, DIAG);
  for (size_t i = 0; i <= n; i++) {
    s.push_back(rint);
    l.push_back(rchar);
  }
  rint.clear();
  rchar.clear();

  s[0][0] = 0;
  for (size_t i = 1; i <= m; i++) {
    s[0][i] = s[0][i - 1] + w[INDEL];
    l[0][i] = LEFT;
  }
  for (size_t i = 1; i <= n; i++) {
    s[i][0] = s[i - 1][0] + w[INDEL];
    l[i][0] = UP;
    for (size_t j = 1; j <= m; j++) {
      int s1 = s[i - 1][j - 1] + w[charMatch(U[I], V[J])];
      int s2 = s[i - 1][j] + w[INDEL];  //CharMatch(U[I], '-']))]
      int s3 = s[i][j - 1] + w[INDEL];  //CharMatch('-', V[J])
      pair<int, char> charMatchResult = max(s1, s2, s3);
      s[i][j] = charMatchResult.first;
      l[i][j] = charMatchResult.second;
    }
  }

  int p = n, q = m;
  while (p > 0 && q > 0) {  // trace back from s[n][m] to s[0][0]
    if (l[p][q] == DIAG) {
      rU.push_back(U[P]);
      rV.push_back(V[Q]);
      p = p - 1;
      q = q - 1;
    } else if (l[p][q] == UP) {
      rU.push_back(U[P]);
      rV.push_back('-');
      p = p - 1;
    } else if (l[p][q] == LEFT) {
      rU.push_back('-');
      rV.push_back(V[Q]);
      q = q - 1;
    }
  }
  stringReverse (rU);
  stringReverse (rV);
  alignScore = s[n][m];
  int cnt = 0;
  for (size_t t = 0; t < rU.size(); t++) {
    if (rU[t] == rV[t])
      cnt++;
  }
  lfIdentity = (double) cnt / rV.size();
  TimeEnd = clock();
  memory = sizeof(rU[0]) * rU.size() + sizeof(rV[0]) * rV.size()
      + sizeof(U[0]) * U.size() + sizeof(V[0]) * V.size()
      + (sizeof(l[0][0]) + sizeof(s[0][0])) * ((double) m + 1)
          * ((double) n + 1);
}

void CGlobalAlignment::runGlobalAlignment() {
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

  if (n == 0 && m == 0) {
    printf("The two strings are empty. Please check the DNA fasta files.\n");
    ossGlLog << "The two strings are empty. Please check the DNA fasta files."
        << endl;
    return;
  }

  globalAlignAlgorithm();
  resultDisplay();
  outputLog(ossGlLog.str());
}
