#include "Alignment.h"

int main(int argc, char *argv[])
{
	vector<int> w;	//w[0] is match, w[1] is mismatch, w[2] is indel
	w.push_back(2), w.push_back(-1), w.push_back(-2);
	vector<string> seq;
	if (argc < 3 || argc > 6) {
		cout << ErrorInfo << endl;
		return EXIT_FAILURE;
	}

	for (int t = 1; t < argc; t++) {
		if (argv[t][0] == '-') {
			int score = 0;
			if (argv[t][1] == 'm') {
				sscanf(argv[t], "-m%d", &score);
				w[MATCH] = score;
			} else if (argv[t][1] == 's') {
				sscanf(argv[t], "-s%d", &score);
				w[MISMATCH] = -score;
			} else if (argv[t][1] == 'i') {
				sscanf(argv[t], "-i%d", &score);
				w[INDEL] = -score;
			}
		} else {
			if (checkFilePath(argv[t]))
				seq.push_back(argv[t]);
			else {
				printf("Cannot open the file: %s.\n", argv[t]);
				return EXIT_FAILURE;
			}
		}
	}
	if (seq.size() != 2) {
		printf("Please input paths of TWO DNA sequence file!\n");
		cout << ErrorInfo << endl;
		return EXIT_FAILURE;
	}

	CGlobalAlignment gl(seq[0], seq[1], w);
	gl.runGlobalAlignment();
	return EXIT_SUCCESS;
}
