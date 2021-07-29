
#include <cstdarg>
#include "G.h"
#include <limits>
#include <cmath>
#include <sys/time.h>


// Code used in the paper
// Scalable Betweenness Centrality Maximization via Sampling
// KDD 2015
// Ahmad Mahmoody, Charalampos E. Tsourakakis, Eli Upfal



void run_exact(const string, const bool, int, const int, const string);
void run_approx(const string, const bool, const int, const int, const string);
void run_cent_exact(const string, const bool, const int, const int, const string);
void run_cent_approx(const string, const bool, const int, const int, const string);


int main(int argc, char* argv[]) {

	string task, dir, seed_file;
	int k, idx;
	bool undir;


	for (int i=0; i<argc; ++i) {			
		if (! strcmp(argv[i],"-t")) {
			task = argv[i+1];
		} else if (! strcmp(argv[i],"-dir")) {
			dir = argv[i+1];
		} else if (! strcmp(argv[i],"-k")) {
			k = stoi(argv[i+1]);
		} else if (! strcmp(argv[i],"-i")) {
			idx = stoi(argv[i+1]);
		} else if (! strcmp(argv[i],"-seed_file")) {
			seed_file = argv[i+1];
		} else if (! strcmp(argv[i],"-edge")) {
			
			if        (! strcmp(argv[i+1],"undirected")) {
				undir = true;
			} else if (! strcmp(argv[i+1],"directed")) {
				undir = false;
			} else {
				cerr << "direction not mentioned.\n";
				return -1;
			}

		}
	}

	if (task == "exact") {
		run_exact(dir, undir, k, idx, seed_file);

	} else if (task == "approx") {
		run_approx(dir, undir, k, idx, seed_file);

	} else if (task == "cent_exact") {
		run_cent_exact(dir, undir, k, idx, seed_file);

	} else if (task == "cent_approx") {
		run_cent_approx(dir, undir, k, idx, seed_file);
	} 

	return 0;
}

void run_exact(const string dir, const bool undir, int k, const int idx, const string seed_file) {
	Graph g(dir, undir);
	vint_t seed;

	
	if (k == 0) {
		k = g.number_of_nodes;
	}

	const string centrality = "betweenness";
	clock_t t = clock();
	seed = g.bet_exact_set(k);	

	t = clock() - t;

	ofstream fout(seed_file);
	fout << double(t)/CLOCKS_PER_SEC;
	for (auto s : seed) {
		fout << "\t" << s;
	}
	fout << endl;
	fout.close();

}


void run_approx(const string dir, const bool undir, const int k, const int idx, const string seed_file) {
	Graph g(dir, undir);
	double eps = 0.1;

	//double SAMP = log(g.number_of_nodes) * pow(eps, -2);
	double SAMP = 2 * log(2* pow(g.number_of_nodes,3)) * pow(eps,-2);

	const string centrality = "betweenness";
	vint_t seed;
	
	clock_t t = clock();
	seed = g.bet_approx_set(k, SAMP);
	t = clock() - t;

	ofstream fout(seed_file);
	fout << double(t)/CLOCKS_PER_SEC;
	for (auto s : seed) {
		fout << "\t" << s;
	}
	fout << endl;
	fout.close();

}


void run_cent_exact(const string dir, const bool undir, const int k, const int idx, const string seed_file) {
	Graph g(dir, undir);

	ifstream fin(seed_file);	
	if(!fin.is_open())  {
		cerr << "file not found: " << seed_file;
	}

	long double cc;
	fin >> cc;
	vint_t S(k);
	for (int i=0; i < k; ++i) {
		fin >> cc;
		S[i] = cc;
	}
	fin.close();

	const string centrality = "betweenness";
	double cent;
	cent = g.bet_compute(S);

	ofstream fout(seed_file+".cent");
	fout << cent << "\n";
	fout.close();

}


void run_cent_approx(const string dir, const bool undir, const int k, const int idx, const string seed_file) {
	Graph g(dir, undir);

	ifstream fin(seed_file);	
	if(!fin.is_open())  {
		cerr << "file not found: " << seed_file;
	}

	long double cc;
	fin >> cc;
	vint_t S(k);
	for (int i=0; i < k; ++i) {
		fin >> cc;
		S[i] = cc;
	}
	fin.close();

	const string centrality = "betweenness";
	double cent;
	cent = g.bet_compute(S);

	ofstream fout(seed_file + ".cent");
	fout << cent << "\n";
	fout.close();
}
