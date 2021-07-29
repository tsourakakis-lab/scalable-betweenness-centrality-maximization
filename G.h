/*
	Graph Class
*/
#include <iostream>
#include <vector>
#include <deque>
#include <cstring>
#include <unordered_set>
#include <set>
#include <fstream>
#include <algorithm>
#include <utility>
#include <random>


using namespace std;

typedef vector<vector<vector<int>>> prev_t;
typedef vector<vector<int>>	vvint_t;
typedef vector<deque<int>>	vdint_t;
typedef vector<double>		vdoub_t;
typedef unordered_set<int> 	unset_t;
typedef deque<int> 			deint_t;
typedef vector<int>			vint_t;
typedef set<int> 			set_t;
typedef vector<set<int>> 	vsint_t;

/*
	The algorithms for computing betweenness are implemented
	 similar to NetworkX package
*/

class Graph
{
private:
	vvint_t nodes; // here we assume nodes are 0,1,...,N-1

	uniform_int_distribution<> nodes_distr;
	uniform_int_distribution<> nodes_distr_2;

	// The Greedy algorithm
	vint_t maxGreedy(const vvint_t& forward, const vvint_t& backward, 
		const int k);



	void bet_shortest_path_stats(int& source, deint_t &OS, vvint_t& P, vvint_t &subP, 
		vdoub_t& sigma, vdoub_t& subsig, vdoub_t &delta, unset_t& avoid_set);

	vdoub_t	bet_betweenness_cent(unset_t& avoid_set);

	void bet_combine_stats(int& source, vdoub_t &betweenness, deint_t &OS, vvint_t &P, 
		vvint_t &subP, vdoub_t &sigma, vdoub_t& subsig, vdoub_t &delta, 
		unset_t& avoid_set);

	void bet_random_shortest_path(int& s, int& t, 
		vvint_t & forward, vvint_t & backward, int& samp_idx);	



public:
	int number_of_nodes;
	int number_of_edges;

	Graph(const string dataset, bool undir); 
	// betweenness
	vint_t	bet_exact_set(const int k);
	vint_t  bet_approx_set(const int k, const int num_samp);	
	double	bet_compute(vint_t S);
};


