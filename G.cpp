#include "G.h"
#include <cmath>
#include <random>


random_device rd;
mt19937 eng(rd());
discrete_distribution<int> jump;

void sample_two(const int B, int & x, int & y) {
	uniform_int_distribution<> FIRST(0,B-1);
	uniform_int_distribution<> SECOND(0,B-2);

	x = FIRST(eng);
	y = SECOND(eng);
	if (y >= x) {
		++y;
	}
}

vint_t sort_index(const vint_t & v) {
	vector<pair<int,int>> a;
	for(int i=0; i<v.size(); ++i) {
		pair<int,int> tmp;
		tmp.first = -v[i];
		tmp.second = i;
		a.push_back(tmp);
	}
	sort(a.begin(), a.end());
	vint_t output;
	for (auto p : a) {
		output.push_back(p.second);
	}
	return output;
}

vvint_t loadGraph(const string dataset, bool undir, int & num_V, int & num_E) {
	/*
		Each row is an edge.  
		The first row of each file should have the number of nodes and edges, resp.
	*/
	ifstream fin;
	fin.open(dataset+"graph");
	if(!fin.is_open()) {
		cerr << "Cannot open: " << dataset + "graph" << endl;
		exit(1);
	}

	int u, v;
	unset_t V;
	
	fin >> num_V; // number of nodes
	fin >> num_E; // number of edges
	// cout << "V: " << num_V << ", E: " << num_E << endl; cout.flush();
	
	while (fin >> u >> v) {
		V.insert(u);
		V.insert(v);

	}
	fin.close();
	
	vector<unset_t> G(V.size());
	fin.open(dataset+"graph");
	fin >> num_V; 
	fin >> num_E;
	while (fin >> u >> v) {
		G[u].insert(v);
		if (undir) {
			G[v].insert(u);
		}
	}
	fin.close();

	// getting nodes:
	vvint_t nodes(V.size());
	
	for (int i=0; i<V.size(); ++i) {
		for (auto a : G[i]) {
			nodes[i].push_back(a);
		}
	}
	return nodes;
}

vint_t Graph::maxGreedy(const vvint_t& forward, const vvint_t& backward, const int k) {
	vint_t seed = vint_t(k);
	vint_t	 degree(number_of_nodes);
	vector<bool> seen_samp(backward.size(), false);

	for (int i=0; i<number_of_nodes; ++i) {
		degree[i] = forward[i].size();
	}

	int max_idx;
	for(int i=0; i < k; ++i) {
		auto res = max_element(degree.begin(), degree.end());
		max_idx=distance(degree.begin(), res);
		seed[i] = max_idx;		

		for (int samp_idx : forward[max_idx]) {
			if (!seen_samp[samp_idx]) {
				seen_samp[samp_idx] = true;	
				for (int u : backward[samp_idx]) {
					--degree[u];
				}
			}				
		}
		--degree[max_idx];
	}
	return seed;	
}

Graph::Graph(const string dataset, bool undir) {		
	nodes = loadGraph(dataset, undir, number_of_nodes, number_of_edges);

	vector<long double>  wedge_weights(number_of_nodes, 0);
	long double Z = 0;
	for (int i=0; i < number_of_nodes; ++i) {
		sort(nodes[i].begin(), nodes[i].end());
		wedge_weights[i] = nodes[i].size() * (nodes[i].size()-1) / 2;
		Z += wedge_weights[i];
	}
	for (int i=0; i < number_of_nodes; ++i) {
		wedge_weights[i] /= Z;
	}
	
	nodes_distr = uniform_int_distribution<>(0,number_of_nodes-1);
	nodes_distr_2 = uniform_int_distribution<>(0,number_of_nodes-2);

	cout << "Loaded -- number of nodes: " << number_of_nodes << endl;
	cout << "       -- number of edges: " << number_of_edges << endl;
}



// -----------------------------------------------------------
// -----------------------------------------------------------
// betweenness:
// -----------------------------------------------------------
// -----------------------------------------------------------

vint_t Graph::bet_approx_set(const int k, const int num_samp) {	
	int s, t; // source and terminal
	vvint_t forward(number_of_nodes);
	vvint_t backward(num_samp);

	for (int samp_idx=0; samp_idx < num_samp; ++samp_idx) {
		s = nodes_distr(eng);		
		int tmp = nodes_distr_2(eng);
		if (tmp < s) {
			t = tmp;
		} else {
			t = tmp + 1;
		}		
		bet_random_shortest_path(s, t, forward, backward, samp_idx);
	}	
	return maxGreedy(forward, backward, k);
}


vint_t	Graph::bet_exact_set(const int k){
	
	int idx = 0, max_node;	
	vint_t output;
	unset_t avoid_set;
	vdoub_t B;

	while ((idx < k) && (idx < number_of_nodes)) {
		++idx;
		B = bet_betweenness_cent(avoid_set);				
		
		auto mres = max_element(B.begin(), B.end());		
		max_node=distance(B.begin(), mres);		
				
		output.push_back(max_node);		
		avoid_set.insert(max_node);
	}	
	return output;
}


double Graph::bet_compute(vint_t S) {
	int u;
	double centrality = 0;
	vdoub_t B;
	unset_t avoid_set = unset_t();

	while (!S.empty()) {
		B = bet_betweenness_cent(avoid_set);
		u = S.back();		
		S.pop_back();
		avoid_set.insert(u);		
		centrality += B[u];		
	}
	return centrality / (number_of_nodes * (number_of_nodes-1));	
}









/*
	Functions needed for computing exact centralty
*/
vdoub_t Graph::bet_betweenness_cent(unset_t& avoid_set){ 
// computes the centrality of all nodes avoiding a set of nodes.		
	vdoub_t  betweenness(number_of_nodes,0);
	vvint_t P(number_of_nodes);
	vvint_t subP(number_of_nodes);
	deint_t OS;
	vdoub_t sigma(number_of_nodes,0);
	vdoub_t subsig(number_of_nodes,0);
	vdoub_t delta(number_of_nodes,0);
	
	for(int source=0; source< number_of_nodes; ++source) {		
		bet_shortest_path_stats(source, OS, P, subP, sigma, subsig, delta, avoid_set);			
		bet_combine_stats(source, betweenness, OS, P, subP, sigma, subsig, delta, avoid_set);		
		
		fill(P.begin(), P.end(), vector<int>());
		fill(subP.begin(), subP.end(), vector<int>());
		OS.clear(); // TODO: shouldn't be necessary
		fill(sigma.begin(), sigma.end(), 0);
		fill(subsig.begin(), subsig.end(), 0);
		fill(delta.begin(), delta.end(), 0);
	}		
	return betweenness;
}

void Graph::bet_shortest_path_stats(int& source, deint_t &OS, vvint_t& P, vvint_t &subP, 
	vdoub_t& sigma, vdoub_t& subsig, vdoub_t &delta, unset_t &avoid_set) {
	
	int w,v, Dv;		
	double sigmav, subsigv;

	vector<bool> added_Q(number_of_nodes, false);
	vector<int> D(number_of_nodes, -1);
	sigma[source] = 1;
	subsig[source] = 1;

	D[source] = 0;
	deint_t Q;

	Q.push_back(source);	
	while (!Q.empty()) {
		v = Q.front();				
		Q.pop_front();
		OS.push_back(v);
		delta[v] = 0;

		Dv = D[v];
		sigmav = sigma[v];		
		subsigv = subsig[v];		
		for (int w : nodes[v]) {						
			if (D[w] == -1) {				
				if (!avoid_set.count(v) || v==source) {
					Q.push_back(w);
					added_Q[w] = true;
				}
				D[w] = Dv+1;
			}
			if (D[w] == Dv+1){				
				sigma[w] += sigmav;
				
				if (!avoid_set.count(v) || v==source) {
					if (! added_Q[w]) {
						Q.push_back(w);
						added_Q[w] = true;	
					}

					subsig[w] += subsigv;
					subP[w].push_back(v);
				}
				P[w].push_back(v);
			}
		} 
	}	
}

void Graph::bet_combine_stats(int& source, vdoub_t &betweenness, deint_t& OS, vvint_t &P, vvint_t &subP, vdoub_t &sigma, vdoub_t& subsig, vdoub_t &delta, unset_t& avoid_set) {
	int w;
	double coeff;
	while (!OS.empty()) {		
		w = OS.back();
		OS.pop_back();
		coeff = (double) (1.0 / sigma[w] + delta[w]/subsig[w]);		
		for (int v : subP[w]) {			
			delta[v] += subsig[v] * coeff;
		}		
		if (w!= source) {
			betweenness[w] += delta[w];				
		}
	}
}



/*
	Functions needed for approximate centrality
*/
void Graph::bet_random_shortest_path(int & s, int & t, 
	vvint_t & forward, vvint_t & backward, int & samp_idx) {

	int tmp;

	vvint_t P(number_of_nodes);
	vvint_t subP(number_of_nodes);
	deint_t OS;
	vdoub_t sigma(number_of_nodes,0);
	vdoub_t subsig(number_of_nodes,0);
	vdoub_t delta(number_of_nodes,0);
	unset_t avoid_set;

	bet_shortest_path_stats(s, OS, P, subP, sigma, subsig, delta, avoid_set);

	vector<int> current_prev = P[t];

	while (!current_prev.empty()) {

		vector<int> tmp_pr(current_prev.size());
		for (int i=0; i < current_prev.size(); ++i) {
			tmp_pr[i] = sigma[current_prev[i]];
		}

		jump = discrete_distribution<int>(tmp_pr.begin(), tmp_pr.end());
		tmp = jump(eng);

		if (current_prev[tmp] != s) {
			backward[samp_idx].push_back(current_prev[tmp]);
			forward[current_prev[tmp]].push_back(samp_idx);
		}		
		current_prev = P[current_prev[tmp]];
	}	
}