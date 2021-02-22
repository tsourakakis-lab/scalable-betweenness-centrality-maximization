# Readme


This repo provides the C++ implementation of the **Hedge** algorithm proposed in the KDD 2015 paper "Scalable Betweenness Centrality Maximization via Sampling" by  A. Mahmoody, C.E Tsourakakis, E. Upfal. To compile 

- g++ -std=c++11 -Ofast hedge.cpp G.cpp -o hedge



### Input

Each graph should be stored in a folder with its name, and a file named "graph". The first line of the file should be the number of nodes and edges, resp. The other lines, should represent each edge.  The graph may be directed or undirected. 



### Instructions 



*Flags*
-t : task. is a string that can be evaluated by one of the following:
"exact" for running the exhaustive algorithm and find the set of seeds of size k
"approx" same as exact but runs the hedge sampling algorithm for finding the seeds  --> main algorithm Hedge. You can choose how many samples you want by modifying the code.

"cent_exact" computing the exact centrality for a given set of seeds (not scalable)
"cent_approx" computing the centrality for a given set of seeds via sampling (scalable) 
--> For cent_exact and cent_approx the centrality is normalized by 1/n*(n-1)
--> Also, the results will be written in a file named "<seed_file> + .cent"


-dir : directory that contains the graph file, e.g., can be passed "data/CA-GrQc/" and the file with name "graph" is expected to be found in that directory.

-k : number of seeds

-i : just an index for your experiment. you can pass any integer

-seed_file: for tasks "exact" and "approx" it is going to be the output file that stores the k seeds followed by the time in seconds it took for generating that set of seeds.
For "cent_exact" and "cent_approx" it is going to be the input file that keeps the set of seeds.

-edge: direction of the edges. It can be "directed" or "undirected"



### Demos

1. ./hedge -t approx -k 5 -i 0 -dir data/CA-GrQc/ -seed_file data/CA-GrQc/seed_5_0_approx -e undirected
2. ./hedge -t exact -k 5 -i 1 -dir data/CA-GrQc/ -seed_file data/CA-GrQc/seed_5_1_approx -e undirected
3. ./hedge -t cent_exact -k 5 -i 2 -dir data/CA-GrQc/ -seed_file data/CA-GrQc/seed_5_1_approx -e undirected
