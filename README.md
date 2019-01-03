# Likelihood-Based-Directed-Network-Clustering
This code was created to support the paper "On the statistical detection of clusters in directed networks" by Alan Ballard and Marcus B. Perry. 

Given a user-supplied directed network, this code will cluster it into a range of user-specified cluster numbers using a likelihood objective function proposed in the current paper and the directed modularity objective function proposed by Leicht and Newman (2008) in "Community structure in directed networks". 
This clustering is accomplished using simulated annealing and the cluster number range, along with the simulated annealing cooling schedule, are adjustable within the code.

A sample edge list file, hansell_el.txt, is included in this repository. This is an edge list for the directed network included by S. Hansell (1984) in "Cooperative groups, weak ties, and the integration of peer friendships". The input code for this program is not very sophisticated so users should change their network formatting to match the following requirements when inputting:

Vertex numbering must start at zero and no numbers can be skipped. For example, if there are N vertices, they should be numbered {0,1,2,....,(N-1)} with no missing integers in the set. If vertex 0<=L<=(N-1) is missing from the edge list, the program may accept the data but it will load the edges to incorrect positions and produce corrupted results.

This code is for directed networks, so each edge should appear only once in the edge list. If an edge extends from vertex 0 to vertex 1, then there will be one line on the edge list: 0 1. If an edge also extends from vertex 1 to vertex 0, then there will be another line in the edge list: 1 0. Otherwise, this line should not exist. 

The edge list should be sorted least-to-greatest first on the left column and then likewise on the right column.

This code also has the ability to generate and cluster so-called LFR benchmark directed networks proposed by Lancichinetti and Fortunato (2009) in "Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities", along with corresponding ground-truth solutions.
This clustering is accomplished as described above for user-supplied networks and the LFR network generation parameters are also adjustable within the code.







