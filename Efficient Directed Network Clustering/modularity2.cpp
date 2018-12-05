/*
This code uses simulated annealing to find the near-optimal clustering of an unweighted, directed network
using modularity as the objective function.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>	
#include <time.h>		//for calculating run-time
#include <string>		//for outputting file name
#include <iomanip>
#include <math.h>       //need for exponent fct.
#include <algorithm>    //need for mininimum fct.
#include <vector>
#include <cctype>		//for yes/no question
#include <tuple>

using namespace std;

//FUNCTION PROTOTYPES
int initial_pop(vector<int>& group_size, vector<int>& cluster_membership, int k, int N); //create initial membership vector
vector<int> mutate_membership(vector<int> current_membership, vector<int> group_size); //mutate membership vector
long double uniform_rand(int seed); //random (0,1) draw

int modularity2(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success) {

	//Get number of distinct edges in the network
	int E = 0; //Number of edges (in and out) in network
	for (int i = 0; i < links.size(); i++)
	{
		for (int j = 0; j < links[i].size(); j++) {
			if (get<2>(links[i][j]) == 'T') {
				E++;
			}
		}
	}
	//cout << "The number of edges in MAIN are: " << E << endl;

	int STG1_CT = 0;
	int STG2_CT = 0;
	int STG3_CT = 0;

	//Get the number of clusters in which to cluster the network
	int min_k;
	int max_k;
	std::cout << "Enter the lower bound for number of clusters: "; 	std::cin >> min_k;
	std::cout << "Enter the upper bound for number of clusters: ";	std::cin >> max_k;

	double LimitIT = 1.0e-8; //Minimum temperature. Program stops when this number is reached
	long double epsilon = 1.0e-6; //Minimum change in loglikelihood required to accept proposed membership mutation
	int TL_counter = 0; //Number of iterations at the current temperature
	int Success_counter = 0; //Number of proposed membership mutations accepted at the current temperature

	double IT = InitTemp;  //Starting temperature. Reset for every new value of k

						   //Matrices and vectors to store final cluster assignments, BIC, loglikelihood and test statistics for each k
	vector< vector <int >> MEMBERSHIP; //final cluster assignment for each k  
	vector<long double> BIC; //BIC for each k's final cluster assignment
	vector<long double> AIC; //BIC for each k's final cluster assignment
	vector<long double> LOGLIK;  //loglikelihood for each k's final cluster assignment
	vector<long double> TESTSTAT;  //test statistic for each k's final cluster assignment

	int final_k_BIC = 0;
	long double final_BIC = 0;
	long double final_LOGLIK_BIC = 0;
	long double final_TESTSTAT_BIC = 0;
	vector<int> final_membership_BIC;

	int final_k_AIC = 0;
	long double final_AIC = 0;
	long double final_LOGLIK_AIC = 0;
	long double final_TESTSTAT_AIC = 0;
	vector<int> final_membership_AIC;

	time_t SA_start = time(0);//Time until solution

							  //Repeate Simulated Annealing process for each k in the range of specified cluster numbers
	for (int k = min_k; k <= max_k; k++)
	{
		std::cout << "Now evaluating k=" << k << endl;

		int overall_iteration_count = 0;
		long double delta_q = 0; //change in modularity

								 //Initialize POS and OBS matrices (kxk vector of vectors)
		vector< vector <unsigned long long int >> OBS(k, vector<unsigned long long int>(k, 0));
		vector< vector <unsigned long long int >> POS(k, vector<unsigned long long int>(k, 0));

		//Initialize and populate cluster assignment and size vectors
		vector<int> group_size(k, 0);//used in initial pop to store number of nodes in each group 1:k
		vector<int> cluster_membership;	//holds cluster assignments of N vertices

										//vertex to be reassigned, the cluster in which it is currently assigned, and the cluster to which it will be reassigned
		vector<int> reassignment;
		int l_star;		//The vertex to be reassigned (0:[N-1])
		int clust_l;	//L-star's cluster membership before reassignment (1:k)
		int clust_j;	//L-star's cluster membership after reassignment (1:k)

						//Total weight of links to l*, by cluster (These are the Z_(to *) variables in the paper)
		vector<int> Z_to(k);
		//Total weight of links from l*, by cluster (These are the Z_(from *) variables in the paper)
		vector<int> Z_from(k);

		//Flag to allow creation of initial clustering (at the current k) that will be mutated on later runs.
		int first_run = 1;

		//Flag that indicates we should retain the current solution likelihood if we reject the proposed solution, so the current solution likelihood need not be recalculated
		int retain_current_sol = 0;

		//Reset initial temperature at the start of each k-loop
		IT = InitTemp;
		//Reset count length at the start of each k-loop
		TL_counter = 0;
		//Reset success counter at the start of each k-loop
		Success_counter = 0;

		//START Simulated Annealing Algorithm
		while (IT > LimitIT)
		{
			TL_counter++; //Iterate number of attempts at current temperature

						  //After TL attempts at reclustering at the current temperature, or Max_Success # of successes at curren temperature
						  //Slowly decrease temperature until it reaches lower limit
			if (TL_counter >= TL || Success_counter >= Max_Success)
			{
				IT = CR*IT; //decrease current temperature
				TL_counter = 0; //reset length counter
				Success_counter = 0; //reset success counter
			}

			//Generate initial clustering (1st run only). This is the initial "current clustering".
			//On runs >1, the "current clustering" will either be the one generated here, or a proposed clustering that was accepted below.
			if (first_run == 1)
			{

				//Populate cluster assignment and size vectors
				initial_pop(group_size, cluster_membership, k, N);//assigns initial cluster assignments to z, and counts cluster sizes
				first_run = 0;

				//Populate observed OBS edge weight vector. This requires 1 full traversal of adjacency list
				for (int i = 0; i < links.size(); i++)
				{
					for (int j = 0; j < links[i].size(); j++)
					{
						if (get<2>(links[i][j]) == 'T')//edge is FROM src vert links[i] TO dest vert links[i][j]...
						{
							if (cluster_membership[i] == cluster_membership[get<0>(links[i][j])])//src/dest verts belong to same cluster 
							{
								OBS[cluster_membership[i] - 1][cluster_membership[i] - 1]
									= OBS[cluster_membership[i] - 1][cluster_membership[i] - 1]
									+ get<1>(links[i][j]);
							}
							else//src/dest clusters differ
							{
								OBS[cluster_membership[i] - 1][cluster_membership[get<0>(links[i][j])] - 1]
									= OBS[cluster_membership[i] - 1][cluster_membership[get<0>(links[i][j])] - 1]
									+ get<1>(links[i][j]);
							}
						}
					}
				}

				//Populate possible edge count vector. This is done using the cluster size information
				for (int i = 0; i < k; i++)
				{
					POS[i][i] = ((unsigned long long int) group_size[i] * (group_size[i] - 1)); //2*[nchoose2]=2*[N(N-1)/2]=N*N-1
				}
				for (int i = 0; i < k - 1; i++)
				{
					for (int j = i + 1; j < k; j++)
					{
						//add number of edges possible between cluster i and cluster j
						//This is only run once, on the first run. Since each i/j combination is only considered once, can't we just shorten this to?:
						//POS[i][j] = (unsigned long long int) group_size[i] * group_size[j];
						//But keeping it as-is shouldn't hurt anything either...
						POS[i][j] = POS[i][j] + (unsigned long long int) group_size[i] * group_size[j];
						POS[j][i] = POS[i][j];
					}
				}
			}//end 1st-run WHILE loop

			 //Propose a new clustering
			 //Select Random Vertex l* for reassignment and select cluster to which it will be moved
			 //This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
			reassignment = mutate_membership(cluster_membership, group_size);
			l_star = reassignment[0];				//The vertex to be reassigned
			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
			clust_j = reassignment[1];				//l-star's cluster membership after reassignment

													//Zero out Z (to/from cluster size) vectors. Best practice for runs > 1			
			for (int i = 0; i < Z_to.size(); i++) { Z_to[i] = 0; }//Total weight of links to l*, by cluster (These are the Z_(to s) values in the paper)			
			for (int i = 0; i < Z_from.size(); i++) { Z_from[i] = 0; }//Total weight of links from l*, by cluster (These are the Z_(from s) values in the paper)

																	  //Populate to/from cluster size vectors
			for (int j = 0; j < links[l_star].size(); j++)
			{
				//get<0>(links) -Destination vertex
				//get<1>(links) -Weight
				//get<2>(links) -char in {T,F}
				if (get<2>(links[l_star][j]) == 'T')//edge extends FROM vertex l_star TO the jth vertex in l_star's adjacency list row
				{
					Z_to[cluster_membership[get<0>(links[l_star][j])] - 1] = Z_to[cluster_membership[get<0>(links[l_star][j])] - 1] + get<1>(links[l_star][j]);
				}
				else//edge extends FROM jth vertex in l_star's adjacency list row TO vertex l_star
				{
					Z_from[cluster_membership[get<0>(links[l_star][j])] - 1] = Z_from[cluster_membership[get<0>(links[l_star][j])] - 1] + get<1>(links[l_star][j]);
				}
			}
			/*
			std::cout << "Cluster memberships:" << endl;
			for (int i = 0; i < cluster_membership.size(); i++){ std::cout << cluster_membership[i] << "   "; } std::cout << endl;

			for (int i = 1; i <= k; i++){
			std::cout << endl << "The number of edges from vertex " << vert << " to cluster " << i << " is " << Z_to[i - 1] << endl;}
			for (int i = 1; i <= k; i++){
			std::cout << endl << "The number of edges from cluster " << i << " to vertex " << vert << " is " << Z_from[i - 1] << endl;}

			std::cout << endl << "Z_to():" << endl;
			for (int i = 0; i < Z_to.size(); i++){std::cout << "<" << i + 1 << ", " << Z_to[i] << "> ";}

			std::cout << endl << "Z_from():" << endl;
			for (int i = 0; i < Z_from.size(); i++){std::cout << "<" << i + 1 << ", " << Z_from[i] << "> ";}
			*/

			//NEED TO TEMPORARILY UPDATE OBS and POS elements 
			//IF we accept the proposed solution, these will become the new OBS/POS matrices. Otherwise, we'll discare and try another solution
			vector< vector <unsigned long long int >> OBS_TEMP(k, vector<unsigned long long int>(k, 0)); //Observed edge weight matrix for the proposed clustering
			vector< vector <unsigned long long int >> POS_TEMP(k, vector<unsigned long long int>(k, 0)); //Possible edge count matrix for the proposed clustering
																										 //How can I automatically populate these matrices with the contents of OBS/POS rather than first setting
																										 //them to 0 and then copying OBS/POS??
			OBS_TEMP = OBS;
			POS_TEMP = POS;

			//Update OBS and POS matrix elements for proposed clustering strictly associated ONLY with elements IN {l,j}
			OBS_TEMP[clust_l - 1][clust_l - 1] = OBS[clust_l - 1][clust_l - 1] - Z_to[clust_l - 1] - Z_from[clust_l - 1];
			OBS_TEMP[clust_j - 1][clust_j - 1] = OBS[clust_j - 1][clust_j - 1] + Z_to[clust_j - 1] + Z_from[clust_j - 1];
			OBS_TEMP[clust_l - 1][clust_j - 1] = OBS[clust_l - 1][clust_j - 1] - Z_to[clust_j - 1] + Z_from[clust_l - 1];
			OBS_TEMP[clust_j - 1][clust_l - 1] = OBS[clust_j - 1][clust_l - 1] + Z_to[clust_l - 1] - Z_from[clust_j - 1];
			//group_size[i-1] is the size of cluster i
			POS_TEMP[clust_l - 1][clust_l - 1] = (group_size[clust_l - 1] - 1)*(group_size[clust_l - 1] - 2);
			POS_TEMP[clust_j - 1][clust_j - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_j - 1]);
			POS_TEMP[clust_l - 1][clust_j - 1] = (group_size[clust_l - 1] - 1)*(group_size[clust_j - 1] + 1);
			POS_TEMP[clust_j - 1][clust_l - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_l - 1] - 1);

			//Update OBS and POS matrix elements for proposed clustering that are associated with elements {l,j} AND another element NOT IN {l,j}
			for (int star = 1; star <= k; star++)
			{
				if ((star != clust_l) && (star != clust_j))
				{
					OBS_TEMP[star - 1][clust_l - 1] = OBS[star - 1][clust_l - 1] - Z_from[star - 1];
					OBS_TEMP[clust_l - 1][star - 1] = OBS[clust_l - 1][star - 1] - Z_to[star - 1];
					OBS_TEMP[star - 1][clust_j - 1] = OBS[star - 1][clust_j - 1] + Z_from[star - 1];
					OBS_TEMP[clust_j - 1][star - 1] = OBS[clust_j - 1][star - 1] + Z_to[star - 1];

					POS_TEMP[star - 1][clust_l - 1] = (group_size[star - 1])*(group_size[clust_l - 1] - 1);
					POS_TEMP[clust_l - 1][star - 1] = (group_size[clust_l - 1] - 1)*(group_size[star - 1]);
					POS_TEMP[star - 1][clust_j - 1] = (group_size[star - 1])*(group_size[clust_j - 1] + 1);
					POS_TEMP[clust_j - 1][star - 1] = (group_size[clust_j - 1] + 1)*(group_size[star - 1]);
				}
			}


			//Calculate change in modularity for proposed clustering
			//If change-in-mod>epsilon accept. Otherwise reject and return to top to propose a new clustering
			//Change-in-mod formula taken from: https://www.researchgate.net/publication/284869815_Directed_Louvain_maximizing_modularity_in_directed_networks
			//"Directed Louvain : maximizing modularity in directed networks"
			/*
			//node: node that we are moving
			//comm: community to which we are moving node, i.e. C in the paper
			//dnodecomm: Sum of weights of edges from vertex i to cluster C
			//w_degree_out: the sum of the weights of outgoing edges from vertex i
			//w_degree_in: the sum of the weights of incoming edges from vertex i

			//totc_out: number of edges extending from vertices of cluster C,
			//and includes edges FROM vertices of cluster C
			//totc_in: number of edges coming in to vertices of cluster C,
			//and includes edges FROM vertices of cluster C
			double totc_out = tot_out[comm];
			double totc_in = tot_in[comm];
			double degc_out = w_degree_out;
			double degc_in = w_degree_in;
			double m2   = (*g).total_weight;
			double dnc  = dnodecomm;

			return (dnc/m2 - ((degc_out*totc_in + degc_in*totc_out)/(m2*m2)));

			My interpretation:
			node: l_star
			comm: clust_j
			dnc: Z_from[clust_j-1] for i=C
			degc_out = w_degree_out: (NOTE, this was old interpretation: SUM of Z_from[clust_*-1] for i=1:k)
			Note: degc_out=out-degree of vertex c

			*Z_to[i] is the sum of the edge weights FROM vertex l_star TO the set of vertices in cluster i
			*Thus SUM(Z_to[i])=the total edge weights FROM vertex l_star TO all other vertices. In the binomial
			*case, this is equivalent to the out degree of vertex l_star.

			Thus sum(Z_to[i]) in vertex c's row of the adjacency list is degc_out=w_degree_out

			degc_in = w_degree_in: (NOTE, this was old interpretation: SUM of Z_to[clust_*-1] for i=1:k)
			Note: degc_in=in-degree of vertex c

			*Z_from[i] is the sum of the edge weights FROM the set of vertices in cluster i TO vertex l_star.
			*Thus SUM(Z_from[i])=total edge weights FROM all vertices TO vertex l_star. in the binomial case,
			*this is the equivalent of the in degree of vertex l_star.

			Thus, sum(Z_from[i]) in vertex c's row of the adjacency list is degc_in=w_degree_in

			totc_in: sum OBS[*][C]
			totc_out: sum OBS[C][*]
			m2: E

			int degc_out=0;
			int degc_in=0;
			int totc_out=0;
			int totc_in=0;
			for(int i=1;i<=k;i++){degc_out=degc_out+Z_from[i-1];}
			for(int i=1;i<=k;i++){degc_in=degc_in+Z_to[i-1];}
			for(int i=1;i<=k;i++){totc_out=totc_out+OBS[i-1][clust_j-1]}
			for(int i=1;i<=k;i++){totc_in=totc_in+OBS[clust_j-1][i-1]}

			delta_q = ((long double)Z_from[clust_j -1]/E) -((degc_out*totc_in) + (degc_in*totc_out))/(E*E)
			*/
			int degc_out = 0;
			int degc_in = 0;
			int totc_out = 0;
			int totc_in = 0;
			/*NOTE:
			Z_to[i] is the sum of the edge weights FROM vertex l_star TO the set of vertices in cluster i
			Thus SUM(Z_to[i])=the total edge weights FROM vertex l_star TO all other vertices. In the binomial
			case, this is equivalent to the out degree of vertex l_star.

			Z_from[i] is the sum of the edge weights FROM the set of vertices in cluster i TO vertex l_star.
			Thus SUM(Z_from[i])=total edge weights FROM all vertices TO vertex l_star. in the binomial case,
			this is the equivalent of the in degree of vertex l_star.
			*/
			for (int i = 1; i <= k; i++) { degc_out = degc_out + Z_to[i - 1]; }/*OLD: { degc_out = degc_out + Z_from[i - 1]; }*/
			for (int i = 1; i <= k; i++) { degc_in = degc_in + Z_from[i - 1]; }/*OLD: { degc_in = degc_in + Z_to[i - 1]; }*/
			for (int i = 1; i <= k; i++) { totc_out = totc_out + OBS[i - 1][clust_j - 1]; }
			for (int i = 1; i <= k; i++) { totc_in = totc_in + OBS[clust_j - 1][i - 1]; }

			delta_q = ((long double)Z_from[clust_j - 1] / E) - ((long double)(degc_out*totc_in) + (long double)(degc_in*totc_out)) / (E*E);

			//cout << "delta_q=" << delta_q << endl;

			overall_iteration_count++;
			//long double prob_accept_modifier = (1 + log(epsilon / (long double)InitTemp) / (log((long double)CR))*TL) - (long double)overall_iteration_count -1;
			long double prob_accept_modifier = 1;
			long double uni_prob_modifier = 1;

			int STAGE = 0; //Indicates the stage when a proposed recustering was accepted or rejected.

			if (delta_q > epsilon)
			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
			 //	std::cout << "STAGE 1: " << delta_l << "..." << epsilon << endl;
				STG1_CT++;
				STAGE = 1;
				//1. We overwrite current membership with proposed membership
				//2. We overwrite current OBS/POS matrices with those of the proposed membership
				//3. We return to the top of the SA algorithm and calculate likelihood for the (new) current membership

				//current clustering <- proposed clustering
				cluster_membership[l_star] = clust_j;

				//OBS = OBS_TEMP;
				//Update OBS vector elements indexed only by l and j
				OBS[clust_l - 1][clust_l - 1] = OBS_TEMP[clust_l - 1][clust_l - 1];//OBS[clust_l - 1][clust_l - 1] - Z_to[clust_l - 1] - Z_from[clust_l - 1];
				OBS[clust_j - 1][clust_j - 1] = OBS_TEMP[clust_j - 1][clust_j - 1];//OBS[clust_j - 1][clust_j - 1] + Z_to[clust_j - 1] + Z_from[clust_j - 1];
				OBS[clust_l - 1][clust_j - 1] = OBS_TEMP[clust_l - 1][clust_j - 1];//OBS[clust_l - 1][clust_j - 1] - Z_to[clust_j - 1] + Z_from[clust_l - 1];
				OBS[clust_j - 1][clust_l - 1] = OBS_TEMP[clust_j - 1][clust_l - 1];//OBS[clust_j - 1][clust_l - 1] + Z_to[clust_l - 1] - Z_from[clust_j - 1];
																				   //Update POS vector elements indexed only by l and j
																				   //group_size[i-1] is the size of cluster i, n_i
				POS[clust_l - 1][clust_l - 1] = POS_TEMP[clust_l - 1][clust_l - 1];//(group_size[clust_l - 1] - 1)*(group_size[clust_l - 1] - 2);
				POS[clust_j - 1][clust_j - 1] = POS_TEMP[clust_j - 1][clust_j - 1];//(group_size[clust_j - 1] + 1)*(group_size[clust_j - 1]);
				POS[clust_l - 1][clust_j - 1] = POS_TEMP[clust_l - 1][clust_j - 1];// (group_size[clust_l - 1] - 1)*(group_size[clust_j - 1] + 1);
				POS[clust_j - 1][clust_l - 1] = POS_TEMP[clust_j - 1][clust_l - 1];// (group_size[clust_j - 1] + 1)*(group_size[clust_l - 1] - 1);

				//Update OBS and POS vector elements indexed by either l OR j
				for (int star = 1; star <= k; star++)
				{
					if ((star != clust_l) && (star != clust_j))
					{
						OBS[star - 1][clust_l - 1] = OBS_TEMP[star - 1][clust_l - 1];//OBS[star - 1][clust_l - 1] - Z_from[star - 1];
						OBS[clust_l - 1][star - 1] = OBS_TEMP[clust_l - 1][star - 1];//OBS[clust_l - 1][star - 1] - Z_to[star - 1];
						OBS[star - 1][clust_j - 1] = OBS_TEMP[star - 1][clust_j - 1];//OBS[star - 1][clust_j - 1] + Z_from[star - 1];
						OBS[clust_j - 1][star - 1] = OBS_TEMP[clust_j - 1][star - 1];//OBS[clust_j - 1][star - 1] + Z_to[star - 1];
																					 //group_size[i] is the size of cluster i, n_i
						POS[star - 1][clust_l - 1] = POS_TEMP[star - 1][clust_l - 1];// (group_size[star - 1])*(group_size[clust_l - 1] - 1);
						POS[clust_l - 1][star - 1] = POS_TEMP[clust_l - 1][star - 1];// (group_size[clust_l - 1] - 1)*(group_size[star - 1]);
						POS[star - 1][clust_j - 1] = POS_TEMP[star - 1][clust_j - 1];// (group_size[star - 1])*(group_size[clust_j - 1] + 1);
						POS[clust_j - 1][star - 1] = POS_TEMP[clust_j - 1][star - 1];// (group_size[clust_j - 1] + 1)*(group_size[star - 1]);
					}
				}
				//Update group size vector
				group_size[clust_l - 1]--; //cluster l loses 1 vertex
				group_size[clust_j - 1]++; //cluster j gains 1 vertex

				Success_counter++; //Iterate number of success at this temperature

			}
			else
			{
				long double uni_draw = uniform_rand(rand());
				long double min_val = min((long double)1, exp(delta_q / (IT)));
				//long double uni_draw = uni_prob_modifier*uniform_rand(rand());
				//long double min_val = min((long double)1, exp(delta_l / (prob_accept_modifier*IT)));
				//cout<<uni_draw<<"..."<<min_val<<"...";				
				//https://www.mathworks.com/help/gads/how-simulated-annealing-works.html?s_tid=gn_loc_drop
				//long double min_val =  min((long double)0.5, 1/(1+exp(-delta_l/IT)));
				//http://mathworld.wolfram.com/SimulatedAnnealing.html
				// double min_val = exp(-delta_l/IT);
				//				if (min_val >= 1) { std::cout << "delta_L=" << delta_l<<", IT="<<IT <<", delta_l/IT="<<delta_l/IT<< ", exp(delta_l/IT)="<<exp(delta_l/IT)<<", 1/(1+exp(delta_l/IT))="<< 1 / (1 + exp(delta_l / IT))<<endl; std::cin.get(); }

				if (uni_draw < min_val) //Accept proposed clustering
				{
					//	std::cout << "STAGE2: " << uni_draw << "..." << min_val << endl;
					STG2_CT++;
					STAGE = 2;
					//1. We overwrite current membership with proposed membership
					//2. We overwrite current OBS/POS matrices with those of the proposed membership
					//3. We return to the top of the SA algorithm and calculate likelihood for the (new) current membership 

					//current clustering <- proposed clustering
					cluster_membership[l_star] = clust_j;

					//OBS = OBS_TEMP;
					//Update OBS vector elements indexed only by l and j
					OBS[clust_l - 1][clust_l - 1] = OBS_TEMP[clust_l - 1][clust_l - 1];//OBS[clust_l - 1][clust_l - 1] - Z_to[clust_l - 1] - Z_from[clust_l - 1];
					OBS[clust_j - 1][clust_j - 1] = OBS_TEMP[clust_j - 1][clust_j - 1];//OBS[clust_j - 1][clust_j - 1] + Z_to[clust_j - 1] + Z_from[clust_j - 1];
					OBS[clust_l - 1][clust_j - 1] = OBS_TEMP[clust_l - 1][clust_j - 1];//OBS[clust_l - 1][clust_j - 1] - Z_to[clust_j - 1] + Z_from[clust_l - 1];
					OBS[clust_j - 1][clust_l - 1] = OBS_TEMP[clust_j - 1][clust_l - 1];//OBS[clust_j - 1][clust_l - 1] + Z_to[clust_l - 1] - Z_from[clust_j - 1];
																					   //Update POS vector elements indexed only by l and j
																					   //group_size[i-1] is the size of cluster i, n_i
					POS[clust_l - 1][clust_l - 1] = POS_TEMP[clust_l - 1][clust_l - 1];//(group_size[clust_l - 1] - 1)*(group_size[clust_l - 1] - 2);
					POS[clust_j - 1][clust_j - 1] = POS_TEMP[clust_j - 1][clust_j - 1];//(group_size[clust_j - 1] + 1)*(group_size[clust_j - 1]);
					POS[clust_l - 1][clust_j - 1] = POS_TEMP[clust_l - 1][clust_j - 1];// (group_size[clust_l - 1] - 1)*(group_size[clust_j - 1] + 1);
					POS[clust_j - 1][clust_l - 1] = POS_TEMP[clust_j - 1][clust_l - 1];// (group_size[clust_j - 1] + 1)*(group_size[clust_l - 1] - 1);

					//Update OBS and POS vector elements indexed by either l OR j
					for (int star = 1; star <= k; star++)
					{
						if ((star != clust_l) && (star != clust_j))
						{
							OBS[star - 1][clust_l - 1] = OBS_TEMP[star - 1][clust_l - 1];//OBS[star - 1][clust_l - 1] - Z_from[star - 1];
							OBS[clust_l - 1][star - 1] = OBS_TEMP[clust_l - 1][star - 1];//OBS[clust_l - 1][star - 1] - Z_to[star - 1];
							OBS[star - 1][clust_j - 1] = OBS_TEMP[star - 1][clust_j - 1];//OBS[star - 1][clust_j - 1] + Z_from[star - 1];
							OBS[clust_j - 1][star - 1] = OBS_TEMP[clust_j - 1][star - 1];//OBS[clust_j - 1][star - 1] + Z_to[star - 1];
																						 //group_size[i] is the size of cluster i, n_i
							POS[star - 1][clust_l - 1] = POS_TEMP[star - 1][clust_l - 1];// (group_size[star - 1])*(group_size[clust_l - 1] - 1);
							POS[clust_l - 1][star - 1] = POS_TEMP[clust_l - 1][star - 1];// (group_size[clust_l - 1] - 1)*(group_size[star - 1]);
							POS[star - 1][clust_j - 1] = POS_TEMP[star - 1][clust_j - 1];// (group_size[star - 1])*(group_size[clust_j - 1] + 1);
							POS[clust_j - 1][star - 1] = POS_TEMP[clust_j - 1][star - 1];// (group_size[clust_j - 1] + 1)*(group_size[star - 1]);
						}
					}

					//Update group size vector
					group_size[clust_l - 1]--; //cluster l loses 1 vertex
					group_size[clust_j - 1]++; //cluster j gains 1 vertex						

					Success_counter++; //Iterate number of success at this temperature
				}
				else //Reject proposed clustering
				{//Nothing happens here. The current cluster is retained and the process repeats itself.	
				 //	std::cout << "STAGE3: Rejected..." << endl;
					STG3_CT++;

					STAGE = 3;
					//1. We retain the current membership
					//2. We retain the current OBS/POS matrices
					//3. We retain the loglik associated with the current membership
				}
			}
		}//End SA LOOP


		 //if unweighted, all edges have weight=1 and Binomial density is appropriate
		long double sol_loglik = 0;//loglikelihood for the solution at the current value of k
		long double sol_OBS = 0;//Total number of observed edges
		long double sol_POS = 0;//Total number of possible edges
		long double LO = 0;		//loglikelihood under the null hypothesis: k=1
//long double l_density_sol = 0;
//long double j_density_sol = 0;
//long double bw_density_sol = 0;

/*
			//cluster l density
			long double param_l_sol = OBS[clust_l - 1][clust_l - 1] / (long double)POS[clust_l - 1][clust_l - 1];
			if (param_l_sol == 0 || param_l_sol == 1) { l_density_sol = 0; }
			else {
				l_density_sol = (OBS[clust_l - 1][clust_l - 1])*log(param_l_sol) + (POS[clust_l - 1][clust_l - 1] - OBS[clust_l - 1][clust_l - 1])*log(1 - param_l_sol);
			}
			//cluster j density
			long double param_j_sol = OBS[clust_j - 1][clust_j - 1] / (long double)POS[clust_j - 1][clust_j - 1];
			if (param_j_sol == 0 || param_j_sol == 1) { j_density_sol = 0; }
			else {
				j_density_sol = (OBS[clust_j - 1][clust_j - 1])*log(param_j_sol) + (POS[clust_j - 1][clust_j - 1] - OBS[clust_j - 1][clust_j - 1])*log(1 - param_j_sol);
			}
*/
			//Get  loglikelihood of change-in-loglikelihood equation for proposed clustering
			//Add logs of k within-cluster densities
			for (int row = 0; row < k; row++) {
				if (OBS[row][row] == 0 || OBS[row][row] == POS[row][row]) { sol_loglik = sol_loglik; }
				else {
					sol_loglik = sol_loglik + (OBS[row][row])*log(OBS[row][row] / (long double)POS[row][row]) + (POS[row][row] - OBS[row][row])*log(1 - (OBS[row][row] / (long double)POS[row][row]));
				}
			}

			//Just want to see if using a single density for ALL between-cluster edges improves performance
			unsigned long long BW_OBS = 0;
			unsigned long long BW_POS = 0;
			for (int i = 0; i < k; i++)
			{
				for (int j = 0; j < k; j++)
				{
					if (i != j)
					{
						BW_OBS = BW_OBS + OBS[i][j];
						BW_POS = BW_POS + POS[i][j];
					}
				}
			}
			//sol_loglik = l_density_sol + j_density_sol + ((long double)BW_OBS)*log((long double)BW_OBS / BW_POS) + ((long double)BW_POS - BW_OBS)*log(1 - (long double)BW_OBS / BW_POS);
			//Add log-density of between-cluster edges to current log-likelihood
			if (BW_OBS == 0 || BW_OBS == BW_POS) { sol_loglik = sol_loglik; }
			else {
				sol_loglik = sol_loglik + (BW_OBS)*log(BW_OBS / (long double)BW_POS) + (BW_POS - BW_OBS)*log(1 - (BW_OBS / (long double)BW_POS));
			}

			//Calculate likelihood under null model: k=1
			for (int row = 0; row < k; row++) {
				for (int col = 0; col < k; col++) {
					sol_OBS = sol_OBS + OBS[row][col];
				}
			}
			sol_POS = (long double)N*(N - 1); //2*[Nchoose2]=N(N-1) {or N*2-N, b/c no self-edges=NN-N=N(N-1)}
			LO = (sol_OBS)*log(sol_OBS / sol_POS) + (sol_POS - sol_OBS)*log(1 - sol_OBS / sol_POS);
	
		//MEMBERSHIP.push_back(cluster_membership);
		//NOTE: You can uncomment this line to store the solution membership vector for EACH value of k in [min_k, max_k].
		//However, this will consume a lot of memory for large networks and large ranges of k. Consider N=100,000 and k in [1000,10000]
		//MEMBERSHIP would then me a 100,000x9,000 element matrix.
		//Current code only stores the solution membership vector for the k value with minimum BIC
		LOGLIK.push_back(sol_loglik);
		TESTSTAT.push_back(-2 * (LO - sol_loglik)); 	// -2 * (LO - full_loglik_curr)

		//Calculate Bayesian Information Criterion for solution clustering
		//https://en.wikipedia.org/wiki/Bayesian_information_criterion
		//BIC=-2*loglik+(log n)*z
		//BIC=ln(n)z-2ln(L), where 
		//n=# of obs=# of possible edges =2*[Nchoose2] =N(N-1)
		//Or equivalently, N^2 elements in adjacency matrix - N unused diagonal elements = NN - N = N(N-1)
		//z=# of model parameters = 1 between-cluster parameters + k within-cluster parameters = k+1	
		//ln(L)=maximized loglikelihood of model = sol_loglik, which is already the maximized log-likelihood
		BIC.push_back(log(N*(N - 1))*(k + 1) - 2 * sol_loglik);
		//Calculate Akaike Information Criterion for solution clustering
		//AIC=2*z-2*ln(L)
		//https://en.wikipedia.org/wiki/Akaike_information_criterion
		AIC.push_back(2 * (k + 1) - 2 * sol_loglik);

		//Store solutions for first value of k (or ONLY value of k if min_k=max_k)
		if (k == min_k) {
			final_BIC = log(N*(N - 1))*(k + 1) - 2 * sol_loglik;
			final_AIC = 2 * (k + 1) - 2 * sol_loglik;
			final_k_BIC = k;
			final_k_AIC = k;
			final_LOGLIK_BIC = sol_loglik;
			final_TESTSTAT_BIC = -2 * (LO - sol_loglik);
			final_membership_BIC = cluster_membership;
			final_LOGLIK_AIC = sol_loglik;
			final_TESTSTAT_AIC = -2 * (LO - sol_loglik);
			final_membership_AIC = cluster_membership;
		}
		//For all k>min_k...
		//If current BIC is less than all previous BICs, store current solution as the final solution (thus far)
		if (log(N*(N - 1))*(k + 1) - 2 * sol_loglik < final_BIC) {
			final_k_BIC = k;
			final_BIC = log(N*(N - 1))*(k + 1) - 2 * sol_loglik;
			final_LOGLIK_BIC = sol_loglik;
			final_TESTSTAT_BIC = -2 * (LO - sol_loglik);
			final_membership_BIC = cluster_membership;
		}
		//If current AIC is less than all previous AICs, store current solution as the final solution (thus far)
		if (2 * (k + 1) - 2 * sol_loglik < final_AIC) {
			final_k_AIC = k;
			final_AIC = 2 * (k + 1) - 2 * sol_loglik;
			final_LOGLIK_AIC = sol_loglik;
			final_TESTSTAT_AIC = -2 * (LO - sol_loglik);
			final_membership_AIC = cluster_membership;
		}
	}//end for-k loop

	time_t SA_end = time(0);
	double SA_time = difftime(SA_end, SA_start);

	std::cout << fixed;
	std::cout << setprecision(10);
	//Output Solutions
	std::cout << "******************************************************************" << endl;
	std::cout << "Execution time = " << SA_time << " seconds" << "\n"; /*time / 60*/	/*minutes*/
	std::cout << "******************************************************************" << endl << endl;

	std::cout << "***********k values******" << endl;
	for (int k = min_k; k <= max_k; k++) { std::cout << k << " "; }
	std::cout << endl << endl;
	std::cout << "***********BIC***********" << endl;
	for (int i = 0; i < BIC.size(); i++) { std::cout << BIC[i] << " "; }
	std::cout << endl << endl;
	std::cout << "***********AIC***********" << endl;
	for (int i = 0; i < AIC.size(); i++) { std::cout << AIC[i] << " "; }
	std::cout << endl << endl;
	std::cout << "***********LOGLIK***********" << endl;
	for (int i = 0; i < LOGLIK.size(); i++) { std::cout << LOGLIK[i] << " "; }
	std::cout << endl << endl;
	std::cout << "***********TESTSTAT***********" << endl;
	for (int i = 0; i < TESTSTAT.size(); i++) { std::cout << TESTSTAT[i] << " "; }
	std::cout << endl << endl;

	std::cout << "********************************" << endl;
	std::cout << "********************************" << endl;
	std::cout << "Value for solution cluster:" << endl;
	std::cout << "k, under BIC = " << final_k_BIC << endl;
	std::cout << "k, under AIC = " << final_k_AIC << endl;
	std::cout << "BIC = " << final_BIC << endl;
	std::cout << "AIC = " << final_AIC << endl;
	std::cout << "loglik of solution, under BIC = " << final_LOGLIK_BIC << endl;
	std::cout << "loglik of solution, under AIC = " << final_LOGLIK_AIC << endl;
	std::cout << "test statistic of solution, under BIC = " << final_TESTSTAT_BIC << endl;
	std::cout << "test statistic of solution, under AIC = " << final_TESTSTAT_AIC << endl;

/*	ostringstream bestzsol;
	ofstream bestzstorage;
	// If old file exists from a previous run, delete it.
	//remove("MODULARITY_SOLUTION_CLUSTERING.txt");
	bestzsol << "MODULARITY2_SOLUTION_CLUSTERING" << ".txt";
	bestzstorage.open("MODULARITY2_SOLUTION_CLUSTERING.txt", ios::app);
	for (int i = 0; i < final_membership.size(); i++) { bestzstorage << final_membership[i] << ","; }
	bestzstorage << endl;
	bestzstorage.close();
	std::cout << "The membership vector for this clustering has been written to MODULARITY2_SOLUTION_CLUSTERING.txt." << endl;
*/
	ostringstream bestzsol_BIC;
	ofstream bestzstorage_BIC;
	// If old file exists from a previous run, delete it.
	remove("MODULARITY2_SOLUTION_CLUSTERING_BIC.txt");
	bestzsol_BIC << "MODULARITY2_SOLUTION_CLUSTERING_BIC" << ".txt";
	bestzstorage_BIC.open("MODULARITY2_SOLUTION_CLUSTERING_BIC.txt", ios::app);
	//for (int i = 0; i < final_membership_BIC.size(); i++) { bestzstorage_BIC << i << "\t" << final_membership_BIC[i] << endl; }
	for (int i = 0; i < final_membership_BIC.size(); i++) {
		bestzstorage_BIC << final_membership_BIC[i];
		if (i < final_membership_BIC.size() - 1) {
			bestzstorage_BIC << "," << " ";
		}

	}
	bestzstorage_BIC << endl;
	bestzstorage_BIC.close();
	std::cout << "The membership vector for this clustering has been written to MODULARITY2_SOLUTION_CLUSTERING_BIC.txt." << endl;

	ostringstream bestzsol_AIC;
	ofstream bestzstorage_AIC;
	// If old file exists from a previous run, delete it.
	remove("MODULARITY2_SOLUTION_CLUSTERING_AIC.txt");
	bestzsol_AIC << "MODULARITY2_SOLUTION_CLUSTERING_AIC" << ".txt";
	bestzstorage_AIC.open("MODULARITY2_SOLUTION_CLUSTERING_AIC.txt", ios::app);
	//	for (int i = 0; i < final_membership_AIC.size(); i++) { bestzstorage_AIC << i << "\t" << final_membership_AIC[i] << endl; }
	for (int i = 0; i < final_membership_AIC.size(); i++) {
		bestzstorage_AIC << final_membership_AIC[i];
		if (i < final_membership_AIC.size() - 1) {
			bestzstorage_AIC << "," << " ";
		}

	}
	bestzstorage_AIC << endl;
	bestzstorage_AIC.close();
	std::cout << "The membership vector for this clustering has been written to MODULARITY2_SOLUTION_CLUSTERING_AIC.txt." << endl;


	return 0;
}