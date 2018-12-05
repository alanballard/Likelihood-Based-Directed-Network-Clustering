/*
This code uses simulated annealing to find the near-optimal clustering of an unweighted, directed network
using likelihood as the objective function. The likelihood for this code is defined as the product of:
1) the densities of the within-cluster edges
	--No consideration is given to direction within clusters
	--k densities in total
2) the densities of the between-cluster edges, by direction
	--There is one density for the edges extending from cluster s to cluster t,
	--and another density for the edges extending from cluster to to cluster s
	--for all clusters s,t. So, 2*kChoose2=2*(k)*(k-1) densities total

This version is calculated naively, filling the OBS/POS matrices for each proposed reclustering 
and WITHOUT using iterative methods. This is the directed-network equivalent of our original code.
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

int likelihood10(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success) {
	
	//Get the number of clusters in which to cluster the network
	int min_k;
	int max_k;
//	std::cout << "Enter the lower bound for number of clusters: "; 	std::cin >> min_k;
//	std::cout << "Enter the upper bound for number of clusters: ";	std::cin >> max_k;
min_k=5;
max_k = 10;

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

	//Vectors and variables to hold overall "best" cluster assignment and corresponding k value, BIC, loglikelihood and test statistic
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

	long double min_deltal = 0;
	long double max_deltal = 0;

	//Repeat Simulated Annealing process for each k in the range of specified cluster numbers
	for (int k = min_k; k <= max_k; k++)
	{
//		std::cout << "Now evaluating k=" << k << endl;

		long double best_likelihood = 0; //Maximum likelihood over all proposed membership (including rejected memberships)
		vector<int> best_membership; //Proposed membership corresponding to best_likelihood

		long double current_full_likelihood = 0;//Loglikelihood for current solution

		int overall_iteration_count = 0;

		//Running count for each accept/reject stage in the SA algorithm
		int STG1_CT = 0;
		int STG2_CT = 0;
		int STG3_CT = 0;

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

		long double delta_l; //change-in-loglikelihood between current clustering and proposed clustering

		//Flag to allow creation of initial clustering (at the current k) that will be mutated on later runs.
		int first_run = 1;

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
				first_run = 0; //So we don't re-enter this loop again

				//Populate observed OBS edge weight matrix. This requires 1 full traversal of adjacency list
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
						//Note that ever edge is represented twice in the adjacency lists used for this program.
						//One entry is TO ('T') vertex j from vertex i, and one entry is FROM ('F') vertex i to vertex j.
						//It would double count each edge if we considered both of these entries, so we're only counting the ones
						//marked with 'T' above.
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
						//add number of edges possible between cluster i and cluster j. This is only run once, on the first run. 
						//Since each i/j combination is only considered once, can't we just shorten this to?:
						//POS[i][j] = (unsigned long long int) group_size[i] * group_size[j]; //But keeping it as-is shouldn't hurt anything either...
						POS[i][j] = POS[i][j] + (unsigned long long int) group_size[i] * group_size[j];
						POS[j][i] = POS[i][j];
					}
				}


				//Calculate loglikelihood for current clustering, the denominator in the change-in-loglik formula
				current_full_likelihood = 0;
				//Can read observed edge weight sums in the usual manner
				//std::cout << OBS[1][2] << "..." << OBS[2][1]<<endl << endl;
				for (int row = 0; row < k; row++) {
					for (int col = 0; col < k; col++) {
						if (OBS[row][col] == 0 || OBS[row][col] == POS[row][col]) { current_full_likelihood = current_full_likelihood; }
						else {
							current_full_likelihood = current_full_likelihood + (OBS[row][col])*log(OBS[row][col] / (long double)POS[row][col]) + (POS[row][col] - OBS[row][col])*log(1 - (OBS[row][col] / (long double)POS[row][col]));
						}
					}
				}

				//Since this is the first run, take this to be the current *best* solution
				best_membership = cluster_membership;
				best_likelihood = current_full_likelihood;
			}//end 1st-run WHILE loop				


			//Propose a new clustering
			//Select Random Vertex l* for reassignment and select cluster to which it will be moved
			//This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
			reassignment = mutate_membership(cluster_membership, group_size);
			l_star = reassignment[0];				//The vertex to be reassigned
			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
			clust_j = reassignment[1];				//l-star's cluster membership after reassignment

			vector<int> proposed_membership = cluster_membership;
			proposed_membership[l_star] = clust_j;

			//NEED TO TEMPORARILY UPDATE OBS and POS elements required to calculate the denominator in the change-in-loglik formula
			//IF we accept the proposed solution, these will become the new OBS/POS matrices. Otherwise, we'll discare and try another solution
			vector< vector <unsigned long long int >> OBS_TEMP(k, vector<unsigned long long int>(k, 0)); //Observed edge weight matrix for the proposed clustering
			vector< vector <unsigned long long int >> POS_TEMP(k, vector<unsigned long long int>(k, 0)); //Possible edge count matrix for the proposed clustering
						
			vector<int> prop_group_size(k, 0);//used in initial pop to store number of nodes in each group 1:k
			prop_group_size = group_size;
			prop_group_size[clust_l - 1]--;
			prop_group_size[clust_j - 1]++;

			//Populate observed OBS edge weight matrix. This requires 1 full traversal of adjacency list
			for (int i = 0; i < links.size(); i++)
			{
				for (int j = 0; j < links[i].size(); j++)
				{
					if (get<2>(links[i][j]) == 'T')//edge is FROM src vert links[i] TO dest vert links[i][j]...
					{
						if (proposed_membership[i] == proposed_membership[get<0>(links[i][j])])//src/dest verts belong to same cluster 
						{
							OBS_TEMP[proposed_membership[i] - 1][proposed_membership[i] - 1]
								= OBS_TEMP[proposed_membership[i] - 1][proposed_membership[i] - 1]
								+ get<1>(links[i][j]);
						}
						else//src/dest clusters differ
						{
							OBS_TEMP[proposed_membership[i] - 1][proposed_membership[get<0>(links[i][j])] - 1]
								= OBS_TEMP[proposed_membership[i] - 1][proposed_membership[get<0>(links[i][j])] - 1]
								+ get<1>(links[i][j]);
						}
					}
					//Note that ever edge is represented twice in the adjacency lists used for this program.
					//One entry is TO ('T') vertex j from vertex i, and one entry is FROM ('F') vertex i to vertex j.
					//It would double count each edge if we considered both of these entries, so we're only counting the ones
					//marked with 'T' above.
				}
			}

			//Populate possible edge count vector. This is done using the cluster size information
			for (int i = 0; i < k; i++)
			{
				POS_TEMP[i][i] = ((unsigned long long int) prop_group_size[i] * (prop_group_size[i] - 1)); //2*[nchoose2]=2*[N(N-1)/2]=N*N-1
			}
			for (int i = 0; i < k - 1; i++)
			{
				for (int j = i + 1; j < k; j++)
				{
					//add number of edges possible between cluster i and cluster j. This is only run once, on the first run. 
					//Since each i/j combination is only considered once, can't we just shorten this to?:
					//POS[i][j] = (unsigned long long int) group_size[i] * group_size[j]; //But keeping it as-is shouldn't hurt anything either...
					POS_TEMP[i][j] = POS_TEMP[i][j] + (unsigned long long int) prop_group_size[i] * prop_group_size[j];
					POS_TEMP[j][i] = POS_TEMP[i][j];
				}
			}

			//Get FULL loglikelihood for proposed clustering
			long double full_prop_likelihood = 0;
			for (int row = 0; row < k; row++) {
				for (int col = 0; col < k; col++) {
					if (OBS_TEMP[row][col] == 0 || OBS_TEMP[row][col] == POS_TEMP[row][col]) { full_prop_likelihood = full_prop_likelihood; }
					else {
						full_prop_likelihood = full_prop_likelihood + (OBS_TEMP[row][col])*log(OBS_TEMP[row][col] / (long double)POS_TEMP[row][col]) + (POS_TEMP[row][col] - OBS_TEMP[row][col])*log(1 - (OBS_TEMP[row][col] / (long double)POS_TEMP[row][col]));
					}
				}
			}
			
			//If the loglikelihood of the proposed solution is better than any previously seen, retain proposed solution as the *best* solution
			if (full_prop_likelihood>best_likelihood)
			{
				best_likelihood = full_prop_likelihood;
				best_membership = proposed_membership;
			}

			//Calculate change in log-likelihood due to reclustering
			//log(new/old)=log(new)-log(old)

			//If all densities in either the new or the old loglik were discarded (and thus set to zero),
			//it is not possible to evaluate a change-in-loglik
			//set entire change-in-loglik to zero so a new clustering will be considered
			//if (current_full_likelihood == 0 || full_prop_likelihood == 0) { delta_l = 0; }
			//else { delta_l = full_prop_likelihood - current_full_likelihood; }
			delta_l = full_prop_likelihood - current_full_likelihood;
			//std::cout << "The change in log-lik is: " << delta_l << endl;	
			//if delta_l>0, then proposed loglik>current loglik, and the reclustered membership is more "likely"
			//Otherwise, proposed loglik<current loglik and the current membership is more "likely"
//if (delta_l < min_deltal) { min_deltal = delta_l; }
//if (delta_l > max_deltal) { max_deltal = delta_l; }
//cout << "min=" << min_deltal 
//<< ", prop_loglik="<<full_prop_likelihood 
//<<", curr_loglik="<< current_full_likelihood 
//<< ", delta_l=" << delta_l << ", max=" << max_deltal << endl;
//cin.get();

			overall_iteration_count++;

			//long double prob_accept_modifier = (1 + log(epsilon / (long double)InitTemp) / (log((long double)CR))*TL) - (long double)overall_iteration_count -1;
			long double prob_accept_modifier = 1;
			long double uni_prob_modifier = 1;

			int STAGE = 0; //Indicates the stage when a proposed recustering was accepted or rejected.

			if (delta_l > epsilon)
			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
			 //	std::cout << "STAGE 1: " << delta_l << "..." << epsilon << endl;
				STG1_CT++;
				STAGE = 1;

				//current likelihood <- proposed likelihood
				current_full_likelihood = full_prop_likelihood;

				//current clustering <- proposed clustering
				cluster_membership = proposed_membership;

				//current OBS/POS <- proposed OBS/POS
				OBS = OBS_TEMP;
				POS = POS_TEMP;
				//current group-size <- proposed group-size
				group_size = prop_group_size;
				//Equivalent to:
				//group_size[clust_l - 1]--; //cluster l loses 1 vertex
				//group_size[clust_j - 1]++; //cluster j gains 1 vertex

				Success_counter++; //Iterate number of success at this temperature

			}
			else
			{
				long double uni_draw = uniform_rand(rand());
				long double min_val = min((long double)1, exp(delta_l / (IT)));
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

					//current clustering <- proposed clustering
					cluster_membership = proposed_membership;

					//current likelihood <- proposed likelihood
					current_full_likelihood = full_prop_likelihood;

					//current OBS/POS <- proposed OBS/POS
					OBS = OBS_TEMP;
					POS = POS_TEMP;
					//current group-size <- proposed group-size
					group_size = prop_group_size;
					//Equivalent to:
					//group_size[clust_l - 1]--; //cluster l loses 1 vertex
					//group_size[clust_j - 1]++; //cluster j gains 1 vertex

					Success_counter++; //Iterate number of success at this temperature
				}
				else //Reject proposed clustering
				{//Nothing happens here. The current cluster is retained and the process repeats itself.	
				 //	std::cout << "STAGE3: Rejected..." << endl;
					STG3_CT++;
					STAGE = 3;
				}

			}

		//std::cout << STG1_CT << "..." << STG2_CT << "..." << STG3_CT << "..."<< overall_iteration_count<<"..."<<prob_accept_modifier<<"..."<<STAGE<<"..."<<sol_loglik_temp <<"..."<<LO_temp <<endl;

		}//End SA LOOP
		 //At this point, we have:
		 //cluster_membership, which is the solution arrived at by the SA algorithm, and associated current_full_likelihood
		 //best_membership, which is the highest-likelihood solution of all evaluated during the course of the SA algorithm, and associated best_likelihood

		
		long double sol_loglik = 0;//loglikelihood for the final solution at the current value of k
		long double sol_OBS = 0;//Total number of observed edges
		long double sol_POS = 0;//Total number of possible edges
		long double LO = 0;		//loglikelihood under the null hypothesis: k=1	

		//Since the SA algorithm occassionally(?) allows us to accept worse positions than our current one (see Stage 2 above), it is 
		//possible for the solution arrived at by the SA algorithm to be worse than the "best" of all solutions evaluated.

		//If, during the course of the SA algorithm, we found a solution superior to the final solution only to
		//discard it, take the discarded solution as the final solution.
		if (best_likelihood > current_full_likelihood) {
			current_full_likelihood = best_likelihood;
			cluster_membership = best_membership;
//			cout << "Overwriting best solution for k=" << k << endl;
		}
		//Now cluster_membership is the best solution found of all evaluated, and has associated current_full_likelihood.

		//loglikelihood of final solution in the simulated annealing algorithm
		sol_loglik = current_full_likelihood;

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
		//z=# of model parameters = 2*kchoose2 between-cluster parameters + k within-cluster parameters 
		//=2*(k)(k-1)/2 + k = k(k-1)+k
		//ln(L)=maximized loglikelihood of model = sol_loglik, which is already the maximized log-likelihood
		BIC.push_back(log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik);
		//Calculate Akaike Information Criterion for solution clustering
		//AIC=2*z-2*ln(L)
		//https://en.wikipedia.org/wiki/Akaike_information_criterion
		AIC.push_back(2 * (k*(k - 1) + k) - 2 * sol_loglik); 

		//Store solutions for first value of k (or ONLY value of k if min_k=max_k)
		if (k == min_k) {
			final_BIC = log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik;
			final_AIC = 2 * (k*(k - 1) + k) - 2 * sol_loglik;
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
		if (log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik < final_BIC) {
			final_k_BIC = k;
			final_BIC = log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik;
			final_LOGLIK_BIC = sol_loglik;
			final_TESTSTAT_BIC = -2 * (LO - sol_loglik);
			final_membership_BIC = cluster_membership;
		}
		//If current AIC is less than all previous AICs, store current solution as the final solution (thus far)
		if (2 * (k*(k - 1) + k) - 2 * sol_loglik < final_AIC) {
			final_k_AIC = k;
			final_AIC = 2 * (k*(k - 1) + k) - 2 * sol_loglik;
			final_LOGLIK_AIC = sol_loglik;
			final_TESTSTAT_AIC = -2 * (LO - sol_loglik);
			final_membership_AIC = cluster_membership;
		}
	}//end for-k loop

	time_t SA_end = time(0);
	double SA_time = difftime(SA_end, SA_start);

/*
	std::cout << fixed;
	std::cout << setprecision(10);
	//Output Solutions

	std::cout << "******************************************************************" << endl;
	std::cout << "Execution time = " << SA_time << " seconds" << "\n"; //time / 60=minutes
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
*/

	ostringstream stats_name;
	ofstream stats_file;
	stats_name << "LIKELIHOOD10_STATS" << ".txt";
	stats_file.open("LIKELIHOOD10_STATS.txt", ios::app);
	stats_file << final_k_BIC << "," << final_k_AIC << "," << final_BIC << "," << final_AIC << "," << final_LOGLIK_BIC << "," << final_LOGLIK_AIC << "," << final_TESTSTAT_BIC << "," << final_TESTSTAT_AIC;
	stats_file << endl;
	stats_file.close();


	ostringstream bestzsol_BIC;
	ofstream bestzstorage_BIC;
	// If old file exists from a previous run, delete it.
//remove("LIKELIHOOD10_SOLUTION_CLUSTERING_BIC.txt");
	bestzsol_BIC << "LIKELIHOOD10_SOLUTION_CLUSTERING_BIC" << ".txt";
	bestzstorage_BIC.open("LIKELIHOOD10_SOLUTION_CLUSTERING_BIC.txt", ios::app);
	//for (int i = 0; i < final_membership_BIC.size(); i++) { bestzstorage_BIC << i << "\t" << final_membership_BIC[i] << endl; }
	for (int i = 0; i < final_membership_BIC.size(); i++) {
		bestzstorage_BIC << final_membership_BIC[i];
		if (i < final_membership_BIC.size() - 1) {
			bestzstorage_BIC << "," << " ";
		}

	}
	bestzstorage_BIC << endl;
	bestzstorage_BIC.close();
//	std::cout << "The membership vector for this clustering has been written to LIKELIHOOD10_SOLUTION_CLUSTERING_BIC.txt." << endl;

	ostringstream bestzsol_AIC;
	ofstream bestzstorage_AIC;
	// If old file exists from a previous run, delete it.
//remove("LIKELIHOOD10_SOLUTION_CLUSTERING_AIC.txt");
	bestzsol_AIC << "LIKELIHOOD10_SOLUTION_CLUSTERING_AIC" << ".txt";
	bestzstorage_AIC.open("LIKELIHOOD10_SOLUTION_CLUSTERING_AIC.txt", ios::app);
	//	for (int i = 0; i < final_membership_AIC.size(); i++) { bestzstorage_AIC << i << "\t" << final_membership_AIC[i] << endl; }
	for (int i = 0; i < final_membership_AIC.size(); i++) {
		bestzstorage_AIC << final_membership_AIC[i];
		if (i < final_membership_AIC.size() - 1) {
			bestzstorage_AIC << "," << " ";
		}

	}
	bestzstorage_AIC << endl;
	bestzstorage_AIC.close();
//	std::cout << "The membership vector for this clustering has been written to LIKELIHOOD10_SOLUTION_CLUSTERING_AIC.txt." << endl;

	return 0;
}