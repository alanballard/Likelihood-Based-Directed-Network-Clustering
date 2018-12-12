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

int modularity(int network_key, int N, int nb_links, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success, int min_k, int max_k, int k_int) {

	//Get number of distinct edges in the network
	int E = nb_links;
	
	vector<long double> k_i_in(N, 0); //vertex i's in-degree
	vector<long double> k_j_out(N, 0); //vertex j's out-degree
	//Note: It's a little confusing but the k_i_in and k_j_out notation are what's used in "Community Structure in Directed Networks", Leicht and Newman, 2007

	for (int s = 0; s < links.size(); s++)
	{
		for (int t = 0; t < links[s].size(); t++)
		{
			if (get<2>(links[s][t]) == 'F') { k_i_in[s]++; }//if edge extends from t to s, increase s's in degree by 1
			else { k_j_out[s]++; }//else, edge extends from s to t, increase s's out degree by one
		}
	}
	//pg. 135 of "Networks: An Introduction"
	//# of edges in network = sum(in-degrees across all vertices) = sum(out-degrees across all vertices)
	//int in_deg_sum = 0;
	//for (int i = 0; i < N; i++) { in_deg_sum = in_deg_sum + k_i_in[i]; }
	//int out_deg_sum = 0;
	//for (int i = 0; i < N; i++) { out_deg_sum = out_deg_sum + k_i_in[i]; }
	//cout << "The total number of edges is: " << E << endl;
	//cout << "The total in degrees are: " << in_deg_sum << endl;
	//cout << "The total out degrees are: " << out_deg_sum << endl;

	double LimitIT = 1.0e-8; //Minimum temperature. Program stops when this number is reached
	long double epsilon = 1.0e-6; //Minimum change in loglikelihood required to accept proposed membership mutation
	int TL_counter = 0; //Number of iterations at the current temperature
	int Success_counter = 0; //Number of proposed membership mutations accepted at the current temperature

	double IT = InitTemp;  //Starting temperature. Reset for every new value of k
	   
	int final_k = 0;
	long double final_BIC = 0;
	long double final_LOGLIK = 0;
	long double final_TESTSTAT = 0;
	vector<int> final_membership;
	long double final_MODULARITY = 0;
	long double final_TIME = 0;

	int final_iteration_count = 0;
	
	//Repeat Simulated Annealing process for each k in the range of specified cluster numbers
	for (int k = min_k; k <= max_k; k = k + k_int)
	{
		clock_t start_new = clock(); //Time until solution

		long double current_full_modularity = 0;//modularity for current solution
		long double full_prop_modularity = 0; //modularity for proposed solution

		int overall_iteration_count = 0;
		long double delta_q = 0; //change in modularity

//Step: Initialize k+1 current/proposed OBS and POS vectors, and current/proposed membership vector
		vector<unsigned long long int> OBS(k + 1, 0);//elements 0 thru k-1 are within cluster values, element k is between-cluster value
		vector<unsigned long long int> POS(k + 1, 0);//elements 0 thru k-1 are within cluster values, element k is between-cluster value
		//Data types available: http://www.cplusplus.com/doc/tutorial/variables/

		//Initialize and populate cluster assignment and size vectors
		vector<int> group_size(k, 0);//used in initial pop to store number of nodes in each group 1:k
		vector<int> cluster_membership;	//holds cluster assignments of N vertices

		//vertex to be reassigned, the cluster in which it is currently assigned, and the cluster to which it will be reassigned
		vector<int> reassignment;
		int l_star;		//The vertex to be reassigned (0:[N-1])
		int clust_l;	//L-star's cluster membership before reassignment (1:k)
		int clust_j;	//L-star's cluster membership after reassignment (1:k)

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

			}//end 1st-run WHILE loop

			 //Propose a new clustering
			 //Select Random Vertex l* for reassignment and select cluster to which it will be moved
			 //This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
			reassignment = mutate_membership(cluster_membership, group_size);
			l_star = reassignment[0];				//The vertex to be reassigned
			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
			clust_j = reassignment[1];				//l-star's cluster membership after reassignment

			//Then the proposed membership is the same as the current membership, 
			//but vertex l_star's cluster number is changed from clust_l to clust_j
			vector<int> proposed_membership = cluster_membership;
			proposed_membership[l_star] = clust_j;

			if (first_run == 1)
			{
				first_run = 0;

				//Get FULL modularity for current clustering	
				//This is the initial clustering and is only run once.
				//Thus the clustering and associated likelihood are initially assigned to the best_modularity
				//and best_membership variables, which will be updated as we propose future clusterings				

				//Have to traverse entire adjacency list...
				//Formula taken from "Community Structure in Directed Networks", Leicht and Newman, 2007
				//(1/2m)*sum[A_ij - (k_i^in*k_j^out)/2m]
				//A_ij		=1 if edge exists FROM j TO i	= where get<2>(links[i][j]) == 'F'
				//k_i^in	=in degree of vertex i			=k_i_in[i]
				//k_j^out	=out degree of vertex j			=k_j_out[j]
				//m			=# of edges in network			=E
				for (int i = 0; i < links.size(); i++)
				{
					for (int j = 0; j < links[i].size(); j++)
					{
						if (cluster_membership[i] == cluster_membership[get<0>(links[i][j])])
							//i/j belong to same cluster -> this means Kronecker delta function=1
						{
							if (get<2>(links[i][j]) == 'F') //edge extends from j to i -> this means A_ij exists
							{

								current_full_modularity =
									current_full_modularity +
									get<1>(links[i][j])//A_ij
									- ((long double)(k_i_in[i] * k_j_out[j]) / (2 * E));
							}
							else
							{
								current_full_modularity =
									current_full_modularity +
									0//A_ij
									- ((long double)(k_i_in[i] * k_j_out[j]) / (2 * E));
							}
						}

					}
				}
				current_full_modularity = current_full_modularity / (2 * E);

			}//end 1st-run WHILE loop

			 //Get FULL modularity for proposed clustering
			full_prop_modularity = 0;
			for (int i = 0; i < links.size(); i++)
			{
				for (int j = 0; j < links[i].size(); j++)
				{
					if (proposed_membership[i] == proposed_membership[get<0>(links[i][j])])
						//i/j belong to same cluster -> this means Kronecker delta function=1
					{
						if (get<2>(links[i][j]) == 'F') //edge extends from j to i -> this means A_ij exists
						{

							full_prop_modularity =
								full_prop_modularity +
								get<1>(links[i][j])//A_ij
								- ((long double)(k_i_in[i] * k_j_out[j]) / (2 * E));
						}
						else //otherwise A_ij doesn't exist, and the edge weight element is set to zero
						{
							full_prop_modularity =
								full_prop_modularity +
								0//A_ij
								- ((long double)(k_i_in[i] * k_j_out[j]) / (2 * E));
						}
					}

				}
			}
			full_prop_modularity = full_prop_modularity / (2 * E);

			//Calculate change in modularity due to reclustering
			//modularity(new)-modularity(old)
			delta_q = full_prop_modularity - current_full_modularity;
			//std::cout << "The change in modularity is: " << delta_q << endl;	
			//if delta_q>0, then proposed modularity>current modularity, and the reclustered membership is better
			//Otherwise, proposed modularity<current modularity and the current membership is better

			overall_iteration_count++;
			//long double prob_accept_modifier = (1 + log(epsilon / (long double)InitTemp) / (log((long double)CR))*TL) - (long double)overall_iteration_count -1;
			long double prob_accept_modifier = 1;
			long double uni_prob_modifier = 1;

			int STAGE = 0; //Indicates the stage when a proposed recustering was accepted or rejected.

			if (delta_q > epsilon)
			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
				STAGE = 1;

				current_full_modularity = full_prop_modularity;

				//current clustering <- proposed clustering
				cluster_membership = proposed_membership;

				//Update group size vector
				group_size[clust_l - 1]--; //cluster l loses 1 vertex
				group_size[clust_j - 1]++; //cluster j gains 1 vertex

				Success_counter++; //Iterate number of success at this temperature

			}
			else
			{
				long double one = 1;      //used in MIN call to match data types.
				double uni_draw = rand() / double(RAND_MAX);
				double min_val = min(one, exp(delta_q / IT));

				if (uni_draw < min_val) //Accept proposed clustering
				{
					STAGE = 2;

					current_full_modularity = full_prop_modularity;

					//current clustering <- proposed clustering
					cluster_membership = proposed_membership;

					//Update group size vector
					group_size[clust_l - 1]--; //cluster l loses 1 vertex
					group_size[clust_j - 1]++; //cluster j gains 1 vertex						

					Success_counter++; //Iterate number of success at this temperature
				}
				else //Reject proposed clustering
				{//Nothing happens here. The current cluster is retained and the process repeats itself.	

					STAGE = 3;
					//cluster_membership remains unchanged
					//current_full_modularity remains unchanged, and is associated with cluster_membership
				}
			}

		}//End SA LOOP
		 //At this point, we have cluster_membership, which is the solution arrived at by the SA algorithm, and associated current_full_modularity

//CALCULATE TIME REQUIRED TO REACH SOLUTION
		//cout << fixed;
		cout << setprecision(10);
		clock_t end_new = clock();
		double time = (double)(end_new - start_new) / CLOCKS_PER_SEC;

//Step: Populate OBS and POS vectors based on this membership
		//Populate observed OBS edge weight matrix. This requires 1 full traversal of adjacency list
		for (int i = 0; i < links.size(); i++)
		{
			for (int j = 0; j < links[i].size(); j++)
			{
				//Each edge is represented both ways in the adjacency list used in this program. Example, consider the edge from vertex i TO vertex j:
				//In vertex i's row, links[i][j]) == 'T' but in vertex j's row: links[j][i]) == 'F'.
				//Thus the same edge is represented twice. I'm limiting to 'T' (rather than 'F') to avoid double counting of edges
				if (get<2>(links[i][j]) == 'T') //edge is FROM src vert links[i] TO dest vert links[i][j]...
				{
					//If the src and destination vertices belong to the same cluster
					if (cluster_membership[i] == cluster_membership[get<0>(links[i][j])])
					{
						OBS[cluster_membership[i] - 1] = OBS[cluster_membership[i] - 1] + get<1>(links[i][j]);
					}
					//Otherwise they are in different clusters
					else
					{
						OBS[k] = OBS[k] + get<1>(links[i][j]);
					}
				}
			}
		}

//Don't need to divide by 2 in the directed case as the way data is read into OBS is different than the undirected case
						//Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
						//over the adjacency list, the total observed count needs to be divided by 2.
						//for (int i = 0; i < k + 1; i++) { OBS[i] = OBS[i] / 2; }

						//Populate possible edge count vector. This is done using the cluster size information
		for (int i = 0; i < k; i++)
		{
			POS[i] = ((unsigned long long int)group_size[i] * (group_size[i] - 1)); //2*[nchoose2]=2*[N(N-1)/2]=N*N-1
		}
		for (int i = 0; i < k - 1; i++)
		{
			for (int j = i + 1; j < k; j++)
			{
				//add number of edges possible between cluster i and cluster j (in both directions). This is only run once, on the first run.
				POS[k] = POS[k]
					//Adding (possible # of edges from cluster i to cluster j) AND (possible # of edges from cluster j to cluster i), which are equal values
					+ 2 * ((unsigned long long int)group_size[i] * group_size[j])
					;
			}
		}

		long double sol_loglik = 0;//loglikelihood for the solution at the current value of k
		long double sol_OBS = 0;//Total number of observed edges
		long double sol_POS = 0;//Total number of possible edges
		long double LO = 0;		//loglikelihood under the null hypothesis: k=1
		long double BIC = 0;				//Bayesian Information Criterion
		long double TESTSTAT = 0;

//STEP: USING OBS and POS vectors, calculate full loglikelihood of "best" solutions
		double l_0_w = 0; //loglikelihood of within-clusters
		double l_0_b = 0; //loglikelihood of between-clusters
		//Get loglikelihood of within-cluster edges
		for (int i = 0; i < k; i++) {
			if (OBS[i] != 0 && (OBS[i] != POS[i])) {
				l_0_w = l_0_w +
					(
					(OBS[i])*log(OBS[i] / (long double)POS[i]) + (POS[i] - OBS[i])*log(1 - OBS[i] / (long double)POS[i])
						);
			}
		}

		//Get loglikelihood of between-cluster edges
		unsigned long long int TOT_OBS = 0;
		unsigned long long int TOT_POS = 0;
		TOT_OBS = OBS[k];
		TOT_POS = POS[k];

		//Loglikelihood of between cluster edges when modeled with a single parameter
		l_0_b = (TOT_OBS)*log(TOT_OBS / (long double)TOT_POS) + ((long double)TOT_POS - TOT_OBS)*log(1 - (TOT_OBS / (long double)TOT_POS));
		sol_loglik = l_0_w + l_0_b;

//STEP: Calculate likelihood under null hypothesis of k=1 cluster
		long long int LO_OBS = 0;
		long long int LO_POS = 0;
		for (int i = 0; i <= k; i++)
		{
			LO_OBS = LO_OBS + OBS[i];
		}
		LO_POS = N * (N - 1); //2*(N choose 2)
		LO = (LO_OBS)*log(LO_OBS / (long double)LO_POS) + (LO_POS - LO_OBS)*log(1 - (LO_OBS / (long double)LO_POS));

/*
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
*/

//STEP: Calculate Bayesian Information Criterion for solution clustering
		//BIC=-2*loglik+(log n)*(#of parameters)
		//number of parameters=(k+1), number of observations=2*Nchoose2=2*N(N-1)/2=N(N-1), maximized loglikelihood=sol_loglik		
		BIC = -2 * (sol_loglik)+(k + 1)*log(N*(N - 1));
		TESTSTAT = -2 * (LO - sol_loglik);

		//Store solutions for first value of k (or ONLY value of k if min_k=max_k)
		if (k == min_k) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = sol_loglik;
			final_TESTSTAT = -2 * (LO - sol_loglik);
			final_membership = cluster_membership;
			final_MODULARITY = current_full_modularity;
			final_TIME = time;
		}
		//For all k>min_k...
		//If current modularity is greater than all previous modularities, store current solution as the final solution (thus far)
		if (current_full_modularity > final_MODULARITY) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = sol_loglik;
			final_TESTSTAT = -2 * (LO - sol_loglik);
			final_membership = cluster_membership;
			final_MODULARITY = current_full_modularity;
			final_TIME = time;
		}

//STEP: Output solution clustering statistics at each value of k		
		ostringstream atkstring;
		ofstream atkinfo;
		// If old file exists from a previous run, delete it.
		//remove("completion_modularity_all_sols.txt");
		atkstring << "completion_modularity_all_sols" << ".txt";
		atkinfo.open("completion_modularity_all_sols.txt", ios::app);
		atkinfo << network_key << "," << N << "," << nb_links << "," << k << "," << current_full_modularity << "," << sol_loglik << "," << BIC << "," << TESTSTAT << "," << time << endl;
		atkinfo.close();


		//cout << "**********************************" << endl;
	}//end for-k loop

	//STEP: Output best solution across all k
	ostringstream atkstring2;
	ofstream atkinfo2;
	// If old file exists from a previous run, delete it.
	//remove("completion_modularity_final_sol.txt");
	atkstring2 << "completion_modularity_final_sol" << ".txt";
	atkinfo2.open("completion_modularity_final_sol.txt", ios::app);
	atkinfo2 << network_key << "," << N << "," << nb_links << "," << final_k << "," << final_MODULARITY << "," << final_LOGLIK << "," << final_BIC << "," << final_TESTSTAT << "," << final_TIME << endl;
	atkinfo2.close();

	//STEP: Output best solution clustering
	ostringstream bestzsol;
	ofstream bestzstorage;
	// If old file exists from a previous run, delete it.
	//remove("solution_modularity.txt");
	bestzsol << "solution_modularity" << ".txt";
	bestzstorage.open("solution_modularity.txt", ios::app);
	for (int i = 0; i < N; i++) {
		bestzstorage << final_membership[i];
		if (i < N - 1) { bestzstorage << ", "; }
	}
	bestzstorage << endl;
	bestzstorage.close();

	return 0;
}