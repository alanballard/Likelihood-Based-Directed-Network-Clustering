#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>	
#include <time.h>		
#include <string>		
#include <iomanip>
#include <math.h>       
#include <algorithm>   
#include <vector>
#include <cctype>		
#include <tuple>

using namespace std;

//FUNCTION PROTOTYPES
int initial_pop(vector<int>& group_size, vector<int>& cluster_membership, int k, int N); //create initial membership vector
vector<int> mutate_membership(vector<int> current_membership, vector<int> group_size); //mutate membership vector
long double uniform_rand(int seed); //random (0,1) draw

int likelihood(int network_key, int N, vector< vector< tuple<int, int, char> >> links, int nb_links, double InitTemp, double CR, int TL, int Max_Success, int min_k, int max_k, int k_int)
/*
				//Populate observed OBS edge weight matrix. This requires 1 full traversal of adjacency list
				for (int i = 0; i < links.size(); i++)
				{
					for (int j = 0; j < links[i].size(); j++)
					{
						if (get<2>(links[i][j]) == 'T')//edge is FROM src vert links[i] TO dest vert links[i][j]...
						{
							if (cluster_membership[i] == cluster_membership[get<0>(links[i][j])])//src/dest verts belong to same cluster
							{
								//OBS[cluster_membership[i] - 1][cluster_membership[i] - 1]
								//	= OBS[cluster_membership[i] - 1][cluster_membership[i] - 1]
								//	+ get<1>(links[i][j]);
								OBS[cluster_membership[i] - 1]=OBS[cluster_membership[i] - 1] + get<1>(links[i][j]);
							}
							else//src/dest clusters differ
							{
								//OBS[cluster_membership[i] - 1][cluster_membership[get<0>(links[i][j])] - 1]
								//	= OBS[cluster_membership[i] - 1][cluster_membership[get<0>(links[i][j])] - 1]
								//	+ get<1>(links[i][j]);
								OBS[k] = OBS[k] + get<1>(links[i][j]);
							}
						}
						//Note that ever edge is represented twice in the adjacency lists used for this program.
						//One entry is TO ('T') vertex j from vertex i, and one entry is FROM ('F') vertex i to vertex j.
						//It would double count each edge if we considered both of these entries, so we're only counting the ones
						//marked with 'T' above.
					}
				}

//May need to cut this in 1/2, like undirected case:
				//Note: Since the adjacency list contains 2 entries for every edge and the number of observed edges are counted
				//over the adjacency list, the total observed count needs to be divided by 2.
				//for (int i = 0; i < k + 1; i++) { OBS[i] = OBS[i] / 2; }

				//Populate possible edge count vector. This is done using the cluster size information
				for (int i = 0; i < k; i++)
				{
					//POS[i][i] = ((unsigned long long int) group_size[i] * (group_size[i] - 1)); //2*[nchoose2]=2*[N(N-1)/2]=N*N-1
					POS[i] = ((unsigned long long int) group_size[i] * (group_size[i] - 1)); //2*[nchoose2]=2*[N(N-1)/2]=N*N-1
				}
				for (int i = 0; i < k - 1; i++)
				{
					for (int j = i + 1; j < k; j++)
					{
						//add number of edges possible between cluster i and cluster j. This is only run once, on the first run.
						//Since each i/j combination is only considered once, can't we just shorten this to?:
						//POS[i][j] = (unsigned long long int) group_size[i] * group_size[j]; //But keeping it as-is shouldn't hurt anything either...
						//POS[i][j] = POS[i][j] + (unsigned long long int) group_size[i] * group_size[j];
						//POS[j][i] = POS[i][j];
						POS[k] = POS[k] + (unsigned long long int) group_size[i] * group_size[j];
					}
				}
*/
{
	int final_k = 0;
	long double final_BIC = 0;
	long double final_LOGLIK = 0;
	long double final_MODULARITY = 0;
	long double final_TESTSTAT = 0;
	vector<int> final_membership(N, 0);
	long double final_TIME = 0;

	//Repeat Simulated Annealing process for each k in the range of specified cluster numbers
	for (int k = min_k; k <= max_k; k = k + k_int)
	{
cout<<"k="<<k<<endl;
		clock_t start_new = clock(); //Time until solution

//Step: Initialize k+1 current/proposed OBS and POS vectors, and current/proposed membership vector
		vector<unsigned long long int> OBS(k + 1, 0);//elements 0 thru k-1 are within cluster values, element k is between-cluster value
		vector<unsigned long long int> POS(k + 1, 0);//elements 0 thru k-1 are within cluster values, element k is between-cluster value
		//Data types available: http://www.cplusplus.com/doc/tutorial/variables/

		//Initialize and populate cluster assignment and size vectors
		vector<int> group_size(k, 0);//used in initial pop to store number of nodes in each group 1:k
		vector<int> cluster_membership;	//holds cluster assignments of N vertices

		//vertex to be reassigned and the cluster to which it will be reassigned
		vector<int> reassignment;
		int l_star;		//The vertex to be reassigned
		int clust_l;	//L-star's cluster membership before reassignment
		int clust_j;	//L-star's cluster membership after reassignment

		//Total weights and group counts used to develop parameters for proposed clustering
//In paper proofs,  Z_(to/from *) is relative to vertices of cluster *.
//Total weight of links from vertex l* to vertices of cluster * (These are the Z_(to *) variables in the paper). Z_to[0], for example, is the total weight of links from vertex l* to vertices of cluster 1.
		vector<int> Z_to(k);
//Total weight of links from vertices of cluster * to vertex l* (These are the Z_(from *) variables in the paper). Z_from[0], for example, is the total weight of links from vertices of cluster 1 to vertex l*.
		vector<int> Z_from(k);

		//Total weights and group counts used to develop parameters for proposed clustering
		int Z_to_l;		//Sum of edge weights from vertex l* to cluster l
		int Z_from_l;	//Sum of edge weights from cluster l to vertex l*
		int Z_to_j;		//Sum of edge weights from vertex l* to cluster j
		int Z_from_j;	//Sum of edge weights from cluster j to vertex l*

		//Observed edge weights and possible edge counts used to develop parameters for proposed clustering.
		unsigned long long int OBS_l_1;
		unsigned long long int OBS_j_1;
		unsigned long long int OBS_bw_1;

		unsigned long long int POS_l_1;
		unsigned long long int POS_j_1;
		unsigned long long int POS_bw_1;

		//Parameters for current clustering	
		long double param_l_0;
		long double param_j_0;
		long double param_bw_0;

		//Parameters for proposed clustering
		long double param_l_1;
		long double param_j_1;
		long double param_bw_1;

		//Partials loglikelihoods for current and proposed clustering, and the change in loglikelihood.
		long double l_0;
		long double l_1;
		long double delta_l;

		//Densities in the denominator in the change-in-loglik formula [for current clustering]
		long double l_density_0;
		long double j_density_0;
		long double bw_density_0;

		//Densities in the numerator in the change-in-loglik formula [for proposed clustering]
		long double l_density_1;
		long double j_density_1;
		long double bw_density_1;

		//Modularity of solution membership assignment
		long double current_full_modularity = 0; 

		//Flag to allow creation of initial clustering that will be mutated on later runs.
		int first_run = 1;

		double LimitIT = 1.0e-8; //Minimum temperature. Program stops when this number is reached
		long double epsilon = 1.0e-6; //Minimum change in loglikelihood required to accept proposed membership mutation
		int TL_counter = 0; //Reset count length at the start of each k-loop
		int Success_counter = 0; //Reset success counter at the start of each k-loop
		double IT = InitTemp;  		//Reset initial temperature at the start of each k-loop

		//START Simulated Annealing Algorithm
		while (IT > LimitIT)
		{

			TL_counter++; //Iterate number of attempts at current temperature

			//After TL attempts at reclustering at the current temperature, or Max_Success # of successes at curren temperature
			//Decrease temperature
			if (TL_counter >= TL || Success_counter >= Max_Success)
			{
				IT = CR * IT;
				TL_counter = 0;
				Success_counter = 0;
			}

			//Generate initial clustering (1st run only). This is the initial current clustering.
			//On runs >1, the "current clustering" will either be the one generated here, or a proposed clustering that was accepted below.
			while (first_run == 1)
			{
				//Step: Make initial membership
				//Populate cluster assignment and size vectors
				initial_pop(group_size, cluster_membership, k, N); //assigns initial cluster assignments to z, and counts cluster sizes
				first_run = 0;

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
							+ 2*((unsigned long long int)group_size[i] * group_size[j])							
							;
					}
				}

			}//end 1st-run WHILE loop

			//Enter SA Algorithm
//Step: Propose a new clustering by mutating current clustering
			//Select Random Vertex l* for reassignment and select cluster to which it will be moved
			//This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
			reassignment = mutate_membership(cluster_membership, group_size);	//Proposed membership vector due to reassignment
			l_star = reassignment[0];				//The vertex to be reassigned
			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
			clust_j = reassignment[1];				//l-star's cluster membership after reassignment

//Step: Calculate abbreviated loglikelihood of the current membership (Binomial)
			//Get denominator of change-in-loglik equation [l,j,bw for current clustering]
			//Calculate parameters for current clustering.
			param_l_0 = OBS[clust_l - 1] / (long double)POS[clust_l - 1];
			param_j_0 = OBS[clust_j - 1] / (long double)POS[clust_j - 1];
			param_bw_0 = OBS[k] / (long double)POS[k];

			l_density_0 = 0;
			j_density_0 = 0;
			bw_density_0 = 0;
			//Calculate l,j,bw densities for current clustering, the denominator in the change-in-loglik formula
			//if unweighted, all edges have weight=1 and Binomial density is appropriate

				//Note: we are calculating the log of the product of 3 binomial densities each for
				//the numerator and denominator of the change-in-loglik formula.
				//If a parameter estimate is 0 or 1, this results in a density estimate of 0 or INF, which throws an error when the log is taken.
				//Thus, we rewrite the log-product as a sum of logs. If a density has a parameter estimate of 0 or 1, 
				//we force that density's log-value to be zero so that it does not throw an error and allows a modified log-likelihood calculation.

			if (param_l_0 == 0 || param_l_0 == 1) { l_density_0 = 0; }
			else {
				l_density_0 = (OBS[clust_l - 1])*log(param_l_0) + (POS[clust_l - 1] - OBS[clust_l - 1])*log(1 - param_l_0);
			}
			if (param_j_0 == 0 || param_j_0 == 1) { j_density_0 = 0; }
			else {
				j_density_0 = (OBS[clust_j - 1])*log(param_j_0) + (POS[clust_j - 1] - OBS[clust_j - 1])*log(1 - param_j_0);
			}
			if (param_bw_0 == 0 || param_bw_0 == 1) { bw_density_0 = 0; }
			else {
				bw_density_0 = (OBS[k])*log(param_bw_0) + (POS[k] - OBS[k])*log(1 - param_bw_0);
			}
			l_0 = l_density_0 + j_density_0 + bw_density_0;

//Step: Calculate abbreviated loglikelihood of the proposed membership (Binomial)	
			//Get numerator of change-in-loglik equation [l,j,bw for proposed clustering]
			 //Zero out Z (to/from cluster size) vectors. Best practice for runs > 1			
			for (int i = 0; i < Z_to.size(); i++) { Z_to[i] = 0; } //Total weight of links to l*, by cluster (These are the Z_(to s) values in the paper)			
			for (int i = 0; i < Z_from.size(); i++) { Z_from[i] = 0; } //Total weight of links from l*, by cluster (These are the Z_(from s) values in the paper)

			//Populate to/from cluster size vectors
			//Looking at l*'s neighbors in the adjacency list
			for (int j = 0; j < links[l_star].size(); j++)
			{
				//get<0>(links) -Destination vertex number
				//get<1>(links) -Weight
				//get<2>(links) -char in {T,F}
				if (get<2>(links[l_star][j]) == 'T')//edge extends FROM vertex l* TO the jth vertex in l*'s adjacency list row
				{
					Z_to[cluster_membership[get<0>(links[l_star][j])] - 1] = Z_to[cluster_membership[get<0>(links[l_star][j])] - 1] + get<1>(links[l_star][j]);
				}
				else//edge extends FROM jth vertex in l*'s adjacency list row TO vertex l*
				{
					Z_from[cluster_membership[get<0>(links[l_star][j])] - 1] = Z_from[cluster_membership[get<0>(links[l_star][j])] - 1] + get<1>(links[l_star][j]);
				}
			}
/*
				std::cout << "Cluster memberships:" << endl;
				for (int i = 0; i < cluster_membership.size(); i++){ std::cout << cluster_membership[i] << "   "; } std::cout << endl;

				for (int i = 1; i <= k; i++){
				std::cout << endl << "The number of edges from vertex " << l_star << " to cluster " << i << " is " << Z_to[i - 1] << endl;}
				for (int i = 1; i <= k; i++){
				std::cout << endl << "The number of edges from cluster " << i << " to vertex " << l_star << " is " << Z_from[i - 1] << endl;}

				std::cout << endl << "Z_to():" << endl;
				for (int i = 0; i < Z_to.size(); i++){std::cout << "<" << i + 1 << ", " << Z_to[i] << "> ";}

				std::cout << endl << "Z_from():" << endl;
				for (int i = 0; i < Z_from.size(); i++){std::cout << "<" << i + 1 << ", " << Z_from[i] << "> ";}
*/		

			//Calculate parameters for proposed clustering
			//Z_to/from_l/j from paper
			Z_to_l		= Z_to[clust_l - 1]; //Sum of edge weights from vertex l* to cluster l
			Z_from_l	= Z_from[clust_l - 1]; //Sum of edge weights from cluster l to vertex l*
			Z_to_j		= Z_to[clust_j - 1]; //Sum of edge weights from vertex l* to cluster j
			Z_from_j	= Z_from[clust_j - 1]; //Sum of edge weights from cluster j to vertex l*

			//Equations for updating OBS and POS given a reclustering, taken from paper
			OBS_l_1		= OBS[clust_l - 1] - Z_to_l - Z_from_l;
			OBS_j_1		= OBS[clust_j - 1] + Z_to_j + Z_from_j;
			OBS_bw_1	= OBS[k] - Z_to_j  - Z_from_j + Z_to_l + Z_from_l;

			POS_l_1 = ((group_size[clust_l - 1] - 1)*(group_size[clust_l - 1] - 2)) ;
			POS_j_1 = ((group_size[clust_j - 1] + 1)*group_size[clust_j - 1]) ;
			POS_bw_1 = POS[k] + 2*((group_size[clust_l - 1] ) - (group_size[clust_j - 1]) - 1);

			param_l_1 = OBS_l_1 / (long double)POS_l_1;
			param_j_1 = OBS_j_1 / (long double)POS_j_1;
			param_bw_1 = OBS_bw_1 / (long double)POS_bw_1;

			l_density_1 = 0;
			j_density_1 = 0;
			bw_density_1 = 0;

			//Calculate l,j,bw densities for PROPOSED clustering
			//if unweighted, all edges have weight=1 and Binomial density is appropriate
				//Calculate modified new log-lik, the numerator in the change-in-loglik formula
			if (param_l_1 == 0 || param_l_1 == 1) { l_density_1 = 0; }
			else {
				l_density_1 = (OBS_l_1)*log(param_l_1) + (POS_l_1 - OBS_l_1)*log(1 - param_l_1);
			}
			if (param_j_1 == 0 || param_j_1 == 1) { j_density_1 = 0; }
			else {
				j_density_1 = (OBS_j_1)*log(param_j_1) + (POS_j_1 - OBS_j_1)*log(1 - param_j_1);
			}
			if (param_bw_1 == 0 || param_bw_1 == 1) { bw_density_1 = 0; }
			else {
				bw_density_1 = (OBS_bw_1)*log(param_bw_1) + (POS_bw_1 - OBS_bw_1)*log(1 - param_bw_1);
			}
			l_1 = l_density_1 + j_density_1 + bw_density_1;
			
/*for (int i = 0; i < k + 1; i++) { cout << OBS[i] << "..."; }cout << endl;
for (int i = 0; i < k + 1; i++) { cout << POS[i] << "..."; }cout << endl;
cout << param_l_0 << "***" << param_j_0 << "***" << param_bw_0 << endl;
for (int i = 0; i < k; i++) { cout<< group_size[i] << "^^^"; }cout << endl;
cout << "from clust " << clust_l << " to clust " << clust_j << endl;
cout << "*********************" << endl;
{
	int i = l_star;
	cout << "Vertex " << i << ":";
	for (unsigned int j = 0; j < links[i].size(); j++) {
		cout << " " << get<0>(links[i][j]) << " <" << get<1>(links[i][j]) << ", " << get<2>(links[i][j]) << ">";
		//cout << " " <<  links[i][j].first << " <"<< links[i][j].second <<">";
	}
	cout << endl;
}
cout << "*********************" << endl;
cout << OBS_l_1 << "..." << OBS_j_1 << "..." << OBS_bw_1 << endl;
cout << POS_l_1 << "..." << POS_j_1 << "..." << POS_bw_1 << endl;
cout << param_l_1 << "***" << param_l_1 << "***" << param_bw_1 << endl;
cout << "##################################################" << endl;
cin.get();*/
//for (int i = 0; i < k + 1; i++) { cout << OBS[i]/ (long double)POS[i] << "..."; }cout << endl;


//Step: Calculate change in log-likelihood (Binomial) due to reclustering
			//log(new/old)=log(new)-log(old)

			//If all densities in either the new or the old likelihood were discarded (and thus set to zero),
			//it is not possible to evaluate a change-in-likelihood
			//set entire change-in-likelihood to zero so a new clustering will be considered
			//if (l_0 == 0 || l_1 == 0) { delta_l = 0; }
			//else { delta_l = l_1 - l_0; }
			delta_l = l_1 - l_0;
			//if delta_l>0, then l_1>l_0, and the reclustered membership is more "likely"
			//Otherwise, l_1<l_0 and the current membership is more "likely"

			int STAGE = 0; //Will indicate when/if a proposed recustering was accepted or rejected.

			if (delta_l > epsilon)
			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
//Step: Compare the two likelihoods. If accept new likelihood, then:
				//Update current membership vector to reflect l_star's new cluster membership
				//Update 3 elements (each) of OBS/POSvectors 
				//Update group (cluster) size vector
				//Proposed loglikelihood becomes current loglikelihood.	
				STAGE = 1;

				//current clustering <- proposed clustering
				cluster_membership[l_star] = clust_j;

				//Update OBS vector
				OBS[clust_l - 1] = OBS_l_1;
				OBS[clust_j - 1] = OBS_j_1;
				OBS[k] = OBS_bw_1;

				//Update POS vector
				POS[clust_l - 1] = POS_l_1;
				POS[clust_j - 1] = POS_j_1;
				POS[k] = POS_bw_1;

				//Update group size vector
				group_size[clust_l - 1]--; //cluster l loses 1 vertex
				group_size[clust_j - 1]++; //cluster j gains 1 vertex

				l_0 = l_1;
				l_1 = 0;

				Success_counter++; //Iterate number of success at this temperature
			}
			else
			{
				long double one = 1;      //used in MIN call to match data types.

				double uni_draw = rand() / double(RAND_MAX);
				double min_val = min(one, exp(delta_l / IT));

				if (uni_draw < min_val) //Accept proposed clustering
				{
//Step: If accept new loglikelihood, then:
					//Update current membership vector to reflect l_star's new cluster membership
					//Update 3 elements (each) of OBS/POSvectors 
					//Update group (cluster) size vector
					//Proposed loglikelihood becomes current loglikelihood.	
					STAGE = 2;
					//current clustering <- proposed clustering
					cluster_membership[l_star] = clust_j;

					//Update OBS vector
					OBS[clust_l - 1] = OBS_l_1;
					OBS[clust_j - 1] = OBS_j_1;
					OBS[k] = OBS_bw_1;

					//Update POS vector
					POS[clust_l - 1] = POS_l_1;
					POS[clust_j - 1] = POS_j_1;
					POS[k] = POS_bw_1;

					//Update group size vector
					group_size[clust_l - 1]--; //cluster l loses 1 vertex
					group_size[clust_j - 1]++; //cluster j gains 1 vertex

					l_0 = l_1;
					l_1 = 0;

					//Iterate number of success at this temperature
					Success_counter++;
				}
				else //Reject proposed clustering
				{//Nothing happens here. The current cluster is retained and the process repeats itself.					
					STAGE = 3;
				}
			}
			//At this point, the proposed clustering has:
			//1. Been Accepted. membership vector/OBS/POS/cluster counts have beeen updated to reflect new clustering
			//2. Rejected. memberhsip vector/OBS/POS/cluster counts still correspond to the current clustering.


/*			double l_0_w = 0; //loglikelihood of within-clusters
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
			cout << "likelihood is: " << l_0_w + l_0_b << endl;
*/

		}//End SA LOOP

//CALCULATE TIME REQUIRED TO REACH SOLUTION
		//cout << fixed;
		cout << setprecision(10);
		clock_t end_new = clock();
		double time = (double)(end_new - start_new) / CLOCKS_PER_SEC;

		//At this point, a solution clustering has been identified.
		long double full_loglik_curr = 0; //loglikelihood for the full clustering
		long double LO = 0;				//loglikelihood under null hypothesis: # of clusters=1
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
		full_loglik_curr = l_0_w + l_0_b;

//STEP: Calculate likelihood under null hypothesis of k=1 cluster
		long long int LO_OBS = 0;
		long long int LO_POS = 0;
		for (int i = 0; i <= k; i++)
		{
			LO_OBS = LO_OBS + OBS[i];
		}
		LO_POS = N * (N - 1) ; //2*(N choose 2)
		LO = (LO_OBS)*log(LO_OBS / (long double)LO_POS) + (LO_POS - LO_OBS)*log(1 - (LO_OBS / (long double)LO_POS));

//Step: Calculate DIRECTED modularity of solution membership 
				//Have to traverse entire adjacency list...
				//Formula taken from "Community Structure in Directed Networks", Leicht and Newman, 2007
				//(1/2m)*sum[A_ij - (k_i^in*k_j^out)/2m]
				//A_ij		=1 if edge exists FROM j TO i	= where get<2>(links[i][j]) == 'F'
				//k_i^in	=in degree of vertex i			=k_i_in[i]
				//k_j^out	=out degree of vertex j			=k_j_out[j]
				//m			=# of edges in network			=E
		//Get number of distinct edges in the network
		int E = nb_links;

		vector<long double> k_i_in(N, 0); //vertex i's in-degree
		vector<long double> k_j_out(N, 0); //vertex j's out-degree
		//Note: It's a little confusing but the k_i_in and k_j_out notation are what's used in "Community Structure in Directed Networks"
		for (int s = 0; s < links.size(); s++)
		{
			for (int t = 0; t < links[s].size(); t++)
			{
				if (get<2>(links[s][t]) == 'F') { k_i_in[s]++; }//if edge extends from t to s, increase s's in degree by 1
				else { k_j_out[s]++; }//else, edge extends from s to t, increase s's out degree by one
			}
		}

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

//STEP: Calculate Bayesian Information Criterion for solution clustering
		//In old code, BIC = -2 * (LO + 0.5*(-bestfit)) + (k + 1)*log(N*(N - 1) / 2);
		//but bestfit = -(-2 * (LO - full_loglik_curr)), so this can be simplified:		
		//number of parameters=(k+1), number of observations=Nchoose2=N(N-1)/2, maximized loglikelihood=full_loglik_curr		
		BIC = -2 * (full_loglik_curr)+(k + 1)*log(N*(N - 1) / (long double)2);
		TESTSTAT = -2 * (LO - full_loglik_curr);

		//Identify the optimum number of clusters, based on minimum Bayesian Information Criterion
		if (k == min_k) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = full_loglik_curr;
			final_TESTSTAT = TESTSTAT;
			final_membership = cluster_membership;
			final_MODULARITY = current_full_modularity;
			final_TIME = time;
		}
		//if BIC < min BIC thus far, update variables
		if (BIC < final_BIC) {
			final_k = k;
			final_BIC = BIC;
			final_LOGLIK = full_loglik_curr;
			final_TESTSTAT = TESTSTAT;
			final_membership = cluster_membership;
			final_MODULARITY = current_full_modularity;
			final_TIME = time;
		}

//STEP: Output solution clustering statistics at each value of k		
		ostringstream atkstring;
		ofstream atkinfo;
		// If old file exists from a previous run, delete it.
		//remove("completion_likelihood_all_sols.txt");
		atkstring << "completion_likelihood_all_sols" << ".txt";
		atkinfo.open("completion_likelihood_all_sols.txt", ios::app);
		atkinfo << network_key << "," << N << "," << nb_links << "," << k << "," << current_full_modularity << "," << full_loglik_curr << "," << BIC << "," << TESTSTAT << "," << time << endl;
		atkinfo.close();


		//cout << "**********************************" << endl;
	}//end for-k loop

//STEP: Output best solution across all k
	ostringstream atkstring2;
	ofstream atkinfo2;
	// If old file exists from a previous run, delete it.
	//remove("completion_likelihood_final_sol.txt");
	atkstring2 << "completion_likelihood_final_sol" << ".txt";
	atkinfo2.open("completion_likelihood_final_sol.txt", ios::app);
	atkinfo2 << network_key << "," << N << "," << nb_links << "," << final_k << "," << final_MODULARITY << "," << final_LOGLIK << "," << final_BIC << "," << final_TESTSTAT << "," << final_TIME << endl;
	atkinfo2.close();

	//STEP: Output best solution clustering
	ostringstream bestzsol;
	ofstream bestzstorage;
	// If old file exists from a previous run, delete it.
	//remove("solution_likelihood.txt");
	bestzsol << "solution_likelihood" << ".txt";
	bestzstorage.open("solution_likelihood.txt", ios::app);
	for (int i = 0; i < N; i++) {
		bestzstorage << final_membership[i];
		if (i < N - 1) { bestzstorage << ", "; }
	}
	bestzstorage << endl;
	bestzstorage.close();


	return 0;
}

