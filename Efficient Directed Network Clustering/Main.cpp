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
vector< vector< tuple<int, int, char>    >> read_el(char *filename, char weighted); //Read in edge list and convert to adjacency list
void print_el(vector<vector< tuple<int, int, char> >> links); //Print adjacency and/or edge list
int initial_pop(vector<int>& group_size, vector<int>& cluster_membership, int k, int N); //create initial membership vector
vector<int> mutate_membership(vector<int> current_membership, vector<int> group_size); //mutate membership vector
long double uniform_rand(int seed); //random (0,1) draw

int likelihood(int network_key, int N, vector< vector< tuple<int, int, char> >> links, int nb_links, double InitTemp, double CR, int TL, int Max_Success, int min_k, int max_k, int k_int);
int modularity(int network_key, int N, int nb_links, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success, int min_k, int max_k, int k_int);

//Code for generating directed LFR benchmark networks
int main_LFR(
	int num_nodes,				//N		number of nodes
	double average_k,			//k		average in-degree
	int max_degree,				//maxk	maximum in-degree
	double tau,					//t1		minus exponent for the degree sequence
	double tau2,				//t2		minus exponent for the community size distribution
	double mixing_parameter,	//mu		mixing parameter
	int nmin,					//minc	minimum for the community sizes
	int nmax					//maxc	maximum for the community sizes
);

int main()
{
	srand(time(0));          //seeds random number generator

	//Simulated Annealing Algorithm Parameters
	double InitTemp = 1;	//Starting temperature 1
	double CR = 0.99;		//Rate at which temperature is reduced .99

	int TL = 100;			//Maximum number of reclustering attempts at a given temperature 100
	int Max_Success = TL;	//Maximum number of successes allowed at a given temperature 
	//The parameters below are hard-coded in the model programs
	//double LimitIT = ; //Minimum temperature. Program stops when this number is reached
	//long double epsilon =; //Minimum change in loglikelihood required to accept proposed membership mutation

	//Get choice of user input file or LFR file.
	int filesource=0;
	do
	{
		std::cout << "Select network file source:" << endl;
		std::cout << "1=user file" << endl;
		std::cout << "2=generate LFR networks" << endl;
		std::cin >> filesource;
	} while (
		!std::cin.fail()
		&& filesource != 1
		&& filesource != 2
		);
	std::cout << endl;

	char filename[256];
	char weighted = 'N';

	vector< vector< tuple<int, int, char> >> links;
	int N = 0;

	if (filesource == 1) 
	{//if user input file, get name of file containing network and proceed as usual
						  //Get file name
		
						std::cout << "Enter file name, including suffix (ex. myfile.txt): ";
						std::cin.getline(filename, 256);
						//Make sure file exists. otherwise, get file name again.
						int file_exists = 0;
						ifstream ifile(filename);
						file_exists = (bool)ifile;

						while (file_exists == 0)
						{
							std::cout << "Enter file name, including suffix (ex. myfile.txt): ";
							std::cin.getline(filename, 256);
							ifstream ifile(filename);
							file_exists = (bool)ifile;
						}
						std::cout << endl;

						//Input network edge list and convert to adjacency list that includes weight information
						std::cout << "Reading in edge list..." << endl;
						time_t read_start = time(0);//Time until solution
						links = read_el(filename, weighted);

						//Check to see that every vertex in the network has at least 1 edge to another vertex. No singletons allowed.
						int missing = 0;
						for (int i = 0; i < links.size(); i++)
						{
							if (links[i].size() == 0) {
								missing = 1;
								cout << "empty row at vertex " << i << endl;
							}
						}
						if (missing == 1) {
							std::cout << "WARNING WARNING WARNING" << endl;
							std::cout << "*** ADJACENCY LIST IS NOT FULLY POPULATED ***" << endl;
							std::cout << "THE PROGRAM WILL NOT WORK CORRECTLY" << endl;
							std::cin.get();
						}
//print_el(links);

						//GET STATS ON CURRENT NETWORK
						int N = 0;
						int E = 0;
						//Get the number of distinct vertices in the network
						N = links.size();
						//Get number of distinct edges in the network
						for (int i = 0; i < links.size(); i++)
						{
							for (int j = 0; j < links[i].size(); j++) {
								if (get<2>(links[i][j]) == 'T') {
									E++;
								}
							}
						}

						std::cout << "Vertex count: " << links.size() << endl;
						time_t read_end = time(0);
						double read_time = difftime(read_end, read_start);
						std::cout << read_time << " seconds required to read edge list" << endl; //divide by 60 to get minutes

						int min_k = 0; int max_k = 0; int k_int = 0;
						std::cout << "Enter the minimum number of clusters in which to cluster your network (should be >=2):" << endl;
						std::cin >> min_k;
						cout << endl;
						std::cout << "Enter the maximum number of clusters in which to cluster your network (should be <<<N):" << endl;
						std::cin >> max_k;
						std::cout << "Enter cluster count step size (ex. 1 means evaluate every k between the min and max, 2 means evaluate every other k, etc.):" << endl;
						std::cin >> k_int;

						int overall_net_ct = 1;
						//Run all models on current network
						cout << "\t" << "Evaluating Model: likelihood" << endl;
						likelihood(overall_net_ct, N, links, E, InitTemp, CR, TL, Max_Success, min_k, max_k, k_int);
						cout << "\t" << "Evaluating Model: modularity" << endl;
						modularity(overall_net_ct, N, E, links, InitTemp, CR, TL, Max_Success, min_k, max_k, k_int);

	}

	else 
	{//Produce full set of LFR benchmark networks and evaluate both modularity and likelihood on each network

						/*
										int num_nodes;//N		number of nodes
										double average_k;//k		average in-degree
										int max_degree;//maxk	maximum in-degree
										double tau;//t1		minus exponent for the degree sequence
										double tau2;//t2		minus exponent for the community size distribution
										double mixing_parameter;//mu		mixing parameter
										int nmin;//minc	minimum for the community sizes
										int nmax;//maxc	maximum for the community sizes
		
										std::cout << endl;
										std::cout << "Example LFR values: 50, 5, 10, 3, 1, 0.2, 5, 10" << endl;

										std::cout << "Enter the desired number of nodes in the LFR Graph" << endl;
										std::cin >> num_nodes;
										cout<<endl;
										std::cout << "Enter the desired average in-degree in the LFR Graph" << endl;
										std::cin >> average_k;
										cout << endl;
										std::cout << "Enter the desired maximum in-degree in the LFR Graph" << endl;
										std::cin >> max_degree;
										cout << endl;
										std::cout << "Enter the desired exponent for the degree sequence in the LFR Graph" << endl;
										std::cin >> tau;
										cout << endl;
										std::cout << "Enter the desired exponent for the community size distribution in the LFR Graph" << endl;
										std::cin >> tau2;
										cout << endl;
										std::cout << "Enter the desired mixing parameter in the LFR Graph" << endl;
										std::cin >> mixing_parameter;
										cout << endl;
										std::cout << "Enter the desired minimum community size in the LFR Graph" << endl;
										std::cin >> nmin;
										cout << endl;
										std::cout << "Enter the desired maximum community size in the LFR Graph" << endl;
										std::cin >> nmax;

		
										main_LFR(
											num_nodes,
											average_k,
											max_degree,
											tau,
											tau2,
											mixing_parameter,
											nmin,
											nmax
										);

										strcpy_s(filename, "LFR_network.dat");
										ifstream ifile(filename);

										//Input network edge list and convert to adjacency list that includes weight information
										std::cout << "Reading in LFR edge list..." << endl;
										time_t read_start = time(0);//Time until solution
										links = read_el(filename, weighted);

										//Check to see that every vertex in the network has at least 1 edge to another vertex. No singletons allowed.
										int missing = 0;
										for (int i = 0; i < links.size(); i++)
										{
											if (links[i].size() == 0) { missing = 1; }
										}
										if (missing == 1) {
											std::cout << "WARNING WARNING WARNING" << endl;
											std::cout << "*** ADJACENCY LIST IS NOT FULLY POPULATED ***" << endl;
											std::cout << "THE PROGRAM WILL NOT WORK CORRECTLY" << endl;
											std::cin.get();
										}

										//Get the number of distinct vertices in the network
										N = links.size();

										std::cout << "Vertex count: " << links.size() << endl;
										time_t read_end = time(0);
										double read_time = difftime(read_end, read_start);
										std::cout << read_time << " seconds required to read edge list" << endl; //divide by 60 to get minutes
							*/
	
					//PRINT NETWORK as adjacency list, including weights and direction
					//print_el(links); 
					//std::cin.get(); std::cin.get();

				/**********************************/
				int overall_net_ct = 0;
				double overall_time = 0;

				//Generate network at current (fixed) settings
					int num_nodes = 100; //# of nodes in network
					int max_degree = 10;//20 //maximum in-degree
					int nmin = 2;	//minimum for the community sizes
					int nmax = 20;	//maximum for the community sizes
				//NOTE!!!!: average_k, tau, tau2, mixing_paramter and network_ct parameters are set in the FOR loops below

				//int min_cluster_size = num_nodes;
				//int max_cluster_size = 0;

				for(double average_k=5; average_k<=10; average_k=average_k+2.5)//average in-degree //average_k<=10
				{
					for(double tau=2; tau<=3; tau=tau+1)//minus exponent for the degree sequence //tau<=3
					{
						for(double tau2=1; tau2<=2; tau2=tau2+1)//minus exponent for the community size distribution //tau2<=2
						{
							for(double mixing_parameter=0.1; mixing_parameter<=0.95; mixing_parameter=mixing_parameter+0.05)//{0.10,0.15,0.20,0.25...,0.55,0.60} //mixing parameter mixing_parameter<=0.95
							{//was 0.1 to 0.6, changed to 0.65 to 0.95
								for (int network_ct = 1; network_ct <= 10; network_ct++) //Number of networks to create at current parameter settings //network_ct <= 10
								{	

				overall_net_ct = overall_net_ct + 1;
				time_t network_start = time(0);

				cout<<"Generating current network, #"<< overall_net_ct <<endl;
				int number_of_LFR_clusters = 0;
				number_of_LFR_clusters=main_LFR(
										num_nodes,
										average_k,
										max_degree,
										tau,
										tau2,
										mixing_parameter,
										 //These values below are estimated from fig. 10, pg 183, of original pager for N=100 instead of N=500
										nmin,
										nmax
										);
				//if (number_of_LFR_clusters > max_cluster_size) { max_cluster_size = number_of_LFR_clusters; }
				//if (number_of_LFR_clusters < min_cluster_size) { min_cluster_size = number_of_LFR_clusters; }

				////Get distribution of [# of clusters] across networks generated at current pramaters
				//						ostringstream atkstring;
				//						ofstream atkinfo;
				//						// If old file exists from a previous run, delete it.
				//						//remove("nbr_clusters_distribution.txt");
				//						atkstring << "nbr_clusters_distribution" << ".txt";
				//						atkinfo.open("nbr_clusters_distribution.txt", ios::app);
				//						atkinfo << number_of_LFR_clusters << endl;
				//						atkinfo.close();


				cout << "Currently evaluating network with parameters: " << num_nodes << "," << average_k << "," << max_degree << "," << tau << "," << tau2 << "," << mixing_parameter << "," << nmin << "," << nmax << endl;
									ostringstream PARAMETERS_NAME;
									ofstream NETWORK_PARAMETERS;
									PARAMETERS_NAME << "LFR_PARAMETERS" << ".txt";
									NETWORK_PARAMETERS.open("LFR_PARAMETERS.txt", ios::app);
									NETWORK_PARAMETERS << num_nodes << "," << average_k << "," << max_degree << "," << tau << "," << tau2 << "," << mixing_parameter << "," << nmin << "," << nmax;
									NETWORK_PARAMETERS << endl;
									NETWORK_PARAMETERS.close();

									//READ IN CURRENT NETWORK
									char filename[256];		
									strcpy_s(filename, "LFR_network.dat");
									ifstream ifile(filename);

									vector< vector< tuple<int, int, char> >> links;
									char weighted = 'N';
									links = read_el(filename, weighted);

									//GET STATS ON CURRENT NETWORK
									int N = 0;
									int E = 0;
									//Get the number of distinct vertices in the network
									N = links.size();
									//Get number of distinct edges in the network
									for (int i = 0; i < links.size(); i++)
									{
										for (int j = 0; j < links[i].size(); j++) {
											if (get<2>(links[i][j]) == 'T') {
												E++;
											}
										}
									}

									//Range of clusters over which we evaluate each model
									int min_k = 10;
									int max_k = 20;//20
									int k_int = 1;

									//Run all models on current network
				cout << "\t" << "Evaluating Model: likelihood" << endl;
									likelihood(overall_net_ct, N, links, E, InitTemp, CR, TL, Max_Success, min_k, max_k, k_int); 
				cout << "\t" << "Evaluating Model: modularity" <<endl;
									modularity(overall_net_ct, N, E, links, InitTemp, CR,TL, Max_Success, min_k, max_k, k_int);

					time_t network_end = time(0);
					double network_time = difftime(network_end, network_start);
					overall_time = overall_time + network_time;
					std::cout << "Time for current network = " << network_time << " seconds, " << network_time/60 <<" minutes"<<endl; //time / 60=minutes
					std::cout << "Current average time = " << overall_time/ overall_net_ct << " seconds, " << (overall_time / overall_net_ct) / 60 << " minutes" << endl; //time / 60=minutes
								}
							}
						}		
					}
				}

				/**********************************/
				//cout << "min_cluster_size=" << min_cluster_size << endl;
				//cout << "max_cluster_size=" << max_cluster_size << endl;

	}

	std::cout << "Program ended. Hit enter to exit.";
	std::cin.get(); std::cin.get();

}//END main()


//FUNCTION: Returns a U(0,1) variable
long double uniform_rand(int seed)
{

	return seed * (1.0 / (RAND_MAX + 1.0));     //to call the random U(0,1), double r = uniform_rand ( rand() );
}//end uniform_rand function definition
