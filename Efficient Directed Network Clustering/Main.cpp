
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
vector< vector< tuple<int, int, char>    >> read_el(char *filename, char weighted); //Read in edge list and convert to adjacency list
void print_el(vector<vector< tuple<int, int, char> >> links); //Print adjacency and/or edge list
int initial_pop(vector<int>& group_size, vector<int>& cluster_membership, int k, int N); //create initial membership vector
vector<int> mutate_membership(vector<int> current_membership, vector<int> group_size); //mutate membership vector
long double uniform_rand(int seed); //random (0,1) draw

int likelihood1(int N,vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);
int likelihood2(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);
int modularity1(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);
int modularity2(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);
int modularity3(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);
int modularity4(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);

int likelihood10(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);
int likelihood20(int N, vector< vector< tuple<int, int, char> >> links, double InitTemp, double CR, int TL, int Max_Success);

int main_LFR(/*int argc, char * argv[]*/
	int num_nodes,//N		number of nodes
	double average_k,//k		average in-degree
	int max_degree,//maxk	maximum in-degree
	double tau,//t1		minus exponent for the degree sequence
	double tau2,//t2		minus exponent for the community size distribution
	double mixing_parameter,//mu		mixing parameter
	int nmin,//minc	minimum for the community sizes
	int nmax//maxc	maximum for the community sizes
);

int main()
{
	srand(time(0));          //seeds random number generator

/*
	//Get choice of user input file or LFR file.
	int filesource=0;
	do
	{
		std::cout << "Select network file source:" << endl;
		std::cout << "1=user file" << endl;
		std::cout << "2=generate LFR network for use" << endl;
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

	if (filesource == 1) {//if user input file, get name of file containing network and proceed as usual
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
	}

	else {//if LFR, then get parameters. An LFR named "LFR_network.dat" will be created and then read in
		//as if it were the user input file.

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

	}
*/	
	//PRINT NETWORK as adjacency list, including weights and direction
	//print_el(links); 
	//std::cin.get(); std::cin.get();

/**********************************/
/*
int num_nodes;//N		number of nodes
double average_k;//k		average in-degree
int max_degree;//maxk	maximum in-degree
double tau;//t1		minus exponent for the degree sequence -3
double tau2;//t2		minus exponent for the community size distribution -1
double mixing_parameter;//mu		mixing parameter
int nmin;//minc	minimum for the community sizes
int nmax;//maxc	maximum for the community sizes
*/

int overall_net_ct = 0;
double overall_time = 0;

for(double average_k=5;average_k<10.1; average_k=average_k+2.5)//5, 7.5, 10.0
{
	for(double tau=2; tau<3.1; tau=tau+1)//2,3
	{
		for(double tau2=1;tau2<2.1; tau2=tau2+1)//1,2
		{
			for(double mixing_parameter=0.1;mixing_parameter<0.61; mixing_parameter=mixing_parameter+0.05)//{0.10,0.15,0.20,0.25...,0.55,0.60}
			{
				for (int network_ct = 1; network_ct < 51/*101*/; network_ct++)
				{				

//Parameters for testing
//Should make 3 networks at a total of 2*1*1*2 parameter combinations =3*4=12 networks in total
//Each network is evaluated against 6 models

//Generate network at current settings
int num_nodes = 100;
int max_degree = 10;//20
int nmin = 5;//means a maximum of 100/20=5 clusters
int nmax = 20;//means a minimum of 100/25=4 clusters
			  //Thus we can restrict our search to clusters sizes of 5~10
//max_degree, nmin, nmax
//10,5,5 (11,17) (10,17) (11,18) (11,18) (11,16) (10,16) (11,17)
//20,5,5 (6,15)  (7,14)  (6,17)  (7,18)  (8,19)  (9,15)  (7,18)
//40,5,5 (3,13)  (4,14)  (3,13)  (4,13)  (3,16)  (4,13)  (3,16)
//Lesson: as max_degree goes up, k goes down and the spread gets larger
//10,5,20 (11,16) (11,18) (10,18) (11,19) (11,16) (10,18) (10,16) (9,17) (10,17) (10,19)
	//20,10,40 (6,17) (8,16) ... spread is not maintained when values are doubled
//20,5,20 (7,16) (8,18) (7,19) (7,15) (7,16) (6,15) (6,19)
//40,5,20 (3,16) (4,17)  (4,16)  (4,17)  (3,13)  (4,12)  (4,16)
//Lesson: with higher nmax, increased max_degree still reduces k and the spread gets larger than before
//10,10,20 (11,16) (10,16) (11,16) (10,16) (10,17) (11,19) (12,16)
//20,10,20 (7,20) (7,15) (8,20) (8,18) (6,20) (6,16) (7,18)
//40,10,20  (3,14) (3,16) (4,13)
//Lesson: Same as before
//10,10,10 (11,15) (11,16) (10,19) (10,16) (10,17)
//20,20,20 (8,17) (8,16) (7,17)
//10,20,20 (10,16) (10,16) (10,17) (10,16) (10,18) 

//10,25,20 ???? (9,17) (11,17) (11,17) (11,18)
//10,40,20 ???? (10,19) (9,16) 
//10,40,5 ????  (10,17) (10,19) 

//std::streamsize ss = std::cout.precision();
/*
for (double average_k = 5; average_k < 7; average_k=average_k + 1)
{
	for (double tau = 2; tau<4; tau=tau + 1)
	{
		for (double tau2 = 1; tau2 < 3; tau2=tau2 + 1)
		{
			for (double mixing_parameter = 0.1; mixing_parameter < 0.2; mixing_parameter=mixing_parameter + 0.05)
			{
				for (int network_ct = 1; network_ct < 3; network_ct++)
				{
				*/
//std::cout << fixed;
//std::cout.precision(ss);
overall_net_ct = overall_net_ct + 1;
time_t network_start = time(0);
cout<<"Generating current network, #"<< overall_net_ct <<endl;
					main_LFR(
						num_nodes,
						average_k,//avg_degree: 5-10, by 2.5
						max_degree,//max_degree,
						tau,//Gamma: 2-3, by 1
						tau2,//Beta: 1-2, by 1
						mixing_parameter,//myu: 0.1-0.6, by 0.05 
						 //These values below are estimated from fig. 10, pg 183, of original papger for N=100 instead of N=500
						nmin,//minimum community size 100/4~approx 25 clusters
						nmax//maximum community size 100/10~approx 10 clusters
						);

cout << "Currently evaluating network with parameters: " << num_nodes << "," << average_k << "," << max_degree << "," << tau << "," << tau2 << "," << mixing_parameter << "," << nmin << "," << nmax << endl;
					ostringstream PARAMETERS_NAME;
					ofstream NETWORK_PARAMETERS;
					PARAMETERS_NAME << "LFR_PARAMETERS" << ".txt";
					NETWORK_PARAMETERS.open("LFR_PARAMETERS.txt", ios::app);
					NETWORK_PARAMETERS << num_nodes << "," << average_k << "," << max_degree << "," << tau << "," << tau2 << "," << mixing_parameter << "," << nmin << "," << nmax;
					NETWORK_PARAMETERS << endl;
					NETWORK_PARAMETERS.close();

					//Read in current network
					char filename[256];		
					strcpy_s(filename, "LFR_network.dat");
					ifstream ifile(filename);

					vector< vector< tuple<int, int, char> >> links;
					char weighted = 'N';
					links = read_el(filename, weighted);

					int N = 0;
					int E = 0;
					//Get the number of distinct vertices in the network
					N = links.size();
					//Get number of distinct edges in the network
					/*for (int i = 0; i < links.size(); i++)
					{
						for (int j = 0; j < links[i].size(); j++) {
							if (get<2>(links[i][j]) == 'T') {
								E++;
							}
						}
					}*/
					//Simulated Annealing Algorithm Parameters
					double InitTemp = 1;	//Starting temperature
					double CR = 0.99;		//Rate at which temperature is reduced

					int TL = 10;			//Maximum number of reclustering attempts at a given temperature
					int Max_Success = 100;	//Maximum number of successes allowed at a given temperature 
					//double LimitIT = 1.0e-8; //Minimum temperature. Program stops when this number is reached
					//long double epsilon = 1.0e-6; //Minimum change in loglikelihood required to accept proposed membership mutation
//for k=10 to k=20
//note at CR=0.99, TL=300, MX=100, evaluating 4 algos on 1 network took 24 minutes. 24*132 settings*100 reps= 316800 min/60 min= 5280 hrs/24 hrs = 220 days...impossible
//at CR=0.9, TL=300, MX=100, evaluating 4 algos on 1 network took 3 minutes=660 hrs=27.5 days
//note at CR=0.99, TL=100, MX=100, evaluating 4 algos on 1 network took 8.4 min=110880 min=1,848 hrs=77 days
	//2 also would then take 38.5 days. reducing TL=50 should(?) reduce it to around 20 days.
	//Also reducing the number of networks per setting to 50 should then be around 10 days.
//note at CR=0.99, TL=10, MX=100, evaluating 4 algos on 1 network took 1 min=13200 min=220hr=9.1 days
	//2 algos would then take 4.55 days. Reducing the number of networks per setting to 50 should then be about 2.75 days.
//note at CR=0.9, TL=100, MX=100, evaluating 4 algos on 1 network took .8 min=10560 min=176 hrs=7.3 days
//note at CR=0.9, TL=10, MX=100, evaluating 4 algos on 1 network took .13 min. =1716 min = 28.6 hrs <2 days

					//Run all models on current network
//cout << "\t" << "Evaluating Model: likelihood1" <<endl;
//					likelihood1(N, links, InitTemp, CR, TL, Max_Success); 
cout << "\t" << "Evaluating Model: likelihood2" <<endl;
					likelihood2(N, links, InitTemp, CR, TL, Max_Success); 
					//	modularity1(N, links, InitTemp, CR, TL, Max_Success); 
					//	modularity2(N, links, InitTemp, CR, TL, Max_Success);
//cout << "\t" << "Evaluating Model: modularity3" <<endl; 
//					modularity3(N, links, InitTemp, CR, TL, Max_Success); 
cout << "\t" << "Evaluating Model: modularity4" <<endl;
					modularity4(N, links, InitTemp, CR, TL, Max_Success); 


//cout << "\t" << "Evaluating Model: likelihood10" <<endl;
//					likelihood10(N, links, InitTemp, CR, TL, Max_Success);
//cout << "\t" << "Evaluating Model: likelihood20" <<endl; 
//					likelihood20(N, links, InitTemp, CR, TL, Max_Success); 

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

/*
	//Simulated Annealing Algorithm Parameters
	double InitTemp = 1;	//Starting temperature
	double CR = 0.99;		//Rate at which temperature is reduced
	int TL = 300;			//Maximum number of reclustering attempts at a given temperature
	int Max_Success = 300;	//Maximum number of successes allowed at a given temperature 
							//Note: Setting both these values=1 is the equivalent of not using them at all in the SA algorithm.
	double LimitIT = 1.0e-8; //Minimum temperature. Program stops when this number is reached
	long double epsilon = 1.0e-6; //Minimum change in loglikelihood required to accept proposed membership mutation




	std::cout << endl;
	int SA_choice;
	do
	{
		std::cout << "Enter SA cooling schedule or use defaults of:" << endl;
		std::cout << "\t"<<"Initial Temperature=" << InitTemp << endl;
		std::cout << "\t" << "Cooling Rate=" << CR << endl;
		std::cout << "\t" << "Temperature Length=" << TL << endl;
		std::cout << "\t" << "Max # of successes allowed at current temp=" << Max_Success << endl;
		std::cout << "Enter 1 to use default SA cooling schedule or 2 to enter your own:" << endl;
		std::cin >> SA_choice;
		//cout << "SA_choice is: " << SA_choice << endl;
	} while (
		!std::cin.fail()
		&& SA_choice != 1
		&& SA_choice != 2
		);
	if (SA_choice == 1) {
		std::cout << "Using default SA cooling schedule..." << endl;
		cout << endl;
	}
	else if (SA_choice == 2) {
		std::cout << "Enter the initial temperature:" << endl;
		std::cin >> InitTemp;
		std::cout << endl;
		std::cout << "Enter the cooling rate:" << endl;
		std::cin >> CR;
		std::cout << endl;
		std::cout << "Enter the temperature length:" << endl;
		std::cin >> TL;
		std::cout << endl;
		std::cout << "Enter max # of successes allowed at current temp:" << endl;
		std::cin >> Max_Success;
		std::cout << endl;
	}
	// double InitTemp, double CR, int TL, int Max_Success
	//InitTemp, CR, TL, Max_Success
	int method;
	do
	{
		std::cout << "Enter optimization objective function:" << endl;
		std::cout << "1=likelihood (iterative, full): k within-cluster densities and k(k-1) between-cluster densities" << endl<< endl;
		std::cout << "2=likelihood (iterative, full): k within-cluster densities and 1 between-cluster density" << endl << endl;
//		std::cout << "3=modularity (iterative, abbrev.): BIC/AIC evaluated using #1 likelihood model" << endl<< endl;
//		std::cout << "4=modularity (iterative, abbrev.): BIC/AIC evaluated using #2 likelihood model" << endl<< endl;
		std::cout << "5=modularity (full): BIC/AIC evaluated using #1 likelihood model" << endl << endl;
		std::cout << "6=modularity (full): BIC/AIC evaluated using #2 likelihood model" << endl << endl;
		std::cout << "10=likelihood (non-iterative, full): k within-cluster densities and k(k-1) between-cluster densities" << endl << endl;
		std::cout << "20=likelihood (non-iterative, full): k within-cluster densities and 1 between-cluster density" << endl << endl;
		std::cin >> method;
		//cout << "method is: " << method << endl;
	} while (
		!std::cin.fail() 
		&& method != 1
		&& method != 2
//		&& method != 3
//		&& method != 4
		&& method != 5
		&& method != 6
		&& method != 10
		&& method != 20
		);
	if (method == 1) { likelihood1(N, links, InitTemp, CR, TL, Max_Success); }
	else if (method == 2) { likelihood2(N, links, InitTemp, CR, TL, Max_Success); }
//	else if (method == 3) { modularity1(N, links, InitTemp, CR, TL, Max_Success); }
//	else if (method == 4) { modularity2(N, links, InitTemp, CR, TL, Max_Success); }
	else if (method == 5) { modularity3(N, links, InitTemp, CR, TL, Max_Success); }
	else if (method == 6) { modularity4(N, links, InitTemp, CR, TL, Max_Success); }
	else if (method == 10) { likelihood10(N, links, InitTemp, CR, TL, Max_Success); }
	else if (method == 20) { likelihood20(N, links, InitTemp, CR, TL, Max_Success); }
*/

	/*Original Code Removed From Here. See Bottom of this File.*/

	std::cout << "Program ended. Hit enter to exit.";
	std::cin.get(); std::cin.get();

}//END main()


//FUNCTION: Returns a U(0,1) variable
long double uniform_rand(int seed)
{

	return seed * (1.0 / (RAND_MAX + 1.0));     //to call the random U(0,1), double r = uniform_rand ( rand() );
}//end uniform_rand function definition

//Below is the original code that was removed...

//	int STG1_CT = 0;
//	int STG2_CT = 0;
//	int STG3_CT = 0;
 //	//Get the number of clusters in which to cluster the network
 //	int min_k;
 //	int max_k;
 //	std::cout << "Enter the lower bound for number of clusters: "; 	std::cin >> min_k;
 //	std::cout << "Enter the upper bound for number of clusters: ";	std::cin >> max_k;
 //
 //	//Simulated Annealing Algorithm Parameters
 //	double InitTemp = 100;	//Starting temperature
 //	double CR = 0.99;		//Rate at which temperature is reduced
 //	int TL = 300;			//Maximum number of reclustering attempts at a given temperature
 //	int Max_Success = 100;	//Maximum number of successes allowed at a given temperature 
 //	//Note: Setting both these values=1 is the equivalent of not using them at all in the SA algorithm.
 //	double LimitIT = 1.0e-8; //Minimum temperature. Program stops when this number is reached
 //	long double epsilon = 1.0e-6; //Minimum change in loglikelihood required to accept proposed membership mutation
 //	int TL_counter = 0; //Number of iterations at the current temperature
 //	int Success_counter = 0; //Number of proposed membership mutations accepted at the current temperature
 //
 //	double IT = InitTemp;  //Starting temperature. Reset for every new value of k
 //
 //	//Matrices and vectors to store final cluster assignments, BIC, loglikelihood and test statistics for each k
 //	vector< vector <int >> MEMBERSHIP; //final cluster assignment for each k  
 //	vector<long double> BIC; //BIC for each k's final cluster assignment
 //	vector<long double> LOGLIK;  //loglikelihood for each k's final cluster assignment
 //	vector<long double> TESTSTAT;  //test statistic for each k's final cluster assignment
 //
 //	//Vectors and variables to hold overall "best" cluster assignment and corresponding k value, BIC, loglikelihood and test statistic
 //	int final_k = 0;
 //	long double final_BIC = 0;
 //	long double final_LOGLIK = 0;
 //	long double final_TESTSTAT = 0;
 //	vector<int> final_membership;
 //
 //	time_t SA_start = time(0);//Time until solution
 //
 //	//Repeate Simulated Annealing process for each k in the range of specified cluster numbers
 //	for (int k = min_k; k <= max_k; k++)
 //	{
 //		int overall_iteration_count = 0;
 //
 //		//Initialize POS and OBS matrices (kxk vector of vectors)
 //		vector< vector <unsigned long long int >> OBS(k, vector<unsigned long long int>(k, 0));
 //		vector< vector <unsigned long long int >> POS(k, vector<unsigned long long int>(k, 0));
 //
 //		//Initialize and populate cluster assignment and size vectors
 //		vector<int> group_size(k, 0);//used in initial pop to store number of nodes in each group 1:k
 //		vector<int> cluster_membership;	//holds cluster assignments of N vertices
 //
 //		//vertex to be reassigned, the cluster in which it is currently assigned, and the cluster to which it will be reassigned
 //		vector<int> reassignment;
 //		int l_star;		//The vertex to be reassigned (0:[N-1])
 //		int clust_l;	//L-star's cluster membership before reassignment (1:k)
 //		int clust_j;	//L-star's cluster membership after reassignment (1:k)
 //
 //		//Total weight of links to l*, by cluster (These are the Z_(to *) variables in the paper)
 //		vector<int> Z_to(k);
 //		//Total weight of links from l*, by cluster (These are the Z_(from *) variables in the paper)
 //		vector<int> Z_from(k);
 //
 //		//Partial loglikelihood=the components of the loglikelihood that (may) change after a reclustering
 //		long double l_0; //partial loglikelihood for current clustering
 //		long double l_1; //partial loglikelihood for proposed clustering
 //		long double delta_l; //change-in-loglikelihood between current clustering and proposed clustering
 //
 //		//Densities in the denominator in the change-in-loglik formula [for current clustering]
 //		long double l_density_0 = 0; //density of edges in cluster l
 //		long double j_density_0 = 0; //density of edges in cluster j
 //		long double lj_density_0 = 0; //density of edges from cluster l to cluster j
 //		long double jl_density_0 = 0; //density of edges from cluster j to cluster l
 //
 //		//Densities in the denomonator in the change-in-loglike formula [for current clustering]
 //		//these are to/from clusters l and j to cluster "star", where star is an index star=1,2,...k, exclusive of {j,l}
 //		long double starl_density_0 = 0; 
 //		long double lstar_density_0 = 0;
 //		long double starj_density_0 = 0;
 //		long double jstar_density_0 = 0;
 //
 //		//Densities in the numerator in the change-in-loglik formula [for proposed clustering]
 //		long double l_density_1 = 0;
 //		long double j_density_1 = 0;
 //		long double lj_density_1 = 0;
 //		long double jl_density_1 = 0;
 //
 //		//Densities in the numerator in the change-in-loglike formula [for proposed clustering]
 //		//these are to/from clusters l and j to cluster "star", where star is an index star=1,2,...k, exclusive of {j,l}
 //		long double starl_density_1 = 0;
 //		long double lstar_density_1 = 0;
 //		long double starj_density_1 = 0;
 //		long double jstar_density_1 = 0;
 //
 //		//Flag to allow creation of initial clustering (at the current k) that will be mutated on later runs.
 //		int first_run = 1;
 //
 //		//Flag that indicates we should retain the current solution likelihood if we reject the proposed solution, so the current solution likelihood need not be recalculated
 //		int retain_current_sol = 0;
 //
 //		//Reset initial temperature at the start of each k-loop
 //		IT = InitTemp;
 //		//Reset count length at the start of each k-loop
 //		TL_counter = 0;
 //		//Reset success counter at the start of each k-loop
 //		Success_counter = 0;
 //
 //		//START Simulated Annealing Algorithm
 //		while (IT > LimitIT)
 //		{
 //			TL_counter++; //Iterate number of attempts at current temperature
 //
 //			//After TL attempts at reclustering at the current temperature, or Max_Success # of successes at curren temperature
 //			//Slowly decrease temperature until it reaches lower limit
 //			if (TL_counter >= TL || Success_counter >= Max_Success)
 //			{
 //				IT = CR*IT; //decrease current temperature
 //				TL_counter = 0; //reset length counter
 //				Success_counter = 0; //reset success counter
 //			}
 //			
 //			//Generate initial clustering (1st run only). This is the initial "current clustering".
 //			//On runs >1, the "current clustering" will either be the one generated here, or a proposed clustering that was accepted below.
 //			while (first_run == 1)
 //			{
 //
 //				//Populate cluster assignment and size vectors
 //				initial_pop(group_size, cluster_membership, k, N);//assigns initial cluster assignments to z, and counts cluster sizes
 //				first_run = 0;
 //
 //				//Populate observed OBS edge weight vector. This requires 1 full traversal of adjacency list
 //				for (int i = 0; i < links.size(); i++)
 //				{
 //					for (int j = 0; j < links[i].size(); j++)
 //					{
 //						if (get<2>(links[i][j]) == 'T')//edge is FROM src vert links[i] TO dest vert links[i][j]...
 //						{
 //							if (cluster_membership[i] == cluster_membership[get<0>(links[i][j])])//src/dest verts belong to same cluster 
 //							{
 //								OBS[cluster_membership[i] - 1][cluster_membership[i] - 1]
 //									= OBS[cluster_membership[i] - 1][cluster_membership[i] - 1]
 //									+ get<1>(links[i][j]);
 //							}
 //							else//src/dest clusters differ
 //							{
 //								OBS[cluster_membership[i] - 1][cluster_membership[get<0>(links[i][j])] - 1]
 //									= OBS[cluster_membership[i] - 1][cluster_membership[get<0>(links[i][j])] - 1]
 //									+ get<1>(links[i][j]);
 //							}
 //						}
 //					}
 //				}
 //
 //				//Populate possible edge count vector. This is done using the cluster size information
 //				for (int i = 0; i < k; i++)
 //				{
 //					POS[i][i] = ((unsigned long long int) group_size[i] * (group_size[i] - 1)); //2*[nchoose2]=2*[N(N-1)/2]=N*N-1
 //				}
 //				for (int i = 0; i < k - 1; i++)
 //				{
 //					for (int j = i + 1; j < k; j++)
 //					{
 //						//add number of edges possible between cluster i and cluster j
 ////This is only run once, on the first run. Since each i/j combination is only considered once, can't we just shorten this to?:
 ////POS[i][j] = (unsigned long long int) group_size[i] * group_size[j];
 ////But keeping it as-is shouldn't hurt anything either...
 //						POS[i][j] = POS[i][j] + (unsigned long long int) group_size[i] * group_size[j];
 //						POS[j][i] = POS[i][j];
 //					}
 //				}
 //			}//end 1st-run WHILE loop
 //
 //			//Propose a new clustering
 //			//Select Random Vertex l* for reassignment and select cluster to which it will be moved
 //			//This returns a) The vertex to be reassigned and b) the cluster to which it will be moved.
 //			reassignment = mutate_membership(cluster_membership, group_size);	
 //			l_star = reassignment[0];				//The vertex to be reassigned
 //			clust_l = cluster_membership[l_star];	//l-star's cluster membership before reassignment
 //			clust_j = reassignment[1];				//l-star's cluster membership after reassignment
 //
 //			//Get denominator of change-in-loglik equation [l,j,bw for current clustering]
 //
 //			
 //			if (retain_current_sol == 0) //if 0, then calculate likelihood for new solution. if 1, likelihood for previous solution is retained
 //			{
 //				//Calculate l,j,bw densities for current clustering, the denominator in the change-in-loglik formula
 //				//if unweighted, all edges have weight=1 and Binomial density is appropriate
 //				if (weighted == 'N')
 //				{
 //					//Can read observed edge weight sums in the usual manner
 //					//std::cout << OBS[1][2] << "..." << OBS[2][1]<<endl << endl;
 //
 //					//cluster l density
 //					l_density_0 = 0;
 //					long double param_l_0 = OBS[clust_l - 1][clust_l - 1] / (long double)POS[clust_l - 1][clust_l - 1];
 //					if (param_l_0 == 0 || param_l_0 == 1){ l_density_0 = 0; }
 //					else{
 //						l_density_0 = (OBS[clust_l - 1][clust_l - 1])*log(param_l_0) + (POS[clust_l - 1][clust_l - 1] - OBS[clust_l - 1][clust_l - 1])*log(1 - param_l_0);
 //					}
 //					//cluster j density
 //					j_density_0 = 0;
 //					long double param_j_0 = OBS[clust_j - 1][clust_j - 1] / (long double)POS[clust_j - 1][clust_j - 1];
 //					if (param_j_0 == 0 || param_j_0 == 1){ j_density_0 = 0; }
 //					else{
 //						j_density_0 = (OBS[clust_j - 1][clust_j - 1])*log(param_j_0) + (POS[clust_j - 1][clust_j - 1] - OBS[clust_j - 1][clust_j - 1])*log(1 - param_j_0);
 //					}
 //					//b/w cluster density indexed by l and j (row l, col j)
 //					lj_density_0 = 0;
 //					long double param_lj_0 = OBS[clust_l - 1][clust_j - 1] / (long double)POS[clust_l - 1][clust_j - 1];
 //					if (param_lj_0 == 0 || param_lj_0 == 1){ lj_density_0 = 0; }
 //					else{
 //						lj_density_0 = (OBS[clust_l - 1][clust_j - 1])*log(param_lj_0) + (POS[clust_l - 1][clust_j - 1] - OBS[clust_l - 1][clust_j - 1])*log(1 - param_lj_0);
 //					}
 //					//b/w cluster density indexed by l and j (row j, col l)
 //					jl_density_0 = 0;
 //					long double param_jl_0 = OBS[clust_j - 1][clust_l - 1] / (long double)POS[clust_j - 1][clust_l - 1];
 //					if (param_jl_0 == 0 || param_jl_0 == 1){ jl_density_0 = 0; }
 //					else{
 //						jl_density_0 = (OBS[clust_j - 1][clust_l - 1])*log(param_jl_0) + (POS[clust_j - 1][clust_l - 1] - OBS[clust_j - 1][clust_l - 1])*log(1 - param_jl_0);
 //					}
 //
 //					//From here down: "star" is a general index for cluster combination number. For example, starl means *l, for *=1,2,...,(k-1),k
 //					//density of edges from (clusters not indexed by l or j) to cluster l
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by column l  (row star, col l)
 //					starl_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_starl_0 = OBS[star - 1][clust_l - 1] / (long double)POS[star - 1][clust_l - 1];
 //							//If no edges or complete edges between cluster star and cluster l, set density=0 => don't add anything to starl_density_0
 //							if (param_starl_0 == 0 || param_starl_0 == 1){ starl_density_0 = starl_density_0; }
 //							else{
 //								//else we're adding the sum of the logs of the densities
 //								starl_density_0 = starl_density_0 + (OBS[star - 1][clust_l - 1])*log(param_starl_0) + (POS[star - 1][clust_l - 1] - OBS[star - 1][clust_l - 1])*log(1 - param_starl_0);
 //							}
 //						}
 //
 //					}
 //					//density of edges from cluster l to (clusters not indexed by l or j)
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by row l  (row l, col star)
 //					lstar_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_lstar_0 = OBS[clust_l - 1][star - 1] / (long double)POS[clust_l - 1][star - 1];
 //							//If no edges or complete edges between cluster star and cluster l, set density=0 => don't add anything to lstar_density_0
 //							if (param_lstar_0 == 0 || param_lstar_0 == 1){ lstar_density_0 = lstar_density_0; }
 //							else{
 //								lstar_density_0 = lstar_density_0 + (OBS[clust_l - 1][star - 1])*log(param_lstar_0) + (POS[clust_l - 1][star - 1] - OBS[clust_l - 1][star - 1])*log(1 - param_lstar_0);
 //							}
 //						}
 //					}
 //					//density of edges from (clusters not indexed by l or j) to cluster j
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by column j  (row star, col j)
 //					starj_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_starj_0 = OBS[star - 1][clust_j - 1] / (long double)POS[star - 1][clust_j - 1];
 //							//If no edges or complete edges between cluster star and cluster j, set density=0 => don't add anything to starj_density_0
 //							if (param_starj_0 == 0 || param_starj_0 == 1){ starj_density_0 = starj_density_0; }
 //							else{
 //								starj_density_0 = starj_density_0 + (OBS[star - 1][clust_j - 1])*log(param_starj_0) + (POS[star - 1][clust_j - 1] - OBS[star - 1][clust_j - 1])*log(1 - param_starj_0);
 //							}
 //						}
 //					}
 //					//density of edges from cluster j to (clusters not indexed by l or j) 
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by row l  (row j, col star)
 //					jstar_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_jstar_0 = OBS[clust_j - 1][star - 1] / (long double)POS[clust_j - 1][star - 1];
 //							//If no edges or complete edges between cluster star and cluster j, set density=0 => don't add anything to jstar_density_0
 //							if (param_jstar_0 == 0 || param_jstar_0 == 1){ jstar_density_0 = jstar_density_0; }
 //							else{
 //								jstar_density_0 = jstar_density_0 + (OBS[clust_j - 1][star - 1])*log(param_jstar_0) + (POS[clust_j - 1][star - 1] - OBS[clust_j - 1][star - 1])*log(1 - param_jstar_0);
 //							}
 //						}
 //					}
 //
 //					//Denominator of modified change-in-loglikelihood equation
 //					l_0 = starl_density_0 + lstar_density_0 + starj_density_0 + jstar_density_0
 //						+ l_density_0 + j_density_0 + lj_density_0 + jl_density_0;
 //
 ////Just want to see if using a single density for ALL between-cluster edges improves performance		
 ///*
 //					unsigned long long BW_OBS = 0;
 //					unsigned long long BW_POS = 0;
 //					for (int i = 0; i < k; i++) 
 //					{for (int j = 0; j < k; j++) 
 //						{if (i != j) 
 //							{BW_OBS = BW_OBS + OBS[i][j]; BW_POS = BW_POS + POS[i][j];}
 //						}						
 //					}
 //					l_0= l_density_0 + j_density_0 + ((long double)BW_OBS)*log((long double)BW_OBS/ BW_POS) + ((long double)BW_POS - BW_OBS)*log(1 - (long double)BW_OBS / BW_POS);
 //*/
 //				}
 //
 //				//If weighted, all edges are non-negative integers and Poisson density is appropriate
 //				else
 //				{
 //					//Noe: Parameter estimates can - in theory - take on values of zero but poisson only takes parameter >0.
 //					//If a density's parameter estimate is 0, it will ruin the loglik calculation.
 //
 //					//cluster l density
 //					l_density_0 = 0;
 //					long double param_l_0 = OBS[clust_l - 1][clust_l - 1] / (long double)POS[clust_l - 1][clust_l - 1];
 //					if (param_l_0 == 0){ l_density_0 = 0; }
 //					else{
 //						l_density_0 = OBS[clust_l - 1][clust_l - 1] * log(param_l_0) - param_l_0;
 //					}
 //					//cluster j density
 //					j_density_0 = 0;
 //					long double param_j_0 = OBS[clust_j - 1][clust_j - 1] / (long double)POS[clust_j - 1][clust_j - 1];
 //					if (param_j_0 == 0){ j_density_0 = 0; }
 //					else{
 //						j_density_0 = OBS[clust_j - 1][clust_j - 1] * log(param_j_0) - param_j_0;
 //					}
 //					//b/w cluster density indexed by l and j (row l, col j)
 //					lj_density_0 = 0;
 //					long double param_lj_0 = OBS[clust_l - 1][clust_j - 1] / (long double)POS[clust_l - 1][clust_j - 1];
 //					if (param_lj_0 == 0){ lj_density_0 = 0; }
 //					else{
 //						lj_density_0 = OBS[clust_l - 1][clust_j - 1] * log(param_lj_0) - param_lj_0;
 //					}
 //					//b/w cluster density indexed by l and j (row j, col l)
 //					jl_density_0 = 0;
 //					long double param_jl_0 = OBS[clust_j - 1][clust_l - 1] / (long double)POS[clust_j - 1][clust_l - 1];
 //					if (param_jl_0 == 0){ jl_density_0 = 0; }
 //					else{
 //						jl_density_0 = OBS[clust_j - 1][clust_l - 1] * log(param_jl_0) - param_jl_0;
 //					}
 //
 //					//From here down: "star" is a general index for cluster combination number. For example, starl means *l, for *=1,2,...,(k-1),k
 //					//density of edges from (clusters not indexed by l or j) to cluster l
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by column l  (row star, col l)
 //					starl_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_starl_0 = OBS[star - 1][clust_l - 1] / (long double)POS[star - 1][clust_l - 1];
 //							//If no edges between cluster star and cluster j, set density=0 => don't add anything to jstar_density_0
 //							if (param_starl_0 == 0){ starl_density_0 = starl_density_0; }
 //							else{
 //								starl_density_0 = starl_density_0 + OBS[star - 1][clust_l - 1] * log(param_starl_0) - param_starl_0;
 //							}
 //						}
 //					}
 //					//density of edges from cluster l to (clusters not indexed by l or j)
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by row l  (row l, col star)
 //					lstar_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_lstar_0 = OBS[clust_l - 1][star - 1] / (long double)POS[clust_l - 1][star - 1];
 //							//If no edges between cluster star and cluster j, set density=0 => don't add anything to jstar_density_0
 //							if (param_lstar_0 == 0){ lstar_density_0 = lstar_density_0; }
 //							else{
 //								lstar_density_0 = lstar_density_0 + OBS[clust_l - 1][star - 1] * log(param_lstar_0) - param_lstar_0;
 //							}
 //						}
 //					}
 //					//density of edges from c(clusters not indexed by l or j) to cluster j
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by column j  (row star, col j)		
 //					starj_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_starj_0 = OBS[star - 1][clust_j - 1] / (long double)POS[star - 1][clust_j - 1];
 //							//If no edges between cluster star and cluster j, set density=0 => don't add anything to jstar_density_0
 //							if (param_starj_0 == 0){ starj_density_0 = starj_density_0; }
 //							else{
 //								starj_density_0 = starj_density_0 + OBS[star - 1][clust_j - 1] * log(param_starj_0) - param_starj_0;
 //							}
 //						}
 //					}
 //					//density of edges from cluster j to (clusters not indexed by l or j)
 //					//sum of (log densities)= log of (product densities)
 //					//cluster densities indexed by row j  (row j, col star)
 //					jstar_density_0 = 0;
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if (star != clust_l && star != clust_j)
 //						{
 //							long double param_jstar_0 = OBS[clust_j - 1][star - 1] / (long double)POS[clust_j - 1][star - 1];
 //							//If no edges between cluster star and cluster j, set density=0 => don't add anything to jstar_density_0
 //							if (param_jstar_0 == 0){ jstar_density_0 = jstar_density_0; }
 //							else{
 //								jstar_density_0 = jstar_density_0 + OBS[clust_j - 1][star - 1] * log(param_jstar_0) - param_jstar_0;
 //							}
 //						}
 //					}
 //
 //					//Denominator of modified change-in-loglikelihood equation		
 //					l_0 = starl_density_0 + lstar_density_0 + starj_density_0 + jstar_density_0
 //						+ l_density_0 + j_density_0 + lj_density_0 + jl_density_0;
 //				}
 //			}//end if(retain_current_sol == 0)
 //
 //		//Get numerator of modified change-in-loglikelihood equation for proposed clustering
 //
 //			
 //			//Zero out Z (to/from cluster size) vectors. Best practice for runs > 1			
 //			for (int i = 0; i < Z_to.size(); i++){ Z_to[i] = 0; }//Total weight of links to l*, by cluster (These are the Z_(to s) values in the paper)			
 //			for (int i = 0; i < Z_from.size(); i++){ Z_from[i] = 0; }//Total weight of links from l*, by cluster (These are the Z_(from s) values in the paper)
 //
 //			//Populate to/from cluster size vectors
 //			for (int j = 0; j < links[l_star].size(); j++)
 //			{
 //				//get<0>(links) -Destination vertex
 //				//get<1>(links) -Weight
 //				//get<2>(links) -char in {T,F}
 //				if (get<2>(links[l_star][j]) == 'T')//edge extends FROM vertex l_star TO the jth vertex in l_star's adjacency list row
 //				{
 //					Z_to[cluster_membership[get<0>(links[l_star][j])] - 1] = Z_to[cluster_membership[get<0>(links[l_star][j])] - 1] + get<1>(links[l_star][j]);
 //				}
 //				else//edge extends FROM jth vertex in l_star's adjacency list row TO vertex l_star
 //				{
 //					Z_from[cluster_membership[get<0>(links[l_star][j])] - 1] = Z_from[cluster_membership[get<0>(links[l_star][j])] - 1] + get<1>(links[l_star][j]);
 //				}
 //			}
 ///*
 //			std::cout << "Cluster memberships:" << endl;
 //			for (int i = 0; i < cluster_membership.size(); i++){ std::cout << cluster_membership[i] << "   "; } std::cout << endl;
 //
 //			for (int i = 1; i <= k; i++){
 //			std::cout << endl << "The number of edges from vertex " << vert << " to cluster " << i << " is " << Z_to[i - 1] << endl;}
 //			for (int i = 1; i <= k; i++){
 //			std::cout << endl << "The number of edges from cluster " << i << " to vertex " << vert << " is " << Z_from[i - 1] << endl;}
 //
 //			std::cout << endl << "Z_to():" << endl;
 //			for (int i = 0; i < Z_to.size(); i++){std::cout << "<" << i + 1 << ", " << Z_to[i] << "> ";}
 //
 //			std::cout << endl << "Z_from():" << endl;
 //			for (int i = 0; i < Z_from.size(); i++){std::cout << "<" << i + 1 << ", " << Z_from[i] << "> ";}
 //*/
 //			//Calculate l,j,bw densities for current clustering, the denominator in the change-in-loglik formula
 //			//NEED TO TEMPORARILY UPDATE OBS and POS elements required to calculate the denominator in the change-in-loglik formula
 //			//IF we accept the proposed solution, these will become the new OBS/POS matrices. Otherwise, we'll discare and try another solution
 //			vector< vector <unsigned long long int >> OBS_TEMP(k, vector<unsigned long long int>(k, 0)); //Observed edge weight matrix for the proposed clustering
 //			vector< vector <unsigned long long int >> POS_TEMP(k, vector<unsigned long long int>(k, 0)); //Possible edge count matrix for the proposed clustering
 //			//How can I automatically populate these matrices with the contents of OBS/POS rather than first setting
 //			//them to 0 and then copying OBS/POS??
 //			OBS_TEMP = OBS;
 //			POS_TEMP = POS; 
 //
 //			//Update OBS and POS matrix elements for proposed clustering strictly associated ONLY with elements IN {l,j}
 //			OBS_TEMP[clust_l - 1][clust_l - 1] = OBS[clust_l - 1][clust_l - 1] - Z_to[clust_l - 1] - Z_from[clust_l - 1];
 //			OBS_TEMP[clust_j - 1][clust_j - 1] = OBS[clust_j - 1][clust_j - 1] + Z_to[clust_j - 1] + Z_from[clust_j - 1];
 //			OBS_TEMP[clust_l - 1][clust_j - 1] = OBS[clust_l - 1][clust_j - 1] - Z_to[clust_j - 1] + Z_from[clust_l - 1];
 //			OBS_TEMP[clust_j - 1][clust_l - 1] = OBS[clust_j - 1][clust_l - 1] + Z_to[clust_l - 1] - Z_from[clust_j - 1];
 //			//group_size[i-1] is the size of cluster i
 //			POS_TEMP[clust_l - 1][clust_l - 1] = (group_size[clust_l - 1] - 1 )*(group_size[clust_l - 1] - 2);
 //			POS_TEMP[clust_j - 1][clust_j - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_j - 1]);
 //			POS_TEMP[clust_l - 1][clust_j - 1] = (group_size[clust_l - 1] - 1)*(group_size[clust_j - 1] + 1);
 //			POS_TEMP[clust_j - 1][clust_l - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_l - 1] - 1);
 //
 //			//Update OBS and POS matrix elements for proposed clustering that are associated with elements {l,j} AND another element NOT IN {l,j}
 //			for (int star = 1; star <= k; star++)
 //			{
 //				if ((star != clust_l) && (star != clust_j))
 //				{
 //					OBS_TEMP[star - 1][clust_l - 1] = OBS[star - 1][clust_l - 1] - Z_from[star - 1];
 //					OBS_TEMP[clust_l - 1][star - 1] = OBS[clust_l - 1][star - 1] - Z_to[star - 1];
 //					OBS_TEMP[star - 1][clust_j - 1] = OBS[star - 1][clust_j - 1] + Z_from[star - 1];
 //					OBS_TEMP[clust_j - 1][star - 1] = OBS[clust_j - 1][star - 1] + Z_to[star - 1];
 //
 //					POS_TEMP[star - 1][clust_l - 1] = (group_size[star - 1])*(group_size[clust_l - 1] - 1);
 //					POS_TEMP[clust_l - 1][star - 1] = (group_size[clust_l - 1] - 1)*(group_size[star - 1]);
 //					POS_TEMP[star - 1][clust_j - 1] = (group_size[star - 1])*(group_size[clust_j - 1] + 1);
 //					POS_TEMP[clust_j - 1][star - 1] = (group_size[clust_j - 1] + 1)*(group_size[star - 1]);
 //				}
 //			}
 //
 //			//if unweighted, all edges have weight=1 and Binomial density is appropriate
 //			if (weighted == 'N')
 //			{
 //				//cluster l density
 //				l_density_1 = 0;
 //				long double param_l_1 = OBS_TEMP[clust_l - 1][clust_l - 1] / (long double)POS_TEMP[clust_l - 1][clust_l - 1];
 //				if (param_l_1 == 0 || param_l_1 == 1){ l_density_1 = 0; }
 //				else{
 //					l_density_1 = (OBS_TEMP[clust_l - 1][clust_l - 1])*log(param_l_1) + (POS_TEMP[clust_l - 1][clust_l - 1] - OBS_TEMP[clust_l - 1][clust_l - 1])*log(1 - param_l_1);
 //				}
 //				//cluster j density
 //				j_density_1 = 0;
 //				long double param_j_1 = OBS_TEMP[clust_j - 1][clust_j - 1] / (long double)POS_TEMP[clust_j - 1][clust_j - 1];
 //				if (param_j_1 == 0 || param_j_1 == 1){ j_density_1 = 0; }
 //				else{
 //					j_density_1 = (OBS_TEMP[clust_j - 1][clust_j - 1])*log(param_j_1) + (POS_TEMP[clust_j - 1][clust_j - 1] - OBS_TEMP[clust_j - 1][clust_j - 1])*log(1 - param_j_1);
 //				}
 //				//b/w cluster density indexed by l and j (row l, col j)
 //				lj_density_1 = 0;
 //				long double param_lj_1 = OBS_TEMP[clust_l - 1][clust_j - 1] / (long double)POS_TEMP[clust_l - 1][clust_j - 1];
 //				if (param_lj_1 == 0 || param_lj_1 == 1){ lj_density_1 = 0; }
 //				else{
 //					lj_density_1 = (OBS_TEMP[clust_l - 1][clust_j - 1])*log(param_lj_1) + (POS_TEMP[clust_l - 1][clust_j - 1] - OBS_TEMP[clust_l - 1][clust_j - 1])*log(1 - param_lj_1);
 //				}
 //				//b/w cluster density indexed by j and l (row j, col l)
 //				jl_density_1 = 0;
 //				long double param_jl_1 = OBS_TEMP[clust_j - 1][clust_l - 1] / (long double)POS_TEMP[clust_j - 1][clust_l - 1];
 //				if (param_jl_1 == 0 || param_jl_1 == 1){ jl_density_1 = 0; }
 //				else{
 //					jl_density_1 = (OBS_TEMP[clust_j - 1][clust_l - 1])*log(param_jl_1) + (POS_TEMP[clust_j - 1][clust_l - 1] - OBS_TEMP[clust_j - 1][clust_l - 1])*log(1 - param_jl_1);
 //				}
 //
 //				//cluster densities indexed by column l  (row star, col l)
 //				starl_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_starl_1 = OBS_TEMP[star - 1][clust_l - 1] / (long double)POS_TEMP[star - 1][clust_l - 1];
 //						if (param_starl_1 == 0 || param_starl_1 == 1){ starl_density_1 = starl_density_1; }
 //						else{
 //							starl_density_1 = starl_density_1 + (OBS_TEMP[star - 1][clust_l - 1])*log(param_starl_1) + (POS_TEMP[star - 1][clust_l - 1] - OBS_TEMP[star - 1][clust_l - 1])*log(1 - param_starl_1);
 //						}
 //					}
 //				}
 //				//cluster densities indexed by row l  (row l, col star)
 //				lstar_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_lstar_1 = OBS_TEMP[clust_l - 1][star - 1] / (long double)POS_TEMP[clust_l - 1][star - 1];
 //						if (param_lstar_1 == 0 || param_lstar_1 == 1){ lstar_density_1 = lstar_density_1; }
 //						else{
 //							lstar_density_1 = lstar_density_1 + (OBS_TEMP[clust_l - 1][star - 1])*log(param_lstar_1) + (POS_TEMP[clust_l - 1][star - 1] - OBS_TEMP[clust_l - 1][star - 1])*log(1 - param_lstar_1);
 //						}
 //					}
 //				}
 //				//cluster densities indexed by column j  (row star, col j)
 //				starj_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_starj_1 = OBS_TEMP[star - 1][clust_j - 1] / (long double)POS_TEMP[star - 1][clust_j - 1];
 //						if (param_starj_1 == 0 || param_starj_1 == 1){ starj_density_1 = starj_density_1; }
 //						else{
 //							starj_density_1 = starj_density_1 + (OBS_TEMP[star - 1][clust_j - 1])*log(param_starj_1) + (POS_TEMP[star - 1][clust_j - 1] - OBS_TEMP[star - 1][clust_j - 1])*log(1 - param_starj_1);
 //						}
 //					}
 //				}
 //				//cluster densities indexed by row l  (row j, col star)
 //				jstar_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_jstar_1 = OBS_TEMP[clust_j - 1][star - 1] / (long double)POS_TEMP[clust_j - 1][star - 1];
 //						if (param_jstar_1 == 0 || param_jstar_1 == 1){ jstar_density_1 = jstar_density_1; }
 //						else{
 //							jstar_density_1 = jstar_density_1 + (OBS_TEMP[clust_j - 1][star - 1])*log(param_jstar_1) + (POS_TEMP[clust_j - 1][star - 1] - OBS_TEMP[clust_j - 1][star - 1])*log(1 - param_jstar_1);
 //						}
 //					}
 //				}
 //				//Numerator of modified change-in-loglikelihood equation
 //				l_1 = starl_density_1 + lstar_density_1 + starj_density_1 + jstar_density_1
 //					+ l_density_1 + j_density_1 + lj_density_1 + jl_density_1;
 //			}
 //			//If weighted, all edges are non-negative integers and Poisson density is appropriate
 //			else
 //			{
 //				//Noe: Parameter estimates can - in theory - take on values of zero but poisson only takes parameter >0.
 //				//If a density's parameter estimate is 0, it will ruin the loglik calculation.
 //
 //				//cluster l density
 //				l_density_1 = 0;
 //				long double param_l_1 = OBS_TEMP[clust_l - 1][clust_l - 1] / (long double)POS_TEMP[clust_l - 1][clust_l - 1];
 //				if (param_l_1 == 0){ l_density_1 = 0; }
 //				else{
 //					l_density_1 = OBS_TEMP[clust_l - 1][clust_l - 1] * log(param_l_1) - param_l_1;
 //				}
 //
 //				//cluster j density
 //				j_density_1 = 0;
 //				long double param_j_1 = OBS_TEMP[clust_j - 1][clust_j - 1] / (long double)POS_TEMP[clust_j - 1][clust_j - 1];
 //				if (param_j_1 == 0){ j_density_1 = 0; }
 //				else{
 //					j_density_1 = OBS_TEMP[clust_j - 1][clust_j - 1] * log(param_j_1) - param_j_1;
 //				}
 //				//b/w cluster density indexed by l and j (row l, col j)
 //				lj_density_1 = 0;
 //				long double param_lj_1 = OBS_TEMP[clust_l - 1][clust_j - 1] / (long double)POS_TEMP[clust_l - 1][clust_j - 1];
 //				if (param_lj_1 == 0){ lj_density_1 = 0; }
 //				else{
 //					lj_density_1 = OBS_TEMP[clust_l - 1][clust_j - 1] * log(param_lj_1) - param_lj_1;
 //				}
 //				//b/w cluster density indexed by l and j (row j, col l)
 //				jl_density_1 = 0;
 //				long double param_jl_1 = OBS_TEMP[clust_j - 1][clust_l - 1] / (long double)POS_TEMP[clust_j - 1][clust_l - 1];
 //				if (param_jl_1 == 0){ jl_density_1 = 0; }
 //				else{
 //					jl_density_1 = OBS_TEMP[clust_j - 1][clust_l - 1] * log(param_jl_1) - param_jl_1;
 //				}
 //
 //				//cluster densities indexed by column l  (row star, col l)
 //				starl_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_starl_1 = OBS_TEMP[star - 1][clust_l - 1] / (long double)POS_TEMP[star - 1][clust_l - 1];
 //						if (param_starl_1 == 0){ starl_density_1 = starl_density_1; }
 //						else{
 //							starl_density_1 = starl_density_1 + OBS_TEMP[star - 1][clust_l - 1] * log(param_starl_1) - param_starl_1;
 //						}
 //					}
 //				}
 //				//cluster densities indexed by row l  (row l, col star)
 //				lstar_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_lstar_1 = OBS_TEMP[clust_l - 1][star - 1] / (long double)POS_TEMP[clust_l - 1][star - 1];
 //						if (param_lstar_1 == 0){ lstar_density_1 = lstar_density_1; }
 //						else{
 //							lstar_density_1 = lstar_density_1 + OBS_TEMP[clust_l - 1][star - 1] * log(param_lstar_1) - param_lstar_1;
 //						}
 //					}
 //				}
 //				//cluster densities indexed by column j  (row star, col j)
 //				starj_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_starj_1 = OBS_TEMP[star - 1][clust_j - 1] / (long double)POS_TEMP[star - 1][clust_j - 1];
 //						if (param_starj_1 == 0){ starj_density_1 = starj_density_1; }
 //						else{
 //							starj_density_1 = starj_density_1 + OBS_TEMP[star - 1][clust_j - 1] * log(param_starj_1) - param_starj_1;
 //						}
 //					}
 //				}
 //				//cluster densities indexed by row l  (row j, col star)
 //				jstar_density_1 = 0;
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						long double param_jstar_1 = OBS_TEMP[clust_j - 1][star - 1] / (long double)POS_TEMP[clust_j - 1][star - 1];
 //						if (param_jstar_1 == 0){ jstar_density_1 = jstar_density_1; }
 //						else{
 //							jstar_density_1 = jstar_density_1 + OBS_TEMP[clust_j - 1][star - 1] * log(param_jstar_1) - param_jstar_1;
 //						}
 //					}
 //				}
 //				//Numerator of modified change-in-loglikelihood equation
 //				l_1 = starl_density_1 + lstar_density_1 + starj_density_1 + jstar_density_1
 //					+ l_density_1 + j_density_1 + lj_density_1 + jl_density_1;
 //
 ////Just want to see if using a single density for ALL between-cluster edges improves performance
 ///*
 //				unsigned long long BW_OBS = 0;
 //				unsigned long long BW_POS = 0;
 //				for (int i = 0; i < k; i++)
 //				{
 //					for (int j = 0; j < k; j++)
 //					{
 //						if (i != j)
 //						{
 //							BW_OBS = BW_OBS + OBS_TEMP[i][j]; BW_POS = BW_POS + POS_TEMP[i][j];
 //						}
 //					}
 //				}
 //				l_1 = l_density_1 + j_density_1 + ((long double)BW_OBS)*log((long double)BW_OBS / BW_POS) + ((long double)BW_POS - BW_OBS)*log(1 - (long double)BW_OBS / BW_POS);
 //*/
 //
 //			}
 //			//Calculate change in log-likelihood due to reclustering
 //			//log(new/old)=log(new)-log(old)
 //
 //			//If all densities in either the new or the old loglik were discarded (and thus set to zero),
 //			//it is not possible to evaluate a change-in-loglik
 //			//set entire change-in-loglik to zero so a new clustering will be considered
 //			if (l_0 == 0 || l_1 == 0){ delta_l = 0; }
 //			else{ delta_l = l_1 - l_0; }
 //			//std::cout << "The change in log-lik is: " << delta_l << endl;	
 //			//if delta_l>0, then l_1>l_0, and the reclustered membership is more "likely"
 //			//Otherwise, l_1<l_0 and the current membership is more "likely"
 //
 //			overall_iteration_count++;
 ////long double prob_accept_modifier = (1 + log(epsilon / (long double)InitTemp) / (log((long double)CR))*TL) - (long double)overall_iteration_count -1;
 //long double prob_accept_modifier = 1;
 //long double uni_prob_modifier =1;
 //			int STAGE = 0; //Indicates the stage when a proposed recustering was accepted or rejected.
 //
 //			if (delta_l > epsilon)
 //			{//ACCEPT PROPOSED CLUSTERING AS CURRENT (best) CLUSTERING
 //			//	std::cout << "STAGE 1: " << delta_l << "..." << epsilon << endl;
 //				STG1_CT++;
 //				STAGE = 1;
 //				retain_current_sol = 0;
 //				//Setting this to 0 means:
 //				//1. We overwrite current membership with proposed membership
 //				//2. We overwrite current OBS/POS matrices with those of the proposed membership
 //				//3. We return to the top of the SA algorithm and calculate likelihood for the (new) current membership
 //
 //				//current clustering <- proposed clustering
 //				cluster_membership[l_star] = clust_j;
 //
 //				//Update OBS vector elements indexed only by l and j
 //				OBS[clust_l - 1][clust_l - 1] = OBS[clust_l - 1][clust_l - 1] - Z_to[clust_l - 1] - Z_from[clust_l - 1];
 //				OBS[clust_j - 1][clust_j - 1] = OBS[clust_j - 1][clust_j - 1] + Z_to[clust_j - 1] + Z_from[clust_j - 1];
 //				OBS[clust_l - 1][clust_j - 1] = OBS[clust_l - 1][clust_j - 1] - Z_to[clust_j - 1] + Z_from[clust_l - 1];
 //				OBS[clust_j - 1][clust_l - 1] = OBS[clust_j - 1][clust_l - 1] + Z_to[clust_l - 1] - Z_from[clust_j - 1];
 //				//Update POS vector elements indexed only by l and j
 //				//group_size[i-1] is the size of cluster i, n_i
 //				POS[clust_l - 1][clust_l - 1] = (group_size[clust_l - 1] - 1)*(group_size[clust_l - 1] - 2);
 //				POS[clust_j - 1][clust_j - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_j - 1]);
 //				POS[clust_l - 1][clust_j - 1] = (group_size[clust_l - 1] - 1)*(group_size[clust_j - 1] + 1);
 //				POS[clust_j - 1][clust_l - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_l - 1] - 1);
 //
 //				//Update OBS and POS vector elements indexed by either l OR j
 //				for (int star = 1; star <= k; star++)
 //				{
 //					if ((star != clust_l) && (star != clust_j))
 //					{
 //						OBS[star - 1][clust_l - 1] = OBS[star - 1][clust_l - 1] - Z_from[star - 1];
 //						OBS[clust_l - 1][star - 1] = OBS[clust_l - 1][star - 1] - Z_to[star - 1];
 //						OBS[star - 1][clust_j - 1] = OBS[star - 1][clust_j - 1] + Z_from[star - 1];
 //						OBS[clust_j - 1][star - 1] = OBS[clust_j - 1][star - 1] + Z_to[star - 1];
 //						//group_size[i] is the size of cluster i, n_i
 //						POS[star - 1][clust_l - 1] = (group_size[star - 1])*(group_size[clust_l - 1] - 1);
 //						POS[clust_l - 1][star - 1] = (group_size[clust_l - 1] - 1)*(group_size[star - 1]);
 //						POS[star - 1][clust_j - 1] = (group_size[star - 1])*(group_size[clust_j - 1] + 1);
 //						POS[clust_j - 1][star - 1] = (group_size[clust_j - 1] + 1)*(group_size[star - 1]);
 //					}
 //				}
 //				//Update group size vector
 //				group_size[clust_l - 1]--; //cluster l loses 1 vertex
 //				group_size[clust_j - 1]++; //cluster j gains 1 vertex
 //
 //				Success_counter++; //Iterate number of success at this temperature
 //
 //			}
 //			else
 //			{
 //				long double uni_draw = uniform_rand(rand());
 //				long double min_val = min((long double)1, exp(delta_l / (IT)));
 ////long double uni_draw = uni_prob_modifier*uniform_rand(rand());
 ////long double min_val = min((long double)1, exp(delta_l / (prob_accept_modifier*IT)));
 ////cout<<uni_draw<<"..."<<min_val<<"...";				
 //				//https://www.mathworks.com/help/gads/how-simulated-annealing-works.html?s_tid=gn_loc_drop
 //				//long double min_val =  min((long double)0.5, 1/(1+exp(-delta_l/IT)));
 //				//http://mathworld.wolfram.com/SimulatedAnnealing.html
 //				// double min_val = exp(-delta_l/IT);
 ////				if (min_val >= 1) { std::cout << "delta_L=" << delta_l<<", IT="<<IT <<", delta_l/IT="<<delta_l/IT<< ", exp(delta_l/IT)="<<exp(delta_l/IT)<<", 1/(1+exp(delta_l/IT))="<< 1 / (1 + exp(delta_l / IT))<<endl; std::cin.get(); }
 //				
 //				if (uni_draw < min_val) //Accept proposed clustering
 //				{
 //				//	std::cout << "STAGE2: " << uni_draw << "..." << min_val << endl;
 //				STG2_CT++;
 //					STAGE = 2;
 //					retain_current_sol = 0;
 //					//Setting this to 0 means:
 //					//1. We overwrite current membership with proposed membership
 //					//2. We overwrite current OBS/POS matrices with those of the proposed membership
 //					//3. We return to the top of the SA algorithm and calculate likelihood for the (new) current membership 
 //
 //					//current clustering <- proposed clustering
 //					cluster_membership[l_star] = clust_j;
 //					
 //					//Update OBS vector elements indexed only by l and j
 //					OBS[clust_l - 1][clust_l - 1] = OBS[clust_l - 1][clust_l - 1] - Z_to[clust_l - 1] - Z_from[clust_l - 1];
 //					OBS[clust_j - 1][clust_j - 1] = OBS[clust_j - 1][clust_j - 1] + Z_to[clust_j - 1] + Z_from[clust_j - 1];
 //					OBS[clust_l - 1][clust_j - 1] = OBS[clust_l - 1][clust_j - 1] - Z_to[clust_j - 1] + Z_from[clust_l - 1];
 //					OBS[clust_j - 1][clust_l - 1] = OBS[clust_j - 1][clust_l - 1] + Z_to[clust_l - 1] - Z_from[clust_j - 1];
 //					//Update POS vector elements indexed only by l and j
 //					//group_size[i-1] is the size of cluster i, n_i
 //					POS[clust_l - 1][clust_l - 1] = (group_size[clust_l - 1] - 1)*(group_size[clust_l - 1] - 2);
 //					POS[clust_j - 1][clust_j - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_j - 1]);
 //					POS[clust_l - 1][clust_j - 1] = (group_size[clust_l - 1] - 1)*(group_size[clust_j - 1] + 1);
 //					POS[clust_j - 1][clust_l - 1] = (group_size[clust_j - 1] + 1)*(group_size[clust_l - 1] - 1);
 //
 //					//Update OBS and POS vector elements indexed by either l OR j
 //					for (int star = 1; star <= k; star++)
 //					{
 //						if ((star != clust_l) && (star != clust_j))
 //						{
 //							OBS[star - 1][clust_l - 1] = OBS[star - 1][clust_l - 1] - Z_from[star - 1];
 //							OBS[clust_l - 1][star - 1] = OBS[clust_l - 1][star - 1] - Z_to[star - 1];
 //							OBS[star - 1][clust_j - 1] = OBS[star - 1][clust_j - 1] + Z_from[star - 1];
 //							OBS[clust_j - 1][star - 1] = OBS[clust_j - 1][star - 1] + Z_to[star - 1];
 //							//group_size[i] is the size of cluster i, n_i
 //							POS[star - 1][clust_l - 1] = (group_size[star - 1])*(group_size[clust_l - 1] - 1);
 //							POS[clust_l - 1][star - 1] = (group_size[clust_l - 1] - 1)*(group_size[star - 1]);
 //							POS[star - 1][clust_j - 1] = (group_size[star - 1])*(group_size[clust_j - 1] + 1);
 //							POS[clust_j - 1][star - 1] = (group_size[clust_j - 1] + 1)*(group_size[star - 1]);
 //						}
 //					}
 //
 //					//Update group size vector
 //					group_size[clust_l - 1]--; //cluster l loses 1 vertex
 //					group_size[clust_j - 1]++; //cluster j gains 1 vertex						
 //
 //					Success_counter++; //Iterate number of success at this temperature
 //				}
 //				else //Reject proposed clustering
 //				{//Nothing happens here. The current cluster is retained and the process repeats itself.	
 //				//	std::cout << "STAGE3: Rejected..." << endl;
 //				STG3_CT++;
 //				
 //					STAGE = 3;
 //					retain_current_sol = 1;
 //					//Setting this to 1 means:
 //					//1. We retain the current membership
 //					//2. We retain the current OBS/POS matrices
 //					//3. We retain the loglik associated with the current membership
 ////std::cout << "min_val="<<min_val<<", uni_draw="<<uni_draw<<", delta_L=" << delta_l << ", IT=" << IT << ", delta_l/IT=" << delta_l / IT << ", exp(delta_l/IT)=" << exp(delta_l / IT) << endl;
 //				}
 //
 //			}
 ////long double sol_loglik_temp=0;//loglikelihood for the solution at the current value of k
 ////long double sol_OBS_temp = 0;//Total number of observed edges
 ////long double sol_POS_temp = 0;//Total number of possible edges
 ////long double LO_temp = 0;		//loglikelihood under the null hypothesis: k=1
 ////	for (int row = 0; row < k; row++){
 ////		for (int col = 0; col < k; col++){
 ////			if (OBS[row][col] == 0 || OBS[row][col] == POS[row][col]){ sol_loglik_temp = sol_loglik_temp; }
 ////			else{
 ////				sol_loglik_temp = sol_loglik_temp + (OBS[row][col])*log(OBS[row][col] / (long double)POS[row][col]) + (POS[row][col] - OBS[row][col])*log(1 - OBS[row][col] / (long double)POS[row][col]);
 ////			}
 ////			sol_OBS_temp = sol_OBS_temp + OBS[row][col];
 ////		}
 ////	}
 ////	sol_POS_temp = (long double)N*(N - 1); //2*[Nchoose2]=N(N-1) {or N*2-N, b/c no self-edges=NN-N=N(N-1)}
 ////	LO_temp = (sol_OBS_temp)*log(sol_OBS_temp / sol_POS_temp) + (sol_POS_temp - sol_OBS_temp)*log(1 - sol_OBS_temp / sol_POS_temp);
 ////			//At this point, the proposed clustering has:
 ////			//1. Been Accepted. membership vector/OBS/POS/cluster counts have beeen updated to reflect new clustering OR
 ////			//2. Rejected. memberhsip vector/OBS/POS/cluster counts still correspond to the current clustering.
 ////std::cout << STG1_CT << "..." << STG2_CT << "..." << STG3_CT << "..."<< overall_iteration_count<<"..."<<prob_accept_modifier<<"..."<<STAGE<<"..."<<sol_loglik_temp <<"..."<<LO_temp <<endl;
 //		
 //}//End SA LOOP
 //
 //
 //		//if unweighted, all edges have weight=1 and Binomial density is appropriate
 //		long double sol_loglik=0;//loglikelihood for the solution at the current value of k
 //		long double sol_OBS = 0;//Total number of observed edges
 //		long double sol_POS = 0;//Total number of possible edges
 //		long double LO = 0;		//loglikelihood under the null hypothesis: k=1
 //		if (weighted == 'N')
 //		{
 //			for (int row = 0; row < k; row++){
 //				for (int col = 0; col < k; col++){
 //					if (OBS[row][col] == 0 || OBS[row][col] == POS[row][col]){ sol_loglik = sol_loglik; }
 //					else{
 //						sol_loglik = sol_loglik + (OBS[row][col])*log(OBS[row][col] / (long double)POS[row][col]) + (POS[row][col] - OBS[row][col])*log(1 - OBS[row][col] / (long double)POS[row][col]);
 //					}
 //					sol_OBS = sol_OBS + OBS[row][col];
 //				}
 //			}
 //			sol_POS = (long double)N*(N - 1); //2*[Nchoose2]=N(N-1) {or N*2-N, b/c no self-edges=NN-N=N(N-1)}
 //			LO = (sol_OBS)*log(sol_OBS / sol_POS) + (sol_POS - sol_OBS)*log(1 - sol_OBS / sol_POS);
 //		}
 //		//otherwise, Poisson
 //		else{
 //			for (int row = 0; row < k; row++){
 //				for (int col = 0; col < k; col++){
 //					if (OBS[row][col] == 0 ){ sol_loglik = sol_loglik; }
 //					else{
 //						sol_loglik = sol_loglik + OBS[row][col] * log(OBS[row][col] / (long double)POS[row][col]) - OBS[row][col] / (long double)POS[row][col];
 //					}
 //					sol_OBS = sol_OBS + OBS[row][col];
 //				}
 //			}
 //			sol_POS = (long double)N*(N - 1); //2*[Nchoose2]=N(N-1) {or N*2-N, b/c no self-edges=NN-N=N(N-1)}
 //			LO = sol_OBS * log(sol_OBS / sol_POS) - sol_OBS / sol_POS;
 //		}
 //
 //
 //		//MEMBERSHIP.push_back(cluster_membership);
 //			//NOTE: You can uncomment this line to store the solution membership vector for EACH value of k in [min_k, max_k].
 //			//However, this will consume a lot of memory for large networks and large ranges of k. Consider N=100,000 and k in [1000,10000]
 //			//MEMBERSHIP would then me a 100,000x9,000 element matrix.
 //			//Current code only stores the solution membership vector for the k value with minimum BIC
 //		LOGLIK.push_back(sol_loglik);
 //		TESTSTAT.push_back(-2 * (LO - sol_loglik)); 	// -2 * (LO - full_loglik_curr)
 //
 //		//Calculate Bayesian Information Criterion for solution clustering
 //		//https://en.wikipedia.org/wiki/Bayesian_information_criterion
 //		//BIC=-2*loglik+(log n)*z
 //		//BIC=ln(n)z-2ln(L), where 
 //		//n=# of obs=# of possible edges =2*[Nchoose2] =N(N-1)
 //			//Or equivalently, N^2 elements in adjacency matrix - N unused diagonal elements = NN - N = N(N-1)
 //		//z=# of model parameters = 2*kchoose2 between-cluster parameters + k within-cluster parameters 
 //			//2*(k)(k-1)/2 + k = k(k-1)+k
 //		//ln(L)=maximized loglikelihood of model = sol_loglik, which is already the maximized log-likelihood
 //		BIC.push_back(log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik);
 //
 //		//Store solutions for first value of k (or ONLY value of k if min_k=max_k)
 //		if (k == min_k) {
 //			final_BIC = log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik;
 //			final_k = k;
 //			final_LOGLIK = sol_loglik;
 //			final_TESTSTAT = -2 * (LO - sol_loglik);
 //			final_membership = cluster_membership;
 //		}
 //		//If current BIC is less than all previous BICs, store current solution as the final solution (thus far)
 //		if (log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik < final_BIC) {
 //			final_k = k;
 //			final_BIC = log(N*(N - 1))*(k*(k - 1) + k) - 2 * sol_loglik;
 //			final_LOGLIK = sol_loglik;
 //			final_TESTSTAT = -2 * (LO - sol_loglik);
 //			final_membership=cluster_membership;
 //		}
 //	}//end for-k loop
 //
 //	time_t SA_end = time(0);
 //	double SA_time = difftime(SA_end, SA_start);
 //
 //	std::cout << fixed;
 //	std::cout << setprecision(10);
 //	 //Output Solutions
 //	 std::cout << "******************************************************************" << endl;
 //	 std::cout << "Execution time = " << SA_time << " seconds" << "\n"; /*time / 60*/	/*minutes*/
 //	 std::cout << "******************************************************************" << endl << endl;
 //
 //	std::cout << "***********k values******" << endl;
 //	for (int k = min_k; k <= max_k; k++) { std::cout << k << " "; }
 //	std::cout << endl;
 //	std::cout << "***********BIC***********" << endl;
 //	for (int i = 0; i < BIC.size(); i++) { std::cout << BIC[i] << " "; }
 //	std::cout << endl;
 //	std::cout << "***********LOGLIK***********" << endl;
 //	for (int i = 0; i < LOGLIK.size(); i++) { std::cout << LOGLIK[i] << " "; }
 //	std::cout << endl;
 //	std::cout << "***********TESTSTAT***********" << endl;
 //	for (int i = 0; i < TESTSTAT.size(); i++) { std::cout << TESTSTAT[i] << " "; }
 //	std::cout << endl << endl;
 //
 //	std::cout << "********************************" << endl;
 //	std::cout << "********************************" << endl;
 //	std::cout << "Value for solution cluster:" << endl;
 //	std::cout << "k = " << final_k << endl;
 //	std::cout << "BIC = " << final_BIC << endl;
 //	std::cout << "loglik = " << final_LOGLIK << endl;
 //	std::cout << "test statistic = " << final_TESTSTAT << endl;
 //
 //	ostringstream bestzsol;
 //	ofstream bestzstorage;
 //	// If old file exists from a previous run, delete it.
 //	remove("SOLUTION_CLUSTERING.txt");
 //	bestzsol << "SOLUTION_CLUSTERING" << ".txt";
 //	bestzstorage.open("SOLUTION_CLUSTERING.txt", ios::app);
 //	for (int i = 0; i < final_membership.size(); i++) { bestzstorage << final_membership[i]; } // << "\t";}
 //	bestzstorage << endl;
 //	bestzstorage.close();
 //	std::cout << "The membership vector for this clustering has been written to SOLUTION_CLUSTERING.txt." << endl;
 //	