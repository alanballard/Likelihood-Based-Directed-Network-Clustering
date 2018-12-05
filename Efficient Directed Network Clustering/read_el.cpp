#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

#include <tuple>
//Can use tuple<int, int, int> for making triple ints.
//#include <tuple> 
//tuple<int, int, int> my_name (3,5,7)
//get<0>(my_name) //=3
//get<1>(my_name) //=5
//get<2>(my_name) //=7

using namespace std;
//Reads in edge list from text and converts to adjacency list.
//This code assumes a directed edge list.
//Input  text shoud be set up as: vertex 1 vertex 2, if unweighted
//or vertex 1 vertex 2 WEIGHT, if weighted.
//Since this is directed vertex 1 vertex 2 (WEIGHT) implies an edge extends FROM vertex 1 TO vertex 2 (with weight=WEIGHT)
//Vertex numbering must start at 0 and cannot skip numbers

//Directed edge list will be read into an adjacency list that will have two entries for each directed edge
//Assuming there is an edge FROM vertex 1 TO vertex 2 with weight=WEIGHT, there will be two entries in the adjacency list:
//vertex 1's row: a 3-tuple that contains the int 2, the int WEIGHT and an indicator (T) that indicates TO
//vertex 2's row: a 3-type that contains the int 1, the int WEIGHT and an indicator (F) that indicates FROM
//vector< vector<pair<int, int>> > read_el(char *filename, char weighted)
vector< vector< tuple<int, int, char>    >> read_el(char *filename, char weighted)
{
	//cout << "Inputting text graph..." << endl;
	ifstream finput;
	finput.open(filename, fstream::in);
	

	int nb_links = 0;
		
	vector< vector< tuple<int, int, char> >> links;


//**********************
//**********************
//Look at the file Graph.CPP on the thumb drive, C++ Code, Generic Louvain, src folder
//That shows how to convert this to a weighted network
//Need to change input format to test for presence of weights, or ask if a weighted network at the start
//If so, input format should be "src	dest	weight"
//If not, then "src	dest" and the weight should be filled in as one.
	while (!finput.eof()) {
		 int src, dest;		
		 int weight;

		 //string temp;
		 //getline(finput, temp);
		 //if (temp == ""){ cout << endl << "YOYOYO"; cin.get(); cin.get(); }// continue;
		 // if (temp.length()==0/*!temp.empty()*/){ cout << endl << "YOYOYO"; cin.get(); cin.get(); }

		 if (weighted == 'Y') //weighted case
			{
				finput >> src >> dest >> weight;

				if (finput.eof()) break;//Stops input if last line(s) of edge list are blank

				if (finput) {
					//make room for each new src/dest combination ...
					if (links.size() <= max(src, dest) + 1) //if the size of the vector links is less than (the largest id# of the nodes+1)...
					{
						links.resize(max(src, dest) + 1); //resize the vector links to size (largest id# of the nodes+1) elements
					}
					
					//This is an entry in the src-th position indicating the link FROM src TO dest with weight=weight
					//This is an out-link for src
					links[src].push_back(make_tuple(dest,weight,'T')); // 2<3,1,'T'> means a link from vertex 2 TO vertex 3 with weight=1

					nb_links++; //increase number of links by one.
				}

				//This is an entry in the dest-th position indicating the link FROM src TO dest with weight=weight
				//This is an in-link for src
				links[dest].push_back(make_tuple(src, weight, 'F')); // 2<3,1,'F'> means a link to vertex 2 FROM vertex 3 with weight=1
			}
			else//unweighted case
			{
					//read the line
					finput >> src >> dest;

					if (finput.eof()) break;//Stops input if last line(s) of edge list are blank

					if (finput) {
						//make room for each new src/dest combination ...
						if (links.size() <= max(src, dest) + 1) //if the size of the vector links is less than (the largest id# of the nodes+1)...
						{
							links.resize(max(src, dest) + 1); //resize the vector links to size (largest id# of the nodes+1) elements
						}

						//This is an entry in the src-th position indicating the link FROM src TO dest with weight=weight
						//This is an out-link for src
						links[src].push_back(make_tuple(dest, 1, 'T')); // 2<3,1,'T'> means a link from vertex 2 TO vertex 3 with weight=1

						nb_links++; //increase number of links by one.
					}

					//This is an entry in the dest-th position indicating the link FROM src TO dest with weight=weight
					//This is an in-link for dest
					links[dest].push_back(make_tuple(src, 1, 'F')); // 2<3,1,'F'> means a link to vertex 2 FROM vertex 3 with weight=1

			}

	}

//	cout << "***The number of distinct vertices is: " << links.size() << "***" << endl;
	cout << "Distinct edge count: " << nb_links << endl;

	finput.close();

	return links;
	//cout << "Press Enter to Continue." << endl; cin.get();
}
