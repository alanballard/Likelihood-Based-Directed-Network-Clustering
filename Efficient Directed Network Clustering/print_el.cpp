#include <iostream>
#include <iomanip>

#include <vector>
#include <algorithm>

#include <tuple>

using namespace std;

void print_el(vector<vector< tuple<int, int, char> >> links)
{
	
	//get<0>(links) -Destination vertex
	//get<1>(links) -Weight
	//get<2>(links) -char in {T,F}
		cout << "As an Adjacency List..." << endl;
		for (unsigned int i = 0; i < links.size(); i++) {
			cout << "Vertex " << i << ":";
			for (unsigned int j = 0; j < links[i].size(); j++) {
				cout << " " << get<0>(links[i][j]) << " <" << get<1>(links[i][j]) << ", " << get<2>(links[i][j]) <<  ">";
				//cout << " " <<  links[i][j].first << " <"<< links[i][j].second <<">";
			}
			cout << endl;
		}
	
	//cout << "Press Enter to Continue..." << endl; cin.get();
}
