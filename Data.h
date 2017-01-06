
#ifndef HEADERS

#define HEADERS
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <sstream>
#include "time.h"

#endif


using namespace std;


// Class for each element of given time series
// includes the value and the type 
class Data{

	public:
		int type;
		vector<double> data;
		int type2;

		void print(){
			cout<<type<<" | ";
			for(int i=0; i<data.size(); i++)
				cout<<" "<<data[i];
			cout<<endl;
			cout<<endl;
		}

		int size(){
			return data.size();
		}


};