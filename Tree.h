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


class Node{

public:

	vector<double> shapelet ;
	double dist;
	int left_type;
	int right_type;

	Node* left;
	Node* right;

	Node(){
		left_type = -1;
		right_type = -1;

		left = NULL;
		right = NULL;

		dist = 0;
	}

	void print_nodes(){
		cout<<dist<<endl;
		//printVector(shapelet);
		if(left==NULL)
			cout<<"null"<<endl;
		else
			left->print_nodes();

		if(right==NULL)
			cout<<"null"<<endl;
		else
			right->print_nodes();
	}

private:

	void printVector(vector<double> S){
	for(int i=0; i<S.size(); i++)
		cout<<S[i]<<" ";
	cout<<endl;
	}

};

// class Tree{

// public:

// 	Node root;

// 	Tree(){
// 		root.value = 0;
// 		root.left_type = 1;
// 		root.left = null;
// 		root.right = null;
// 	}

// };