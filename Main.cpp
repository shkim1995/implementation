#include "Data.h"
#include "Histogram.h"
#include "Tree.h"
#include <limits>
#include <cstddef>

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


#define ENTROPY_PRUNING
#define DIST_PRUNING


using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////GLOBAL VARIABLES//////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

Node tree;

//Shapelet info
vector <double> Shapelet;
//cutting distance
double Shapelet_dist;
//what is left of cutting distance (0/1)
int Shapelet_left;


int MAXLEN = 60;
int MINLEN = 30;

// int a_total=0;
// int b_total=0;
// int N=0;

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////TRAINING  FUNCTIONS/////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


//FOR DEBUGGINH
void printVector(vector<double> S){
	for(int i=0; i<S.size(); i++)
		cout<<S[i]<<" ";
	cout<<endl;
}


//get Distance of sequence T and subsequence S from starting index l
double SubsequenceDist(vector<double> T, vector<double> S, int l, double min_dist){
	if(T.size()<S.size()+l)
		return -1;
	double dist = 0;
	for(int i=0; i<S.size(); i++){
		double S_val = S[i];
		double T_val = T[i+l];
		dist += (S_val-T_val)*(S_val-T_val);

#ifdef DIST_PRUNING
		if(min_dist>0 && dist>min_dist){
			dist = numeric_limits<double>::max();
			break;
		}
#endif

	}

	return dist;
}

//get minimum Distance of sequence T and subsequence S
double SubsequenceDist(vector<double> T, vector<double> S){
	int l = 0;
	double min_Dist = -1;
	while(true){
		double dist = SubsequenceDist(T, S, l, min_Dist);
		if(dist<0)
			break;
		if(min_Dist<0)
			min_Dist = dist;
		else{
			if(min_Dist>dist)
				min_Dist = dist;
		}
		l++;
	}
	return min_Dist;
}

//Entropy before division
// a, b : number of each elements 
double Entropy(int  a, int b){
	if(a==0) return 0;
	if(b==0) return 0;
	double pa = (a+0.0)/(a+b+0.0);
	double pb = (b+0.0)/(a+b+0.0);
	return -(pa*log(pa)+pb*log(pb));
}


//Entropy after division
double Entropy(int a1, int b1, int a2, int b2){
	return ((a1+b1)*Entropy(a1, b1)+(a2+b2)*Entropy(a2, b2))/(a1+a2+b1+b2);
}


double CalculateInfoGain(Histogram* hist, int a_total, int b_total){
	list<HistData>::iterator it = hist->data.begin();

	int a1 = 0;
	int b1 = 0;
	int a2 = a_total;
	int b2 = b_total;

	it = hist->data.begin();
	int cut_pos = 0;
	int pos = 0;
	double minEntropy = Entropy(a1, b1, a2, b2);
	double totalEntropy = minEntropy;
	while(it!=hist->data.end()){

		double entropy = Entropy(a1, b1, a2, b2);

		// cout<<entropy<<" for ";
		// cout<<"("<<a1<<", "<<b1<<"), ("<<a2<<", "<<b2<<")"<<endl;

		//update if current entropy is minimum
		if(entropy<minEntropy){
			cut_pos = pos;
			minEntropy = entropy;
		}

		if(it->type==0){
			a1++;
			a2--;
		}
		else{
			b1++;
			b2--;
		}
		it++;
		pos++;
	}

	//cut between cut_pos and cut_pos-1
	return totalEntropy-minEntropy;
	
}



bool HistogramPruning(Histogram* hist, double max_gain, int a_total, int b_total){
	list<HistData>::iterator it = hist->data.begin();

	int min_d = 0;
	int max_d = 100000;

	int a_num = 0;  //count for 0
	int b_num = 0;  //count for 1
 
	while(it!=hist->data.end()){
		if(it->type==0)
			a_num++;
		else
			b_num++;
		it++;
	}

	//needed 0 and 1 values
	a_num = a_total-a_num;
	b_num = b_total-b_num;

	
	it = hist->data.begin();

	Histogram h1;
	Histogram h2;

	while(it!=hist->data.end()){
		
		h1.insert(it->type, it->value);
		h2.insert(it->type, it->value);
		it++;
	}
	
	

	for(int i=0; i<a_num; i++){
		h1.insert(0, min_d);
		h2.insert(0, max_d);
	}


	for(int i=0; i<b_num; i++){
		h1.insert(1, max_d);
		h2.insert(1, min_d);
	}
		

	if(CalculateInfoGain(&h1, a_total, b_total)>max_gain){

		//cout<<CalculateInfoGain(&h1)<<" "<<max_gain<<endl;
		return false;
	}

	if(CalculateInfoGain(&h2, a_total, b_total)>max_gain){
		//cout<<CalculateInfoGain(&h2)<<" "<<max_gain<<endl;
		return false;
	}

	return true;

}

double CheckCandidate(vector<double> S, double max_gain, vector<Data>* dataVector, int a_total, int b_total){
	//cout<<max_gain;
	Histogram hist;
	for(int i=0; i<dataVector->size(); i++){
		double distance = SubsequenceDist((*dataVector)[i].data, S);
		hist.insert((*dataVector)[i].type, distance);

#ifdef ENTROPY_PRUNING
		if(HistogramPruning(&hist, max_gain, a_total, b_total)==true)
		  	return -1;
#endif

	}
	return CalculateInfoGain(&hist, a_total, b_total);
}


//Set node value
double SetShapelet(int k, int i, int j, vector<Data>* dataVector, int a_total, int b_total, Node* node){
	Data temp = (*dataVector)[k];

	for(int p=i; p<=j; p++){
		node->shapelet.push_back(temp.data[p]);
	}

	Histogram hist;
	for(int i=0; i<dataVector->size(); i++){
		double distance = SubsequenceDist((*dataVector)[i].data, node->shapelet);
		hist.insert((*dataVector)[i].type, distance);
	}

	list<HistData>::iterator it = hist.data.begin();

	int a1 = 0;
	int b1 = 0;
	int a2 = a_total;
	int b2 = b_total;

	it = hist.data.begin();
	int cut_pos = 0;
	int pos = 0;
	double minEntropy = Entropy(a1, b1, a2, b2);
	double totalEntropy = minEntropy;

	while(it!=hist.data.end()){

		double entropy = Entropy(a1, b1, a2, b2);

		// cout<<entropy<<" for ";
		// cout<<"("<<a1<<", "<<b1<<"), ("<<a2<<", "<<b2<<")"<<endl;

		//update if current entropy is minimum
		if(entropy<minEntropy){
			cut_pos = pos;
			minEntropy = entropy;
		}

		if(it->type==0){
			a1++;
			a2--;
		}
		else{
			b1++;
			b2--;
		}
		it++;
		pos++;
	}
	
	//SET Shapelet_left
	int count=0;
	int cur = 0;
	it = hist.data.begin();
	while(cur<cut_pos){
		if(it->type==0)
			count++;
		cur++;
		it++;
	}
	if(count>cut_pos/2){
		node->left_type = 0;
		node->right_type = 1;
	}
	else{
		node->left_type = 1;
		node->right_type = 0;
	}


	//SET Shapelet_dist
	double dist1 = it->value;
	it++;
	double dist2 = it->value;
	node->dist = (dist1+dist2)/2;

	//cout<<Shapelet_dist<<endl;	

	return totalEntropy-minEntropy;
}

bool isDataVectorPure(vector<Data> S, int type){
	for(int i=0; i<S.size(); i++)
		if(S[i].type!=type)
			return false;
	return true;
}

void PrintDataType(vector <Data> S){
	for(int i=0; i<S.size(); i++)
		cout<<S[i].type<<" ";
	cout<<endl;
}

void FindingShapeleteBF(vector<Data>* dataVector, Node* node){

	int a_total = 0;
	int b_total = 0;

	for(int k=0; k<dataVector->size(); k++){
		if((*dataVector)[k].type==0)
			a_total++;
		else
			b_total++;
	}


	double gain = 0;
	int shapelet_k;
	int shapelet_i;
	int shapelet_j;



	//get candidate
	for(int k=0; k<dataVector->size(); k++){
		Data d = (*dataVector)[k];

		for(int i=0; i<d.size(); i++){
			for(int j=i+MINLEN; j<fmin(i+MAXLEN, d.size()); j++){
				
				//candidate
				vector<double> candidate;
				for(int t=i; t<=j; t++)
					candidate.push_back(d.data[t]);
				//candidate_done
				

				double temp = CheckCandidate(candidate, gain, dataVector, a_total, b_total);
				//cout<<"("<<k<<", "<<i<<", "<<j<<") : gain = "<<temp<<", min = "<<gain<<endl;
				if(temp>gain){
					gain = temp;
					shapelet_k = k;
					shapelet_i = i;
					shapelet_j = j;
				}

			}
		}
		cout<<"Round "<<k<<" : "<<shapelet_k<<" "<<shapelet_i<<" "<<shapelet_j<<" gain = "<<gain<<endl;
				
	}
	cout<<"SHAPLET is " <<shapelet_k<<" "<<shapelet_i<<" "<<shapelet_j<<endl;

	SetShapelet(shapelet_k, shapelet_i, shapelet_j, dataVector, a_total, b_total, node);


	//FOR Recursion
	vector <Data> dataVector_left;
	vector <Data> dataVector_right;

	for(int k=0; k<dataVector->size(); k++){

		double distance = SubsequenceDist((*dataVector)[k].data, node->shapelet);
		if(distance<node->dist)
			dataVector_left.push_back((*dataVector)[k]);
		else
			dataVector_right.push_back((*dataVector)[k]);
	}

	// cout<<node->left_type<<endl;
	// cout<<node->right_type<<endl;
	// cout<<endl;

	// cout<<isDataVectorPure(dataVector_left, node->left_type)<<endl;
	// cout<<isDataVectorPure(dataVector_right, node->right_type)<<endl;



	cout<<"_________________________________________________"<<endl;
	cout<<endl;

	if(isDataVectorPure(dataVector_left, node->left_type)==false){
		node->left = new Node();
		FindingShapeleteBF(&dataVector_left, node->left);
	}

	if(isDataVectorPure(dataVector_right, node->right_type)==false){
		node->right = new Node();
		FindingShapeleteBF(&dataVector_right, node->right);
	}


}

//Normalize dataVector(Training Set)
void Normalize(vector<Data>* dataVector){
	for(int k=0; k<dataVector->size(); k++){
		Data* d = &(*dataVector)[k];
		double minValue = d->data[0];
		for(int i=0; i<d->size(); i++){
			if(d->data[i]<minValue)
				minValue = d->data[i];
		}
		for(int i=0; i<d->size(); i++){
			d->data[i] = d->data[i]- minValue;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////TESTING FUNCTIONS//////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

void SetType(Data* data, Node* node){
	double dist = SubsequenceDist(data->data, node->shapelet);
	
	if(dist<node->dist){
		if(node->left==NULL){
			data->type2 = node->left_type;
			return;
		}
		else{
			SetType(data, node->left);
			return;
		}
		
	}
	else{
		if(node->right==NULL){	
			data->type2 = node->right_type;
			return;
		}
		else{
			SetType(data, node->right);
			return;
		}
		
	}	
}

//Classify test set using Shapelet, Shapelet_dist and Shapelet_left
void ClassifyTestSet(vector <Data>* dataVector2, Node* node){

	for(int k=0; k<dataVector2->size(); k++){
	
		SetType(&(*dataVector2)[k], node);
	}
}

double GetAccuracy(vector <Data>* dataVector2){
	int total = 0;
	int correct = 0;
	for(int k=0; k<dataVector2->size(); k++){
		total++;
		if((*dataVector2)[k].type==(*dataVector2)[k].type2)
			correct++;
	}
	return correct*100.0/total;
}



///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////MAIN  FUNCTION///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


int main(){

	clock_t before;

	before = clock();

	vector <Data> dataVector;
 	ifstream inFile;
	inFile.open("gun_train_2", ios::in);
	string s;

	while(getline(inFile, s)){
		stringstream ss(s);
		double test;
		int t;
		Data temp;

		ss>>t;
		temp.type = t;
		
		while(ss>>test){
			temp.data.push_back(test);
		}
		dataVector.push_back(temp);
	}

	

	Normalize(&dataVector);

	FindingShapeleteBF(&dataVector, &tree);


	//GET TIME
	double result = (double)(clock()-before)/CLOCKS_PER_SEC;
	cout<<"Time = "<<result<<endl;


	// cout<<endl;
	// cout<<endl;
	// cout<<"_________________________________________________"<<endl;
	// cout<<"______________________TREE_______________________"<<endl;
	
	ofstream outFile("output_no_pruning.txt");
	tree.print_nodes(&outFile);

	// cout<<"_________________________________________________"<<endl;
	// cout<<"_________________________________________________"<<endl;
	// cout<<endl;
	// cout<<endl;



	/////////TESTING//////////
	
	vector <Data> dataVector2;	

	ifstream inFile2;
	inFile2.open("gun_test", ios::in);
	string s2;

	while(getline(inFile2, s2)){
		stringstream ss(s2);
		double test;
		int t;
		Data temp;

		ss>>t;
		temp.type = t;

		while(ss>>test){
			temp.data.push_back(test);
		}
		dataVector2.push_back(temp);
	}


	Normalize(&dataVector2);

	ClassifyTestSet(&dataVector2, &tree);
	cout<<"The Accuracy of this Algorithm is : "<<GetAccuracy(&dataVector2)<<"%"<<endl;

	//for Debugging//


	// vector<double> T;
	// vector<double> S;

	// T.push_back(2);
	// T.push_back(2);
	// T.push_back(3);
	// T.push_back(2);
	// T.push_back(3);
	// T.push_back(2);
	// T.push_back(1);
	// T.push_back(2);
	// T.push_back(2);
	// T.push_back(1);
	// T.push_back(5);
	// T.push_back(8);
	// T.push_back(7);
	// T.push_back(2);
	// T.push_back(3);
	// T.push_back(2);
	// T.push_back(3);
	// T.push_back(1);
	// T.push_back(3);
	// T.push_back(1);
	// T.push_back(2);

	// S.push_back(5.2);
	// S.push_back(7.9);
	// S.push_back(7.1);

	// cout<<SubsequenceDist(T, S);


	// Data temp1;
	// temp1.data.push_back(5);
	// temp1.data.push_back(4);
	// temp1.data.push_back(3);
	// temp1.data.push_back(2);
	// temp1.data.push_back(1);
	// temp1.data.push_back(0);
	// temp1.data.push_back(-1);
	// temp1.data.push_back(-2);
	// temp1.data.push_back(-1);
	// temp1.data.push_back(1);
	// temp1.type = 0;


	// Data temp2;
	// temp2.data.push_back(5);
	// temp2.data.push_back(4.2);
	// temp2.data.push_back(3.2);
	// temp2.data.push_back(1.8);
	// temp2.data.push_back(1);
	// temp2.data.push_back(0);
	// temp2.data.push_back(1);
	// temp2.data.push_back(0);
	// temp2.data.push_back(0);
	// temp2.data.push_back(1);
	// temp2.type = 1;


	// Data temp3;
	// temp3.data.push_back(1);
	// temp3.data.push_back(1);
	// temp3.data.push_back(1);
	// temp3.data.push_back(2);
	// temp3.data.push_back(1);
	// temp3.data.push_back(0);
	// temp3.data.push_back(-1);
	// temp3.data.push_back(-2.5);
	// temp3.data.push_back(-1);
	// temp3.data.push_back(1);
	// temp3.type = 0;


	// Data temp4;
	// temp4.data.push_back(1);
	// temp4.data.push_back(1);
	// temp4.data.push_back(1);
	// temp4.data.push_back(2);
	// temp4.data.push_back(1);
	// temp4.data.push_back(0);
	// temp4.data.push_back(1);
	// temp4.data.push_back(-1);
	// temp4.data.push_back(0);
	// temp4.data.push_back(1);
	// temp4.type = 1;

	// dataVector.push_back(temp1);
	// dataVector.push_back(temp2);
	// dataVector.push_back(temp3);
	// dataVector.push_back(temp4);
	

	// Histogram hist;
	// hist.insert(1, 10);
	// hist.insert(0, 11);
	// hist.insert(0, 4);
	// hist.insert(1, 13);
	// hist.insert(1, 15);
	// hist.insert(0, 3);
	// hist.insert(1, 17);
	// hist.insert(0, 4);

	// hist.print();

	// double pos = CalculateInfoGain(&hist);

	// cout<<endl<<pos<<endl;
	//////////

}
