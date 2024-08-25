// basic
#include <iostream>
#include <NTL/RR.h>
#include <cmath>
#include "PolyUpdate.h"
#include "program.h"
#include "seal/seal.h"
#include "SEALcomp.h"

using namespace seal;
using namespace std;
using namespace NTL;


int main() {

	// SEAL User setting
	long level = 14;		// number of levels L. L = sum ceil(log2(di)) for comparison operation and L = sum ceil(log2(di))+1 for max & ReLU function. di is degree of component polynomials.
	long alpha = 13;			// precision parameter alpha
	long comp_no = 3;		// number of compositions
	long scalingfactor = 45;		// log2 of scaling factor
	vector<int> deg = {15,15,27};		// degrees of component polynomials
	double eta = pow(2.0,-15);		// margin
	double scaled_val = 1.7;		// scaled_val: the last scaled value
	double max_factor = 16;		// max_factor = 1 for comparison operation. max_factor > 1 for max or ReLU function
	vector<Tree> tree;		// structure of polynomial evaluation
	evaltype eval_type = evaltype::oddbaby;
	RR::SetOutputPrecision(25);

	// generate tree
	for(int i=0; i<comp_no; i++) 
	{
		Tree tr;
		if(eval_type == evaltype::oddbaby) upgrade_oddbaby(deg[i], tr);
		else if(eval_type == evaltype::baby) upgrade_baby(deg[i], tr);
		else std::invalid_argument("evaluation type is not correct");
		tree.emplace_back(tr);
		tr.print();
	}
    testReLU(59);



}