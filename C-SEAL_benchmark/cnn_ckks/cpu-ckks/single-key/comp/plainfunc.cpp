#include "plainfunc.h"

// void MultipleAdd_SEAL(Evaluator &evaluator, Ciphertext &cipher, Ciphertext& result, long long n) {
// 	long long k, abs_n;
// 	long long *binary;
// //	Ciphertext temp;
// 	if(n>=0) abs_n = n;
// 	else abs_n = -n;

// 	for(k=1; k<100; k++) {
// 		if(abs_n < pow2(k)) break;
// 	}

// 	binary = new long long[k];
// 	for(long i=0; i<k; i++) {
// 		binary[i] = (abs_n / pow2(i)) % 2;
// 	}

// 	evaluator.add(cipher, cipher, result);
// 	if(binary[k-2] == 1) evaluator.add_inplace(result, cipher);

// 	for(long i=k-3; i>=0; i--) {
// 		evaluator.add_inplace(result, result);
// 		if(binary[i] == 1) evaluator.add_inplace(result, cipher);
// 	}

// 	if(n<0) evaluator.negate_inplace(result);
// }


// // In fact, basis_type is meaningless
// void geneT0T1(Encryptor &encryptor, Evaluator &evaluator, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& T0, Ciphertext& T1, Ciphertext& cipher)
// {
// 	double scale = cipher.scale();
// 	long n = cipher.poly_modulus_degree() / 2;
// //	vector<double> m_one(n), m_scaled(n);
// 	vector<double> m_one(n);

// 	// ctxt_1
// 	for(int i=0; i<n; i++) m_one[i] = 1.0; 
// 	Plaintext plain_1;
// 	encoder.encode(m_one, scale, plain_1);
// 	Ciphertext ctxt_1;
// 	encryptor.encrypt(plain_1, ctxt_1);

// 	T0 = ctxt_1;
// 	T1 = cipher;
// }


// void evalT(Evaluator &evaluator, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& Tmplusn, const Ciphertext& Tm, const Ciphertext& Tn, const Ciphertext& Tmminusn)
// {
// 	Ciphertext temp;
// 	evaluator.multiply_reduced_error(Tm, Tn, relin_keys, temp);
// 	evaluator.add_inplace_reduced_error(temp, temp); 
// 	evaluator.rescale_to_next_inplace(temp); 
// 	evaluator.sub_reduced_error(temp, Tmminusn, Tmplusn); 
// }


void evalT_plain(vector<double>& Tmplusn, const vector<double>& Tm, const vector<double>& Tn, const vector<double>& Tmminusn)
{
	vector<double> temp;
	// evaluator.multiply_reduced_error(Tm, Tn, relin_keys, temp);
	temp = mul_v(Tm, Tn);
	// evaluator.add_inplace_reduced_error(temp, temp);
	temp = add_v(temp, temp);
	// evaluator.rescale_to_next_inplace(temp); 
	// evaluator.sub_reduced_error(temp, Tmminusn, Tmplusn); 
	Tmplusn = sub_v(temp, Tmminusn);
}


void eval_polynomial_integrate(vector<double>& res, vector<double>& cipher, long deg, const vector<double> &decomp_coeff, Tree &tree)
{
	// parameter setting and variables

	long n = res.size();
	long total_depth = ceil_to_int(log(static_cast<double>(deg+1))/log(2.0));		// required minimum depth considering both scalar and nonscalar multiplications
	// Ciphertext temp1, temp2, state, ctxt_zero;
	vector<double> temp1, temp2, state, ctxt_zero;
	evaltype eval_type = tree.type;
	vector<long> decomp_deg(pow2(tree.depth+1), -1);
	vector<long> start_index(pow2(tree.depth+1), -1);

	vector<std::unique_ptr<vector<double>>> T(100);
	vector<std::unique_ptr<vector<double>>> pt(100);
	for(size_t i=0; i<100; i++) T[i] = nullptr;
	for(size_t i=0; i<100; i++) pt[i] = nullptr;

	T[0] = std::make_unique<vector<double>>();
	T[1] = std::make_unique<vector<double>>();
	
	// generation of zero ciphertext 
	vector<double> m_coeff(n), m_zero(n, 0.0);


	// set start temp_index
	long num = 0, temp_index;
	if(eval_type == evaltype::oddbaby) temp_index = 1;
	else if(eval_type == evaltype::baby) temp_index = 0;

	// evaluate decompose polynomial degrees
	decomp_deg[1] = deg;
	for(int i=1; i<=tree.depth; i++)
	{
		for(int j=pow2(i); j<pow2(i+1); j++)
		{
			if(j>=static_cast<int>(decomp_deg.size())) throw std::invalid_argument("invalid index");
			if(j%2 == 0) decomp_deg[j] = tree.tree[j/2] - 1;
			else if(j%2 == 1) decomp_deg[j] = decomp_deg[j/2] - tree.tree[j/2];
		}
	}

	// compute start index. 
	for(int i=1; i<pow2(tree.depth+1); i++)
	{
		if(tree.tree[i] == 0)
		{
			start_index[i] = temp_index;
			temp_index += (decomp_deg[i]+1);
		}
	}	

	// generate T0, T1

	// Set T0 to vector of all 1's
	// Set T1 to vector to cipher 
	vector<double> m_one(n);
	for(int i=0; i<n; i++) m_one[i] = 1.0; 
	*T[0] = m_one;
	*T[1] = cipher;

	if(eval_type == evaltype::oddbaby)
	{
		// i: depth stage
		for(int i=1; i<= total_depth+1; i++)
		{
			// cout << "////////////// stage : " << i << endl;

			// depth i computation. all end points. 
			for(int j=1; j<pow2(tree.depth+1); j++)
			{
				if(tree.tree[j] == 0 && total_depth + 1 - num_one(j) == i) 	// depth i stage end points. j: index
				{
					int temp_idx = start_index[j];
					// cout << "pt: " << j << endl;
					pt[j] = std::make_unique<vector<double>>();
					*pt[j] = mult_const(*T[1], decomp_coeff[temp_idx]);

					temp_idx += 2;
					for(int k=3; k<=decomp_deg[j]; k+=2)
					{
						temp1 = mult_const(*T[k], decomp_coeff[temp_idx]);
						*pt[j] = add_v(*pt[j], temp1);
						temp_idx += 2;
					}

				}
			}	

			// depth i computation. all intersection points.
			long inter[40];
			long inter_num = 0;

			for(int j=1; j<pow2(tree.depth+1); j++)
			{
				if(tree.tree[j] > 0 && total_depth + 1 - num_one(j) == i && j%2 == 1) 	// depth i stage intersection points
				{
					long k = j;	
					// cout << "pt: " << j << endl;
					pt[j] = std::make_unique<vector<double>>();
					*pt[j] = mul_v(*T[tree.tree[k]], *pt[2*k+1]);
					k*=2;
					while(1)
					{
						if(tree.tree[k]==0) break;
						temp1 = mul_v(*T[tree.tree[k]], *pt[2*k+1]);
						*pt[j] = add_v(*pt[j], temp1);
						k*=2;
					}

					*pt[j] = add_v(*pt[j], *pt[k]);

					// cout<<"\n\n\npt["<<j;
					// for(auto ptr:*pt[j] ){cout<<ptr<<" ";}
					// cout<<endl;
				}

			}

			// Ti evaluation
			if(i<=tree.m-1) 
			{
				T[pow2(i)] = std::make_unique<vector<double>>();
				evalT_plain(*T[pow2(i)], *T[pow2(i-1)], *T[pow2(i-1)], *T[0]);
			}
			
			if(i<=tree.l)
			{
				for(int j=pow2(i-1)+1; j<=pow2(i)-1; j+=2)		// T1 is not computed. other odd Tis are computed.
				{
					// cout << "T: " << j << endl;
					T[j] = std::make_unique<vector<double>>();
					evalT_plain(*T[j], *T[pow2(i-1)], *T[j-pow2(i-1)], *T[pow2(i)-j]);

				}
			}

		}
		res = *pt[1];
	}
	else if(eval_type == evaltype::baby)
	{
		// i: depth stage
		for(int i=1; i<= total_depth; i++)
		{
			// cout << "////////////// stage : " << i << endl;

			// depth i computation. all end points. 
			for(int j=1; j<pow2(tree.depth+1); j++)
			{
				if(tree.tree[j] == 0 && total_depth + 1 - num_one(j) == i) 	// depth i stage end points. j: index
				{
					int temp_idx = start_index[j];
					// cout << "pt: " << j << endl;
					pt[j] = std::make_unique<vector<double>>();

					*pt[j] = ctxt_zero;
				}
			}

			// depth i computation. all intersection points.
			long inter[40];
			long inter_num = 0;

			for(int j=1; j<pow2(tree.depth+1); j++)
			{
				if(tree.tree[j] > 0 && total_depth + 1 - num_one(j) == i) 	// depth i stage intersection points
				{
					int temp = j;
					bool no_execute = false;
					for(int k=0; k<inter_num; k++)
					{
						while(1)
						{
							if(temp == inter[k]){
								no_execute = true;
								break;
							} 
							if(temp%2 == 0) temp/=2;
							else break;
						}
					}
					
					if(no_execute == false)
					{
						inter[inter_num] = j;
						inter_num += 1;

						long k = j;
						
						// cout << "pt: " << j << endl;
						pt[j] = std::make_unique<vector<double>>();
						if(T[tree.tree[k]] == nullptr) throw std::runtime_error("T[tree.tree[k]] is not set");
						if(pt[2*k+1] == nullptr) throw std::runtime_error("pt[2*k+1] is not set");
						*pt[j] = mul_v(*T[tree.tree[k]], *pt[2*k+1]);
						k*=2;

						while(1)
						{
							if(tree.tree[k]==0) break;
							if(T[tree.tree[k]] == nullptr) throw std::runtime_error("T[tree.tree[k]] is not set");
							if(pt[2*k+1] == nullptr) throw std::runtime_error("pt[2*k+1] is not set");
							temp1 = mul_v(*T[tree.tree[k]], *pt[2*k+1]);
							*pt[j] = add_v(*pt[j], temp1);
							k*=2;
						}

						*pt[j] = add_v(*pt[j], *pt[k]);

					}
				}
			}

			// Ti evaluation
			for(int j=2; j<=tree.b; j++)
			{
				int g = j;
				if(pow2(i-1)<g && g<=pow2(i))
				{
					// cout << "T: " << g << endl;
					T[g] = std::make_unique<vector<double>>();
					if(g%2 == 0) 
					{
						if(T[g/2] == nullptr) throw std::runtime_error("T[g/2] is not set");
						if(T[0] == nullptr) throw std::runtime_error("T[0] is not set");
						evalT_plain(*T[g], *T[g/2], *T[g/2], *T[0]);
					}
					else
					{
						if(T[g/2] == nullptr) throw std::runtime_error("T[g/2] is not set");
						if(T[(g+1)/2] == nullptr) throw std::runtime_error("T[(g+1)/2] is not set");
						if(T[0] == nullptr) throw std::runtime_error("T[0] is not set");
						evalT_plain(*T[g], *T[g/2], *T[(g+1)/2], *T[1]);
					}
				}
			}
			for(int j=1; j<=tree.m-1; j++)
			{
				int g = pow2(j)*tree.b;
				if(pow2(i-1)<g && g<=pow2(i))
				{
					// cout << "T: " << g << endl;
					T[g] = std::make_unique<vector<double>>();
					if(g%2 == 0) 
					{
						if(T[g/2] == nullptr) throw std::runtime_error("T[g/2] is not set");
						if(T[0] == nullptr) throw std::runtime_error("T[0] is not set");
						evalT_plain(*T[g], *T[g/2], *T[g/2], *T[0]);
					}
					else
					{
						if(T[g/2] == nullptr) throw std::runtime_error("T[g/2] is not set");
						if(T[(g+1)/2] == nullptr) throw std::runtime_error("T[(g+1)/2] is not set");
						if(T[0] == nullptr) throw std::runtime_error("T[0] is not set");
						evalT_plain(*T[g], *T[g/2], *T[(g+1)/2], *T[1]);
					}
				}
			}

		}
		res = *pt[1];

	}
}

long coeff_number(long deg, Tree& tree) 
{

	long num = 0;
	long* decomp_deg = new long[pow2(tree.depth+1)];
	decomp_deg[1] = deg;
	for(int i=1; i<=tree.depth; i++)
	{
		for(int j=pow2(i); j<pow2(i+1); j++)
		{
			if(j%2 == 0) decomp_deg[j] = tree.tree[j/2] - 1;
			else if(j%2 == 1) decomp_deg[j] = decomp_deg[j/2] - tree.tree[j/2];
		}
	}

	for(int i=0; i<pow2(tree.depth+1); i++)
	{
		if(tree.tree[i] == 0) 
		{
			num += (decomp_deg[i]+1);
		}
	}
	delete decomp_deg;
	return num;
}

// long ShowFailure_ReLU(Decryptor &decryptor, CKKSEncoder &encoder, Ciphertext& cipher, vector<double>& x, long precision, long n) 
// {
// 	long failure = 0;
// 	double bound = pow(2.0,static_cast<double>(-precision));
// 	Plaintext plain_out;
// 	vector<double> output;
// 	decryptor.decrypt(cipher, plain_out);
// 	encoder.decode(plain_out, output);

// 	for (int i = 0; i < n; ++i) if(abs(ReLU(x[i]) - output[i]) > bound) failure++;

// 	cout << "-------------------------------------------------" << endl;
// 	cout << "failure : " << failure << endl;
// 	cout << "-------------------------------------------------" << endl;
// 	return failure;
// }


vector<double> mult_const(vector<double> a, double c) {
	for (int i = 0; i < a.size(); i++) {
		a[i] *= c;
	}
	return a;
}

vector<double> add_v(vector<double> a, vector<double> b) {
	if(a.size() != b.size()) throw std::invalid_argument("add - len(a) != len(b)");
	vector<double> out(a.size(), 0.0);
	for (int i = 0; i < a.size(); i++) {
		out[i] = a[i] + b[i];
	}
	return out;
}

vector<double> sub_v(vector<double> a, vector<double> b) {
	if(a.size() != b.size()) throw std::invalid_argument("add - len(a) != len(b)");
	vector<double> out(a.size(), 0.0);
	for (int i = 0; i < a.size(); i++) {
		out[i] = a[i] - b[i];
	}
	return out;
}

vector<double> mul_v(vector<double> a, vector<double> b) {
	if(a.size() != b.size()) throw std::invalid_argument("add - len(a) != len(b)");
	vector<double> out(a.size(), 0.0);
	for (int i = 0; i < a.size(); i++) {
		out[i] = a[i] * b[i];
	}
	return out;
}


