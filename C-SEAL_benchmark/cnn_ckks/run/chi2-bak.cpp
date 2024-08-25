/***
 * Copyright (c) 2020-2022 Duality Technologies, Inc.
 * Licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License <https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode>
 * See the LICENSE.md file for the full text of the license.
 * If you share the Licensed Material (including in modified form) you must include the above attribution in the copy you share.
 ***/
/*

Implementation for the Chi-Square GWAS solution described in
"Secure large-scale genome-wide association studies using homomorphic encryption"
by Marcelo Blatt, Alexander Gusev, Yuriy Polyakov, and Shafi Goldwasser

Command to execute the prototype
./demo-chi2 --SNPdir "../data" --SNPfilename "random_sample" --pvalue "pvalue.txt" --runtime "result.txt" --samplesize="200" --snps="16384"

*/

#include <getopt.h>
#include <numeric>
#include <cmath>

// main.cpp
#include "aes_enc.h"
#include <iostream>
#include <omp.h>
#include <algorithm>
#include "Bootstrapper.h"
#include "ModularReducer.h"
using namespace std;
using namespace NTL;
using namespace seal;
using namespace chrono;
using namespace AES;

// 定义一个方便的类型别名以持有时间点
using TimeVar = std::chrono::high_resolution_clock::time_point;

// 定义两个宏，用于开始时间记录 (TIC) 和计算时间消耗 (TOC) 
#define TIC(start) TimeVar start = std::chrono::high_resolution_clock::now()
#define TOC(start) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - (start)).count()


#ifndef GLOBAL_STRUCT_H
#define GLOBAL_STRUCT_H

// 声明一个全局匿名结构体实例
struct {  // SEAL and bootstrapping setting
    long boundary_K = 16; // \approx 1.81 \sqrt(h), e.g. [ 14 \approx 1.81 \sqrt(64) ]
    long boot_deg = 31;
    long scale_factor = 3;
    long inverse_deg = 1; 
    long logN = 13;
    long loge = 10; 
    long logn = 15; 	// full slots
    long logn_1 = 14;	// sparse slots
    long logn_2 = 13;
    long logn_3 = 12;
    long sparse_slots = (1 << logn_1);
    int logp = 42;
    int logq = 58;
    int log_special_prime = 60;
	double scale = pow(2.0, logp);
    int log_integer_part = logq - logp - loge + 5;
    // int log_integer_part = logq - logp;
    int remaining_level = 4; // Calculation required
    int boot_level = 0; // 
    int n_special_prime = 3; //
    int total_level = remaining_level + boot_level;
    size_t secret_key_hamming_weight = 32;
    size_t slot_count;
} SEAL_context_params;

#endif // GLOBAL_STRUCT


SEALContext create_context(int seclevel) {
    if (seclevel != 128) throw std::runtime_error("Security Level not supported");
    // SEAL and bootstrapping setting
    long logN = SEAL_context_params.logN;
    long logn = SEAL_context_params.logn;		// full slots
    long logn_1 = SEAL_context_params.logn_1;	// sparse slots
    long sparse_slots = (1 << logn_1);
    int logp = SEAL_context_params.logp;
    cout << "logp = "<<logp<<endl;
    int logq = SEAL_context_params.logq;
    cout << "logq = "<<logq<<endl;
    int log_special_prime = SEAL_context_params.log_special_prime;
    cout << "log_special_prime = "<<log_special_prime<<endl;
    int n_special_prime = SEAL_context_params.n_special_prime;
    cout << "n_special_prime = "<<n_special_prime<<endl;
    int remaining_level = SEAL_context_params.remaining_level; // Calculation required
    cout << "remaining_level = "<<remaining_level<<endl;
    int boot_level = SEAL_context_params.boot_level; // 
    cout << "boot_level = "<<boot_level<<endl;
    size_t secret_key_hamming_weight = SEAL_context_params.secret_key_hamming_weight;
    cout << "secret_key_hamming_weight = "<<secret_key_hamming_weight<<endl;
    
    vector<int> coeff_bit_vec;
    coeff_bit_vec.push_back(logq);
    for (int i = 0; i < remaining_level; i++) coeff_bit_vec.push_back(logp);
    for (int i = 0; i < boot_level; i++) coeff_bit_vec.push_back(logq);
    for (int i = 0; i < n_special_prime; i++) coeff_bit_vec.push_back(log_special_prime);

    cout << "Setting Parameters" << endl;
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = (size_t)(1 << logN);
    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(
        {0x3ffffffffbe0001, 0x3ffffe80001, 0x3ffffd20001,
         0x3ffffca0001, 0x3ffffbe0001, 0x3ffff4e0001,
         0x3fffefa0001, 0x3fffee60001, 0x3fffe880001,
         0x3fffe820001, 0x3fffe800001, 0x3fffe580001,
         0x3fffe560001, // first modulu and remaining modulus
         0x3fffffffdd80001, 0x3ffffffff3a0001, 0x3ffffffff040001,
         0x3fffffffed60001, 0x3fffffffed00001, 0x3fffffffeb00001,
         0x3fffffffea00001, 0x3fffffffe800001, 0x3fffffffe440001,
         0x3fffffffe320001, 0x3fffffffe2c0001, 0x3fffffffdfe0001,
         0x7ffffffffcc0001, 0x7ffffffffba0001, 0x7ffffffffb00001,
         0xffffffffffc0001, 0xfffffffff840001});  //

    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec)); 
    
    // modified SEAL
    parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
    parms.set_n_special_primes(n_special_prime);
    // parms.set_sparse_slots(sparse_slots);
    return SEALContext(parms, true, sec_level_type::none);
}

const double EPSILON = 1.0E-08;

void RunChi2(const string &SNPDir, const string &SNPFileName,
		const string &pValue, const string &Runtime, const string &SampleSize, const string &SNPs);

Ciphertext BinaryTreeAdd(std::vector<Ciphertext> &vector, Evaluator &evaluator);

void ReadSNPFile(vector<string>& headers, std::vector<std::vector<double>> & dataColumns, std::vector<double> &y, string dataFileName, size_t N, size_t M);

double BS(double z) {
	double y = exp(-z*z/2);
	return sqrt(1-y) * (31*y/200 - 341*y*y/8000) / sqrt(M_PI);
}

double normalCFD(double value) { return 0.5 * erfc(-value * M_SQRT1_2); }

double sf(double value) { return 1 - normalCFD(value); }

static bool Equal(double a, double b) {   return (EPSILON > fabs(a-b)); }

static bool Less(double a, double b) {  return ((a-b) < (-EPSILON)); }

static bool Greater(double a, double b) {   return ((a-b) > EPSILON); }

double IncompleteGamma(double val, double p);

int main(int argc, char **argv) {

	int opt;

	static struct option long_options[] =
	  {
		/* These options dont set a flag.
		   We distinguish them by their indices. */
		{"SNPdir",  	required_argument, 			0, 'S'},
		{"SNPfilename",  	required_argument, 			0, 's'},
		{"pvalue",  	required_argument, 			0, 'p'},
		{"runtime",  	required_argument, 			0, 'r'},
		{"samplesize",  	required_argument, 			0, 'N'},
		{"snps",  	required_argument, 			0, 'M'},
		{0, 0, 0, 0}
	  };

	/* getopt_long stores the option index here. */
	int option_index = 0;

	string SNPDir;
	string SNPFileName;
	string pValue;
	string Runtime;
	string SampleSize;
	string SNPs;

	while ((opt = getopt_long(argc, argv, "S:s:p:r:N:M", long_options, &option_index)) != -1) {
		switch (opt)
		{
			case 'S':
				SNPDir = optarg;
				break;
			case 's':
				SNPFileName = optarg;
				break;
			case 'p':
				pValue = optarg;
				break;
			case 'r':
				Runtime = optarg;
				break;
			case 'N':
				SampleSize = optarg;
				break;
			case 'M':
				SNPs = optarg;
				break;
			default: /* '?' */
			  std::cerr<< "Usage: "<<argv[0]<<" <arguments> " <<std::endl
				   << "arguments:" <<std::endl
				   << "  -S --SNPdir SNP file directory"  <<std::endl
				   << "  -s --SNPfilename SNP file name"  <<std::endl
				   << "  -o --pvalue p-values file"  <<std::endl
				   << "  -r --runtime runtime output file name"  <<std::endl
				   << "  -N --samplesize number of individuals"  <<std::endl
			   	   << "  -M --snps number of SNPs"  <<std::endl;
			  exit(EXIT_FAILURE);
		}
	}

	RunChi2(SNPDir, SNPFileName, pValue, Runtime, SampleSize, SNPs);

	return 0;

}

void RunChi2(const string &SNPDir,
		const string &SNPFileName, const string &pValue, const string &Runtime, const string &SampleSize, const string &SNPs) {

	

	TIC(tAll);

	double keyGenTime(0.0);
	double encryptionTime(0.0);
	double computationTime(0.0);
	double decryptionTime(0.0);
	double endToEndTime(0.0);

	std::cout << "\n======CHI-SQUARE SOLUTION========\n" << std::endl;

	vector<string> headers1;
	vector<string> headersS;

	std::vector<double> yData;
	std::vector<std::vector<double>> sData;

	size_t N = std::stoi(SampleSize);
	size_t M = std::stoi(SNPs);

	ReadSNPFile(headersS,sData,yData,SNPDir + "/" + SNPFileName,N,M);

	size_t m = 16384;

	size_t init_size = 4;
	double scalingFactor = 2.5e-6;

	//Gen SEAL Context
    omp_set_num_threads(64);
    SEALContext context = create_context(128);//create a context
    long loge = SEAL_context_params.loge;
    long logn = SEAL_context_params.logn; 
    long logn_1 = SEAL_context_params.logn_1;
    long logN = SEAL_context_params.logN;
    int total_level = SEAL_context_params.total_level;
    int remaining_level = SEAL_context_params.remaining_level;
    double scale = SEAL_context_params.scale;
    long boundary_K = SEAL_context_params.boundary_K;
    long boot_deg = SEAL_context_params.boot_deg;
    long scale_factor = SEAL_context_params.scale_factor;
    long inverse_deg = SEAL_context_params.inverse_deg;

	TIC(kg_time);

	// key generate
	KeyGenerator keygen(context);
    PublicKey he_pk;
	keygen.create_public_key(he_pk);
	auto he_sk = keygen.secret_key();
    RelinKeys he_rk;
	keygen.create_relin_keys(he_rk);
	GaloisKeys he_gk;

	keyGenTime = TOC(kg_time);
	

	CKKSEncoder encoder(context);
	Encryptor encryptor(context, he_pk);
	// Evaluator evaluator(context);
	Evaluator evaluator(context, encoder);
	Decryptor decryptor(context, he_sk);
	size_t slot_count = encoder.slot_count();

    vector<uint8_t> sk = {
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00
    };

	// AES_ENC tester(sk, remaining_level, loge, logn, logN - 1, total_level, scale, boundary_K, boot_deg, scale_factor, inverse_deg,
    //  context, keygen, encoder, encryptor, decryptor, evaluator, he_rk, he_gk);
	// Gen SEAL Context Finished
	// vector<Ciphertext> ct = tester.debug_test(sk, 128);
	// decoder of the HE Ciphertext
	// DECODER decoder( data_len, context, keygen, encoder, encryptor, decryptor, evaluator, he_rk, he_gk  );
	// vector<Ciphertext> ct = decoder.decode(ct, 128);
	// decode finished


	TIC(enc_time);

	std::cout << "\nNumber of Individuals = " << sData.size() << std::endl;
	std::cout << "Number of SNPs = " << sData[0].size() << std::endl;


	size_t sizeS = (size_t)std::ceil((double)sData[0].size()/(m/4));

	std::vector<std::vector<std::vector<std::complex<double>>>> sDataArray(sizeS);

	for(size_t s = 0; s < sizeS; s++)
		sDataArray[s] = std::vector<std::vector<std::complex<double>>>(sData.size());

	for (size_t i=0; i < sData.size(); i++){

		for(size_t s = 0; s < sizeS; s++)
			sDataArray[s][i] = std::vector<std::complex<double>>(m/4);

		size_t counter = 0;

		for (size_t j=0; j<sData[i].size(); j++) {
			if ((j>0) && (j%(m/4)==0))
				counter++;
			sDataArray[counter][i][j%(m/4)] = sData[i][j];
		}
	}

	//std::cerr << " sData = " << sDataArray << std::endl;

	std::vector<std::vector<Ciphertext>> S(sizeS);
	std::vector<Ciphertext> Y(N);

	for (size_t i = 0; i < sizeS; i++)
		S[i] = std::vector<Ciphertext>(N);

	//Encryption of single-integer ciphertexts
#pragma omp parallel for
	for (size_t i=0; i<N; i++){
		for (size_t s=0; s < sizeS; s++){
			Plaintext sTemp;
			encoder.encode(sDataArray[s][i], scale, sTemp);
			encryptor.encrypt(sTemp, S[s][i]);			
		}
		Plaintext sTemp2;
		encoder.encode(std::vector<std::complex<double>>(m/4,yData[i]), scale, sTemp2);
		encryptor.encrypt(sTemp2, Y[i]);	
	}
	encryptionTime = TOC(enc_time);

	TIC(cpt_time);

	Plaintext d;
	encoder.encode( std::vector<std::complex<double>>(m/4,2*N), scale, d );

	Plaintext dScaled;
	encoder.encode( std::vector<std::complex<double>>(m/4,2*N*scalingFactor), scale, dScaled );

	std::vector<Ciphertext> ySum = Y;

	auto yU = BinaryTreeAdd(ySum, evaluator);

	std::vector<Ciphertext> chiD(sizeS);
	std::vector<Ciphertext> chiN(sizeS);

	std::vector<Ciphertext> orD(sizeS);
	std::vector<Ciphertext> orN(sizeS);

	for (size_t s = 0; s < sizeS; s++) {

		std::vector<Ciphertext> ySMult(N);

#pragma omp parallel for
		for(size_t i = 0; i < N; i++) {
			evaluator.multiply(S[s][i], Y[i], ySMult[i]);
			evaluator.relinearize_inplace(ySMult[i], he_rk);
			evaluator.rescale_to_next_inplace(ySMult[i]);
		}
		auto n11 = BinaryTreeAdd(ySMult, evaluator);

		auto c1 = BinaryTreeAdd(S[s], evaluator);
		Ciphertext r1 = yU;	
		evaluator.add_inplace_reduced_error(r1, yU);

		Ciphertext r1Scaled;
		evaluator.multiply_const(r1, scalingFactor, r1Scaled);
		evaluator.rescale_to_next_inplace(r1Scaled);

		Ciphertext c1Scaled;
		evaluator.multiply_const(c1, scalingFactor, c1Scaled);
		evaluator.rescale_to_next_inplace(c1Scaled);

		// Chi2 computation
		// numerator
		Ciphertext mult1;
		Plaintext tmp;
		evaluator.mod_switch_to(dScaled, n11.parms_id(), tmp);
		evaluator.multiply_plain(n11, tmp, mult1);
		evaluator.rescale_to_next_inplace(mult1);
		Ciphertext mult2;
		evaluator.multiply_reduced_error(c1, r1Scaled, he_rk, mult2);
		evaluator.rescale_to_next_inplace(mult2);
		Ciphertext chiN1;
		evaluator.sub_reduced_error(mult1, mult2, chiN1);
		chiN[s] = chiN1;
		evaluator.multiply_inplace_reduced_error(chiN[s], chiN1, he_rk);
		evaluator.rescale_to_next_inplace(chiN[s]);
		// denominator
		Ciphertext chiD1;
		evaluator.mod_switch_to(dScaled, c1Scaled.parms_id(), tmp);
		c1Scaled.scale() = tmp.scale();
		evaluator.sub_plain(c1Scaled, tmp, chiD1);
		evaluator.multiply_inplace_reduced_error(chiD1, c1, he_rk);
		evaluator.rescale_to_next_inplace(chiD1);

		Ciphertext chiD2;
		evaluator.mod_switch_to(dScaled, r1Scaled.parms_id(), tmp);
		r1Scaled.scale() = dScaled.scale();
		evaluator.sub_plain(r1Scaled, tmp, chiD2);
		evaluator.multiply_inplace_reduced_error(chiD2, r1, he_rk);
		evaluator.rescale_to_next_inplace(chiD2);
		evaluator.multiply_reduced_error(chiD1, chiD2, he_rk, chiD[s]);
		evaluator.rescale_to_next_inplace(chiD[s]);
		// Odds Ratio Computation
		Ciphertext n11Scaled;
		evaluator.multiply_const(n11, scalingFactor, n11Scaled);
		evaluator.rescale_to_next_inplace(n11Scaled);

		// denominator

		Ciphertext or2;
		evaluator.sub_reduced_error(c1, n11, or2);

		Ciphertext or3;
		evaluator.sub_reduced_error(r1Scaled, n11Scaled, or3);
		evaluator.multiply_reduced_error(or2, or3, he_rk, orD[s]);
		evaluator.rescale_to_next_inplace(orD[s]);
		// numerator
		Ciphertext or1;
		evaluator.sub_reduced_error(n11Scaled, r1Scaled, or1);
		evaluator.sub_inplace_reduced_error(or1, c1Scaled);

		evaluator.mod_switch_to(dScaled, or1.parms_id(), tmp);
		or1.scale() = tmp.scale();
		evaluator.add_plain_inplace(or1, tmp);
		evaluator.multiply_reduced_error(n11, or1, he_rk, orN[s]);
		evaluator.rescale_to_next_inplace(orN[s]);
	}

	computationTime = TOC(cpt_time);

	TIC(dec_time);

	std::vector<Plaintext> pN(sizeS);
	std::vector<Plaintext> pD(sizeS);
	std::vector<Plaintext> oddN(sizeS);
	std::vector<Plaintext> oddD(sizeS);

	for (size_t s = 0; s < sizeS; s++) {
		
		decryptor.decrypt( chiN[s], pN[s] );
		decryptor.decrypt( chiD[s], pD[s] );
		decryptor.decrypt( orN[s], oddN[s] );
		decryptor.decrypt( orD[s], oddD[s] );
	}

	decryptionTime = TOC(dec_time);

	std::vector<double> chival(headersS.size());
	std::vector<double> pval(headersS.size());
	std::vector<double> odds(headersS.size());

	for (size_t s = 0; s < sizeS; s++) {
		vector<complex<double> > msg_pN;
		vector<complex<double> > msg_pD;
		vector<complex<double> > msg_oddN;
		vector<complex<double> > msg_oddD;
		encoder.decode(pN[s], msg_pN);
		encoder.decode(pD[s], msg_pD);
		encoder.decode(oddN[s], msg_oddN);
		encoder.decode(oddD[s], msg_oddD);
		
		for (size_t i = 0; i < m/4; i++) {
			if (s*m/4 + i < headersS.size()) {
				
				chival[s*m/4 + i] = msg_pN[i].real()*2*N / msg_pD[i].real();
				if (chival[s*m/4 + i] < 0)
					chival[s*m/4 + i] = 0;
				pval[s*m/4 + i] = (double)1-IncompleteGamma(chival[s*m/4 + i]/2,0.5);
				if (pval[s*m/4 + i] < 0)
					pval[s*m/4 + i] = 1e-15;
				else
					if (pval[s*m/4 + i]==0)
						pval[s*m/4 + i] = BS(sqrt(chival[s*m/4 + i]));

				odds[s*m/4 + i] = msg_oddN[i].real() / msg_oddD[i].real();

			}
		}
	}

    ofstream myfile;
    myfile.open(SNPDir + "/" + pValue);
    myfile.precision(10);
    for(uint32_t i = 0; i < headersS.size(); i++) {
    	myfile << headersS[i] << "\t" << pval[i] << std::endl;
    }
    myfile.close();

    ofstream myfile2;
    myfile2.open(SNPDir + "/" + "odds.txt");
    myfile2.precision(10);
    for(uint32_t i = 0; i < headersS.size(); i++) {
    	myfile2 << headersS[i] << "\t" << odds[i] << std::endl;
    }
    myfile2.close();

    ofstream myfile3;
    myfile3.open(SNPDir + "/" + "chi2.txt");
    myfile3.precision(10);
    for(uint32_t i = 0; i < headersS.size(); i++) {
    	myfile3 << headersS[i] << "\t" << chival[i] << std::endl;
    }
    myfile3.close();

	std::cout << "\nKey Generation Time: \t\t" << keyGenTime/1000 << " s" << std::endl;
	std::cout << "Encoding and Encryption Time: \t" << encryptionTime/1000 << " s" << std::endl;
	std::cout << "Computation Time: \t\t" << computationTime/1000 << " s" << std::endl;
	std::cout << "Decryption & Decoding Time: \t" << decryptionTime/1000 << " s" << std::endl;

	endToEndTime = TOC(tAll);

    std::cout << "\nEnd-to-end Runtime: \t\t" << endToEndTime/1000 << " s" << std::endl;

    ofstream myfileRuntime;
    myfileRuntime.open(SNPDir + "/" + Runtime);
    myfileRuntime << "Key Generation Time: \t\t" << keyGenTime/1000 << " s" << std::endl;
    myfileRuntime << "Encoding and Encryption Time: \t" << encryptionTime/1000 << " s" << std::endl;
    myfileRuntime << "Computation Time: \t\t" << computationTime/1000 << " s" << std::endl;
    myfileRuntime << "Decryption & Decoding Time: \t" << decryptionTime/1000 << " s" << std::endl;
    myfileRuntime << "End-to-end Runtime: \t\t" << endToEndTime/1000 << " s" << std::endl;
    myfileRuntime.close();

}

Ciphertext BinaryTreeAdd(std::vector<Ciphertext> &vector, Evaluator &evaluator) {

	for(size_t j = 1; j < vector.size(); j=j*2) {
		for(size_t i = 0; i<vector.size(); i = i + 2*j) {
			if ((i+j)<vector.size())
				evaluator.add_inplace_reduced_error(vector[i],vector[i+j]);
		}
	}
	return vector[0];
}

double IncompleteGamma(double val, double p)
{
    if( !Greater(val, 0) || !Greater(p, 0) )
        return 0;

    double expValue = p*log(val) - val - lgamma(p);
    if( Less(expValue, log(1.0E-37)) ) // undeflow
        return 0;

    double factor = exp(expValue);
    if( !Greater(val, 1) || Less(val, p) )
    {
        double igamma = 1;
        double term = 1;

        for( int i = 1; Greater(term, EPSILON); ++i )
        {
            term *= (val/(p+i));
            igamma += term;
        }

        return (igamma*factor/p);
    }

    double pn[6] = { 1, val, val+1, val*(2+val-p) };
    double upperIncGamma = pn[2]/pn[3];

    for( int j = 1; ; ++j )
    {
        double a = (j+1)*2 + val- p;
        double b = (1 + j - p)*j;
        pn[4] = a*pn[2] - b*pn[0];
        pn[5] = a*pn[3] - b*pn[1];

        if( !Equal(pn[5], 0) )
        {
            double rn = pn[4]/pn[5];
            double diff = fabs(upperIncGamma - rn);
            if( !Greater(diff, EPSILON) && !Greater(diff, (EPSILON*rn)) )
                return (1 - factor*upperIncGamma);

            upperIncGamma = rn;
        }

        for( int i = 0; i < 4; i++ )
            pn[i] = pn[i+2];

        if( !Greater(1.0E+37, fabs(pn[4])) ) // overflow
        {
            for( int i = 0; i < 4; i++ )
                pn[i] = pn[i] / 1.0E+37;
        }
    }

    return 0;
}

void ReadSNPFile(vector<string>& headers, std::vector<std::vector<double>> & dataColumns, std::vector<double> &y,
                 string dataFileName, size_t N, size_t M)
{

	uint32_t cols = 0;

	string fileName = dataFileName + ".csv";

	std::cerr << "file name = " << fileName << std::endl;

	ifstream file(fileName);
	string line, value;

	if(file.good()) {

		getline(file, line);
		cols = std::count(line.begin(), line.end(), ',') + 1;
		stringstream ss(line);
		vector<string> result;

		size_t tempcounter = 0;

		for(uint32_t i = 0; i < cols; i++) {
			string substr;
			getline(ss, substr, ',');
			if ((substr != "") && (i>4) && (i<M+5)) {
				headers.push_back(substr);
				tempcounter++;
			}
		}

		cols = tempcounter;

	}

	size_t counter = 0;
	while((file.good()) && (counter < N)) {
		getline(file, line);
		uint32_t curCols = std::count(line.begin(), line.end(), ',') + 1;
		if (curCols > 2) {
			stringstream ss(line);
			for(uint32_t i = 0; i < 5; i++) {
				string substr;
				getline(ss, substr, ',');
				if ((i==1))
					y.push_back(std::stod(substr));
			}
			std::vector<double> row(cols);
			for(uint32_t i = 5; i < cols + 5; i++) {
				string substr;
				getline(ss, substr, ',');
				if (i < M+5)
				{
					double val;
					val = std::stod(substr);
					row[i-5] = val;
				}
			}
			dataColumns.push_back(row);
		}
		counter++;
	}

	file.close();

    std::cout << "Read in data: ";
	std::cout << dataFileName << std::endl;
}

