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
./chi2 --SNPdir "../data" --SNPfilename "random_sample" --pvalue "pvalue.txt" --runtime "result.txt" --samplesize="200" --snps="16384"

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

using TimeVar = std::chrono::high_resolution_clock::time_point;

#define TIC(start) TimeVar start = std::chrono::high_resolution_clock::now()
#define TOC(start) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - (start)).count()


#ifndef GLOBAL_STRUCT_H
#define GLOBAL_STRUCT_H

struct {  // SEAL and bootstrapping setting
    long boundary_K = 16; // \approx 1.81 \sqrt(h), e.g. [ 14 \approx 1.81 \sqrt(64) ]
    long boot_deg = 31;
    long scale_factor = 3;
    long inverse_deg = 1; 
    long logN = 16;
    long loge = 10; 
    long logn = 15; 	// full slots
    long logn_1 = 14;	// sparse slots
    long logn_2 = 13;
    long logn_3 = 12;
    long sparse_slots = (1 << logn_1);
    int logp = 42;
    int logq = 50;
    int log_special_prime = 50;
	double scale = pow(2.0, logp);
    int log_integer_part = logq - logp - loge + 5;
    // int log_integer_part = logq - logp;
    int remaining_level = 12; // Calculation required
    int boot_level = 0; // 
    int n_special_prime = 5; //
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

    // parms.set_coeff_modulus(
    //     {0x3ffffffffbe0001, 0x3ffffe80001, 0x3ffffd20001,
    //      0x3ffffca0001, 0x3ffffbe0001, 0x3ffff4e0001,
    //      0x3fffefa0001, 0x3fffee60001, 0x3fffe880001,
    //      0x3fffe820001, 0x3fffe800001, 0x3fffe580001,
    //      0x3fffe560001, // first modulu and remaining modulus
    //      0x3fffffffdd80001, 0x3ffffffff3a0001, 0x3ffffffff040001,
    //      0x3fffffffed60001, 0x3fffffffed00001, 0x3fffffffeb00001,
    //      0x3fffffffea00001, 0x3fffffffe800001, 0x3fffffffe440001,
    //      0x3fffffffe320001, 0x3fffffffe2c0001, 0x3fffffffdfe0001,
    //      0x7ffffffffcc0001, 0x7ffffffffba0001, 0x7ffffffffb00001,
    //      0xffffffffffc0001, 0xfffffffff840001
	// 	 });  //

    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec)); 
    
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

	for(int i=0;i<10;i++){
		cout <<headersS[i]<<endl;
	}

	size_t m = 16384;

	// size_t init_size = 4;
	double appdata_encode_scale = pow(2, 30);
	double scalingFactor = 0.1*pow(N, -2);

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
	cout << "slot count: "<<slot_count<<endl;

	vector<int> gal_steps_vector = { (int)slot_count/2 };
	keygen.create_galois_keys(gal_steps_vector, he_gk);

    vector<uint8_t> sk = {
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00
    };

	TIC(enc_time);

	std::cout << "\nNumber of Individuals = " << sData.size() << std::endl;
	std::cout << "Number of SNPs = " << sData[0].size() << std::endl;

	// //print the input data for debug
	// for(int i=0;i<N;i++){
	// 	cout << "Individual[" <<i<<"]: SNP: "<<endl;
	// 	for(int j=0;j<3;j++){
	// 		cout << "SNP["<<j<<"]: "<<sData[i][j]<<endl;
	// 	}
	// }
	// //print the input data for debug
	// for(int i=0;i<N;i++){
	// 	cout << "yData["<<i<<"]: "<<yData[i]<<endl;
		
	// }

	// block number = N*M/4, due to one 128bit block can store four 32-bit float value
	size_t block_num = (N/2) * M / 4;
	// formated_input_data[SNP][N]
	std::vector< std::vector<uint8_t> > formated_input_data(block_num * 2);
	// half of the data block, first slots / 2
	for(size_t i=0;i<block_num;i++){
		vector<uint8_t> data_block;
		for(size_t j=0;j<4;j++){//data is uint_8 form 128/8=16
			// 32bit scaled float data
			uint32_t SNP_value = (uint32_t) appdata_encode_scale * sData[j][i];
			//transform the 32bit data into 8bit
			for(size_t k=0;k<4;k++){
				uint8_t tmp = SNP_value & 0xff;
				SNP_value = SNP_value >> 8;
				data_block.push_back(tmp);
			}
		}
		formated_input_data[i] = data_block;
	}

	// half of the data block, last slots / 2
	for(size_t i=0;i<block_num;i++){
		vector<uint8_t> data_block;
		for(size_t j=0;j<4;j++){//data is uint_8 form 128/8=16
			// 32bit scaled float data
			uint32_t SNP_value = (uint32_t) appdata_encode_scale * sData[j+4][i];
			//transform the 32bit data into 8bit
			for(size_t k=0;k<4;k++){
				uint8_t tmp = SNP_value & 0xff;
				SNP_value = SNP_value >> 8;
				data_block.push_back(tmp);
			}
		}
		formated_input_data[block_num+i] = data_block;
	}
	AES_ENC transcipher(sk, remaining_level, loge, logn, logN - 1, total_level, scale, boundary_K, boot_deg, scale_factor, inverse_deg,
     context, keygen, encoder, encryptor, decryptor, evaluator, he_rk, he_gk);
#if 0

	// Gen SEAL Context Finished

	// decoder of the HE Ciphertext
	//reformat the input data
	vector<uint8_t> pt_encode;
	//formated_input_data[SNP+SNP][N+2N]
	for(size_t i=0;i<formated_input_data.size();i++){
		for(size_t j=0;j<formated_input_data[0].size();j++){
			// if (i<3) {
			// 	uint8_t pt = formated_input_data[i][j];
			// 	// cout << "formated_: SNP "<<i<<"N: "<<j<<"value: "<<(int) pt<<endl;
			// }
			pt_encode.push_back( formated_input_data[i][j] );
		}
	}
	transcipher.encode_ciphertext(pt_encode, block_num*2);
	// vector<Ciphertext> ct_txt = transcipher.HE_decrypt(sk, 128*block_num*2);
	vector<Ciphertext> ct_txt = transcipher.get_encoded_ct();

	std::vector<Ciphertext> S;
	for(size_t i=0;i<4;i++){
		vector<Ciphertext> bits( ct_txt.begin()+32*i, ct_txt.begin()+32*(i+1) );
		Ciphertext rtncipher = transcipher.construct_number_from_bits(bits, 10);
		evaluator.multiply_const_inplace(rtncipher, 1/appdata_encode_scale);//scale back to input
		evaluator.rescale_to_next_inplace(rtncipher);
		S.push_back(rtncipher);
	}

	for(int i=0;i<S.size();i++){
		transcipher.debugPrint(S[i], "S: ");
	}
	
	// vector<Ciphertext> ct = decoder.decode(ct, 128);
	// decode finished
#else
	encryptionTime = TOC(enc_time);
	std::vector<Ciphertext> S(N/2);
	
	for(size_t i=0;i<N/2;i++){
		Plaintext sTemp2;
		vector<complex<double> > individual;
		for(size_t j=0;j<M;j++){
			individual.push_back(sData[i][j]);
		}
		for(size_t j=0;j<M;j++){
			individual.push_back(sData[N/2+i][j]);
		}
		encoder.encode(individual, scale, sTemp2);
		encryptor.encrypt(sTemp2, S[i]);
	}
#endif
	std::vector<Ciphertext> Y(N/2);
	for(size_t i=0;i<N/2;i++){
		Plaintext sTemp2;
		std::vector<std::complex<double>> tmp(m, yData[i]);
		std::vector<std::complex<double>> tmp2(m, yData[N/2+i]);
		tmp.insert(tmp.end(), tmp2.begin(), tmp2.end());
		encoder.encode(tmp, scale, sTemp2);
		encryptor.encrypt(sTemp2, Y[i]);
	}

	TIC(cpt_time);

	Plaintext d;
	encoder.encode( std::vector<std::complex<double>>(m,2*N), scale, d );

	Plaintext dScaled;
	encoder.encode( std::vector<std::complex<double>>(m,2*N*scalingFactor), scale, dScaled );

	std::vector<Ciphertext> ySum = Y;

	auto yU = BinaryTreeAdd(ySum, evaluator);


	Ciphertext chiD;
	Ciphertext chiN;

	Ciphertext orD;
	Ciphertext orN;

	std::vector<Ciphertext> ySMult(N/2);

#pragma omp parallel for
	for(size_t i = 0; i < N/2; i++) {//two individual are packed in one ct
		evaluator.multiply_reduced_error(S[i], Y[i], he_rk, ySMult[i]);
		evaluator.rescale_to_next_inplace(ySMult[i]);
	}
	auto n11 = BinaryTreeAdd(ySMult, evaluator);
	auto c1 = BinaryTreeAdd(S, evaluator);
	//rotate and add
	int step = slot_count/2;
	Ciphertext tmp_n11 = n11;
	evaluator.rotate_vector(n11, step, he_gk, tmp_n11);
	evaluator.add_inplace_reduced_error(n11, tmp_n11);
	
	Ciphertext tmp_c1 = c1;
	evaluator.rotate_vector(c1, step, he_gk, tmp_c1);
	evaluator.add_inplace_reduced_error(c1, tmp_c1);

	Ciphertext tmp_yU = yU;
	evaluator.rotate_vector(yU, step, he_gk, tmp_yU);
	evaluator.add_inplace_reduced_error(yU, tmp_yU);

	transcipher.debugPrint(n11, "n11: ");
	transcipher.debugPrint(c1, "c1: ");
	transcipher.debugPrint(yU, "yU: ");

	//r1 = 2 yU
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
	tmp.scale() = n11.scale();
	evaluator.multiply_plain(n11, tmp, mult1);
	evaluator.rescale_to_next_inplace(mult1);

	Ciphertext mult2;
	evaluator.multiply_reduced_error(c1, r1Scaled, he_rk, mult2);
	evaluator.rescale_to_next_inplace(mult2);

	Ciphertext chiN1;
	evaluator.sub_reduced_error(mult1, mult2, chiN1);

	chiN = chiN1;
	evaluator.multiply_inplace_reduced_error(chiN, chiN1, he_rk);
	evaluator.rescale_to_next_inplace(chiN);


	// denominator
	Ciphertext chiD1, negete_c1Scaled;
	evaluator.negate(c1Scaled, negete_c1Scaled);
	evaluator.mod_switch_to(dScaled, negete_c1Scaled.parms_id(), tmp);
	tmp.scale() = negete_c1Scaled.scale();
	evaluator.add_plain(negete_c1Scaled, tmp, chiD1);
	evaluator.multiply_inplace_reduced_error(chiD1, c1, he_rk);
	evaluator.rescale_to_next_inplace(chiD1);

	Ciphertext chiD2, negete_r1Scaled;
	evaluator.negate(r1Scaled, negete_r1Scaled);
	evaluator.mod_switch_to(dScaled, negete_r1Scaled.parms_id(), tmp);
	tmp.scale() = negete_r1Scaled.scale();
	evaluator.add_plain(negete_r1Scaled, tmp, chiD2);
	evaluator.multiply_inplace_reduced_error(chiD2, r1, he_rk);

	evaluator.rescale_to_next_inplace(chiD2);
	evaluator.multiply_reduced_error(chiD1, chiD2, he_rk, chiD);
	evaluator.rescale_to_next_inplace(chiD);
	// Odds Ratio Computation
	Ciphertext n11Scaled;
	evaluator.multiply_const(n11, scalingFactor, n11Scaled);
	evaluator.rescale_to_next_inplace(n11Scaled);
	// denominator
	Ciphertext or2;
	evaluator.sub_reduced_error(c1, n11, or2);

	Ciphertext or3;
	evaluator.sub_reduced_error(r1Scaled, n11Scaled, or3);
	evaluator.multiply_reduced_error(or2, or3, he_rk, orD);
	evaluator.rescale_to_next_inplace(orD);
	// numerator
	Ciphertext or1;
	evaluator.sub_reduced_error(n11Scaled, r1Scaled, or1);
	evaluator.sub_inplace_reduced_error(or1, c1Scaled);

	evaluator.mod_switch_to(dScaled, or1.parms_id(), tmp);
	tmp.scale() = or1.scale();
	evaluator.add_plain_inplace(or1, tmp);
	evaluator.multiply_reduced_error(n11, or1, he_rk, orN);
	evaluator.rescale_to_next_inplace(orN);

	transcipher.debugPrint(chiN, "chiN: ");
	transcipher.debugPrint(chiD, "chiD: ");
	transcipher.debugPrint(orN, "orN: ");
	transcipher.debugPrint(orD, "orD: ");
	
	computationTime = TOC(cpt_time);

	TIC(dec_time);

	Plaintext pN;
	Plaintext pD;
	Plaintext oddN;
	Plaintext oddD;

	decryptor.decrypt( chiN, pN );
	decryptor.decrypt( chiD, pD );
	decryptor.decrypt( orN, oddN );
	decryptor.decrypt( orD, oddD );

	decryptionTime = TOC(dec_time);

	std::vector<double> chival(headersS.size());
	std::vector<double> pval(headersS.size());
	std::vector<double> odds(headersS.size());

	vector<complex<double> > msg_pN;
	vector<complex<double> > msg_pD;
	vector<complex<double> > msg_oddN;
	vector<complex<double> > msg_oddD;
	encoder.decode(pN, msg_pN);
	encoder.decode(pD, msg_pD);
	encoder.decode(oddN, msg_oddN);
	encoder.decode(oddD, msg_oddD);
		
	for (size_t i = 0; i < m; i++) {
		if (i < headersS.size()) {
			
			chival[i] = msg_pN[i].real()*2*N / msg_pD[i].real();
			if (chival[i] < 0)
				chival[i] = 0;
			pval[i] = (double)1-IncompleteGamma(chival[i]/2,0.5);
			if (pval[i] < 0)
				pval[i] = 1e-15;
			else
				if (pval[i]==0)
					pval[i] = BS(sqrt(chival[i]));
			odds[i] = msg_oddN[i].real() / msg_oddD[i].real();
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

