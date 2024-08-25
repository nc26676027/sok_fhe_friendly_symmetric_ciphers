// main.cpp
#include "gift64_enc.h"
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <thread>
#include <chrono>
#include "Bootstrapper.h"
#include "ModularReducer.h"
using namespace std;
using namespace NTL;
using namespace seal;
using namespace chrono;
using namespace GIFT64;

#ifndef GLOBAL_STRUCT_H
#define GLOBAL_STRUCT_H

// 声明一个全局匿名结构体实例
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
    int logq = 58;
    int log_special_prime = 60;
	double scale = pow(2.0, logp);
    int log_integer_part = logq - logp - loge + 5;
    // int log_integer_part = logq - logp;
    int remaining_level = 12; // Calculation required
    int boot_level = 12; // 
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
    coeff_bit_vec.push_back(log_special_prime);


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
    //      0xffffffffffc0001, 0xfffffffff840001});  //
    // parms.set_coeff_modulus(
    //     {0x3ffffffffbe0001, 0x3ffffe80001, 0x3ffffd20001,
    //      0x3ffffca0001, 0x3ffffbe0001, 0x3ffff4e0001,
    //      0x3fffefa0001, 0x3fffee60001, 0x3fffe880001,
    //      0x3fffe820001, 0x3fffe800001, 
    //     //  0x3fffe580001,
    //     //  0x3fffe560001, // first modulu and remaining modulus
    //      0x3fffffffdd80001, 0x3ffffffff3a0001, 0x3ffffffff040001,
    //      0x3fffffffed60001, 0x3fffffffed00001, 0x3fffffffeb00001,
    //      0x3fffffffea00001, 0x3fffffffe800001, 0x3fffffffe440001,
    //      0x3fffffffe320001, 0x3fffffffe2c0001, 0x3fffffffdfe0001,
    //      0x7ffffffffcc0001, 0x7ffffffffba0001,
    //      0xffffffffffc0001});  //
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec)); 
    
    // modified SEAL
    parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
    // parms.set_sparse_slots(sparse_slots);
    auto sec = seal::sec_level_type::none;
    return SEALContext(parms, true, sec);
}


int main() {
    omp_set_num_threads(1);
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

	// key generate
	KeyGenerator keygen(context);
    PublicKey he_pk;
	keygen.create_public_key(he_pk);
	auto he_sk = keygen.secret_key();
    RelinKeys he_rk;
	keygen.create_relin_keys(he_rk);
	GaloisKeys he_gk;
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

	GIFT64_ENC tester(sk, remaining_level, loge, logn, logN - 1, total_level, scale, boundary_K, boot_deg, scale_factor, inverse_deg,
     context, keygen, encoder, encryptor, decryptor, evaluator, he_rk, he_gk);
    auto  start  =  std::chrono::high_resolution_clock::now();
    vector<Ciphertext> ptxt = tester.debug_test(sk, 128);
    // vector<Ciphertext> ptxt = tester.HE_decrypt(sk, 128);
    auto  end  =  std::chrono::high_resolution_clock::now();
    auto  duration  =  std::chrono::duration_cast<std::chrono::milliseconds>(end  -  start);
     //  输出结果
    std::cout  <<  "代码执行时间："  <<  duration.count()/1000  <<  "  秒 :: " << duration.count()%1000<< "毫秒"  <<  std::endl;
    std::cout  <<  "reEnc 代码执行时间："  <<  tester.time_BTSreEnc.count()/1000  <<  "  秒 :: " << tester.time_BTSreEnc.count()%1000<< "毫秒"  <<  std::endl;
    cout << "scale: " << ptxt[0].scale() << endl;
    // for(int i=0;i<8;i++){
    //     tester.debugPrint(ptxt[i], "first sbox");
    // }
    return 0;
}
