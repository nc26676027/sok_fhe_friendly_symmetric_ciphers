#pragma	once
#include "seal/seal.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>
#include "Bootstrapper.h"
#include <omp.h>
#include <NTL/RR.h>
 
using namespace std;
using namespace seal;


class SEALHelper: public Bootstrapper {
protected:
    uint64_t mod_degree;
    double scale;
    int total_level;
    int remaining_level;
    std::vector<uint8_t> secret_key;
    std::vector<Ciphertext> secret_key_encrypted;
    std::vector<Ciphertext> input_encrypted;
    std::vector<Ciphertext> encoded_ct;

public:
    SEALHelper(std::vector<uint8_t> key, long _loge, long _logn, long _logNh, long _L, long _rL, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys);

    SEALHelper(std::vector<uint8_t> key, long _loge, long _logn, long _logNh, long _L, long _rL, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys,
    bool initialed);

    virtual ~SEALHelper() = default;

    void addScalar(Ciphertext &x, double scalar);
    void subScalar(Ciphertext &x, double scalar);
    //boolean function constructed from Heaan using CKKS scheme

    void XOR_inplace(Ciphertext &x, Ciphertext y);
    Ciphertext XOR(Ciphertext x, Ciphertext y);

    void NOT(Ciphertext &x);

    void AND_inplace(Ciphertext &x, Ciphertext y);
    Ciphertext AND(Ciphertext x, Ciphertext y);

    void OR_inplace(Ciphertext &x, Ciphertext y);
    Ciphertext OR(Ciphertext x, Ciphertext y);

    void cleanTensor(Ciphertext& x);
    void vector_left_rotation_inplace(vector<Ciphertext> &v, int rot_num);
    void vector_right_rotation_inplace(vector<Ciphertext> &v, int rot_num);
    
    void bootstrap_print(Ciphertext &cipher);
    void bootstrap_cipher(Ciphertext &cipher);
    void print_vector(std::vector<double> vec, std::size_t print_size, int prec);
    void print_vector_trunc(std::vector<double> vec, std::size_t print_size, size_t trunc, int prec);
    void debugPrint(Ciphertext encrypted_result, string str);
    void print_parameters();
    int get_chain_index(Ciphertext encrypted_data);

    std::vector<Ciphertext> get_encoded_ct();

    Ciphertext BinaryTreeAdd(std::vector<Ciphertext> &vector);

    Ciphertext construct_number_from_bits(vector<Ciphertext> ciphers, int error);
};
