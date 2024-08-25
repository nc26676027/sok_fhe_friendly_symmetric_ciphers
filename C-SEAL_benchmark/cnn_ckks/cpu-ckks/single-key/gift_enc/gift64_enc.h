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
#include "seal_helper.h"
#include <bitset>
#include <omp.h>
#include <NTL/RR.h>
 
using namespace std;
using namespace seal;

namespace GIFT64{

//debug mode, input all zero encrypted data
#define all_zero_in 1

constexpr unsigned blocksize = 64;  // Block size in bits
constexpr unsigned keysize = 128;  // Key size in bits
constexpr unsigned rounds = 28;                // Number of rounds

typedef std::bitset<blocksize> block;  // Store messages and states
typedef std::bitset<keysize> keyblock;

block ctr(const block& iv, uint64_t ctr);

class GIFT64_ENC: public SEALHelper {
protected:
    vector<vector<int> > bit_index;
    vector<vector<int> > bit_sbox;

	// A permutation with 16 entry, which is a function y=S[x] 
	uint8_t SBox[16] = {
		0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe
	};

	//A bit permutation with size 128
	size_t BitPerm_64[64] = {
		0, 17, 34, 51, 48, 1, 18, 35, 32, 49, 2, 19, 16, 33, 50, 3,
		4, 21, 38, 55, 52, 5, 22, 39, 36, 53, 6, 23, 20, 37, 54, 7,
		8, 25, 42, 59, 56, 9, 26, 43, 40, 57, 10, 27, 24, 41, 58, 11,
		12, 29, 46, 63, 60, 13, 30, 47, 44, 61, 14, 31, 28, 45, 62, 15
	};

	//A bit permutation with size 128
	size_t BitPerm_128[128] = {
		0, 33, 66, 99, 96, 1, 34, 67, 64, 97, 2, 35, 32, 65, 98, 3,
		4, 37, 70, 103, 100, 5, 38, 71, 68, 101, 6, 39, 36, 69, 102, 7,
		8, 41, 74, 107, 104, 9, 42, 75, 72, 105, 10, 43, 40, 73, 106, 11,
		12, 45, 78, 111, 108, 13, 46, 79, 76, 109, 14, 47, 44, 77, 110, 15,
		16, 49, 82, 115, 112, 17, 50, 83, 80, 113, 18, 51, 48, 81, 114, 19,
		20, 53, 86, 119, 116, 21, 54, 87, 84, 117, 22, 55, 52, 85, 118, 23,
		24, 57, 90, 123, 120, 25, 58, 91, 88, 121, 26, 59, 56, 89, 122, 27,
		28, 61, 94, 127, 124, 29, 62, 95, 92, 125, 30, 63, 60, 93, 126, 31
	};

    //GIFT utils
	Ciphertext IsEq(vector<Ciphertext> x, int index);
	void subbyte_lut(vector<Ciphertext> &x);
	void bit_permutation(vector<Ciphertext> &x);
    void add_round_key(vector<Ciphertext> &state, vector<Ciphertext> key);
    //TODO
    vector<vector<Ciphertext> > key_expansion(vector<int> MK);

public:
    GIFT64_ENC(std::vector<uint8_t> key, int remaining_level, long _loge, long _logn, long _logNh, long _L, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys);

    virtual ~GIFT64_ENC() = default;

    static chrono::microseconds time_BTSreEnc;
    void init_offline_data();
    void encrypt_key();
    void encrypt_input(block iv, size_t num_block);
	void encode_ciphertext( std::vector<uint8_t>& ciphertexts, size_t num_block );

    // add white key before first round (confirm!)
    void round_function(vector<Ciphertext> &state, vector<Ciphertext> round_key, size_t rnd);
    void last_round(vector<Ciphertext> &state, vector<Ciphertext> round_key);
    std::vector<Ciphertext> HE_decrypt(std::vector<uint8_t>& ciphertext, size_t bits);
    
	vector<Ciphertext> debug_test(std::vector<uint8_t>& ciphertexts, size_t bits);
};

}