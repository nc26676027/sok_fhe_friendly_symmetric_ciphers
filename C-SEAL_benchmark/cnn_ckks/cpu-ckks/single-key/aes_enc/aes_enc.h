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
#include "poly_sbox.h"
#include <bitset>
#include <omp.h>
#include <NTL/RR.h>
 
using namespace std;
using namespace seal;
// 确保头文件只被包含一次
#ifndef AES_ENC_H
#define AES_ENC_H

namespace AES{

//debug mode, input all zero encrypted data
#define all_zero_in 0

constexpr unsigned blocksize = 128;  // Block size in bits
constexpr unsigned keysize = 128;  // Key size in bits
constexpr unsigned rounds = 10;                // Number of rounds

typedef std::bitset<blocksize> block;  // Store messages and states
typedef std::bitset<keysize> keyblock;

block ctr(const block& iv, uint64_t ctr);

// A permutation with 256 entry, which is a function y=S[x] 
static uint8_t SBox[256] = {
	0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5,
	0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
	0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0,
	0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
	0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC,
	0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
	0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A,
	0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
	0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0,
	0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
	0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B,
	0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
	0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85,
	0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
	0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5,
	0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
	0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17,
	0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
	0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88,
	0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
	0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C,
	0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
	0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9,
	0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
	0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6,
	0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
	0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E,
	0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
	0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94,
	0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
	0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68,
	0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
};



class AES_ENC: public SEALHelper {
protected:
    vector<vector<int> > bit_index;
    vector<bitset<8> > bit_sbox;
    vector<bitset<8> > sbox_monomial_order;

    //AES utils
	Ciphertext IsEq(vector<Ciphertext> x, int index);
	void aes_subbyte_lut(vector<Ciphertext> &x);
    void gf_2pow8_mult(vector<Ciphertext> &x);
    void column_gf_2pow8_mult(vector<Ciphertext> &x);
    // ### AES Mixcolumn Operation, apply the simplified equation:
    //  $D_0$ = $x\cdot(b_0 \oplus b_1) \oplus b_1 \oplus b_2 \oplus b_3$\
    //  $D_1$ = $x\cdot(b_1 \oplus b_2) \oplus b_2 \oplus b_3 \oplus b_0$\
    //  $D_2$ = $x\cdot(b_2 \oplus b_3) \oplus b_3 \oplus b_0 \oplus b_1$\
    //  $D_3$ = $x\cdot(b_3 \oplus b_0) \oplus b_0 \oplus b_1 \oplus b_2$,\
    //  where one element $\in GF(2^8)$ multi x can be simplified to:\
    //  $(a_7, a_6, a_5, a_4, a_3, a_2, a_1, a_0)$\
    //  =$( a_6, a_5, a_4, a_3 \oplus a_7, a_2 \oplus a_7, a_1, a_0 \oplus a_7, a_7 )$
    void mixcolumn(vector<Ciphertext> &x);
    // Operate ShiftRow operation over whole state, where the byte order
    // is refferd to FIPS 197:
    //         x0, x4,  x8, x12
    //         x1, x5,  x9, x13
    //         x2, x6, x10, x14
    //         x3, x7, x11, x15
    void shift_row(vector<Ciphertext> &x);
    void add_round_key(vector<Ciphertext> &state, vector<Ciphertext> key);
	vector<Ciphertext> add_white_key(vector<Ciphertext> state, vector<Ciphertext> key);
    //TODO
    vector<vector<Ciphertext> > key_expansion(vector<int> MK);

public:
    AES_ENC(std::vector<uint8_t> key, int remaining_level, long _loge, long _logn, long _logNh, long _L, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys);

    AES_ENC(std::vector<uint8_t> key, int remaining_level, long _loge, long _logn, long _logNh, long _L, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys,
    bool initialed);

    virtual ~AES_ENC() = default;

    static chrono::microseconds time_BTSreEnc;
    //sbox polynomial generation
    vector<Ciphertext> combine(vector<Ciphertext> var1, vector<Ciphertext> var2);
    vector<Ciphertext> layered_combine( vector<Ciphertext> variables );
    Ciphertext coefficient_mult_monomial(vector<Ciphertext> mon, const int *coeff_arr, int pos);

    void encrypt_key();
    void encrypt_input(block iv, size_t num_block);
	void encode_ciphertext( std::vector<uint8_t>& ciphertexts, size_t num_block );

    // add white key before first round (confirm!)
    void aes_round_function(vector<Ciphertext> &state, vector<Ciphertext> round_key);
    void aes_two_round_function(vector<Ciphertext> &state, vector<Ciphertext> round_key);

    void aes_last_round(vector<Ciphertext> &state, vector<Ciphertext> round_key);
    void aes_last_two_round(vector<Ciphertext> &state, vector<Ciphertext> round_key);

    std::vector<Ciphertext> HE_decrypt(std::vector<uint8_t>& ciphertext, size_t bits);
    
	vector<Ciphertext> debug_test(std::vector<uint8_t>& ciphertexts, size_t bits);
};

}
#endif