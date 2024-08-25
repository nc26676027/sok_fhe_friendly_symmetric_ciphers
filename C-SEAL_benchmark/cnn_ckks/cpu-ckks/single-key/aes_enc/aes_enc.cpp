// aes_enc.cpp
#include "aes_enc.h"
using namespace std;
namespace AES{

std::chrono::microseconds AES_ENC::time_BTSreEnc(0);

block ctr(const block& iv, uint64_t ctr) {
  block out = iv;
  for (size_t i = blocksize - 64; i < blocksize; i++) {
    out[i] = iv[i] ^ ((ctr >> (blocksize - i - 1)) & 1);
  }
  return out;
}

AES_ENC::AES_ENC(std::vector<uint8_t> key, int _remaining_level, long _loge, long _logn, long _logNh, long _L, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys)
    :   SEALHelper(key, _loge, _logn, _logNh, _L, _remaining_level, _final_scale, _boundary_K, _sin_cos_deg,
                   _scale_factor, _inverse_deg, _context, _keygen, _encoder, _encryptor, 
                   _decryptor, _evaluator, _relin_keys, _gal_keys)
{
    vector<bitset<8> > monomial_order;
    for (int i=0;i<8;i++){
        bitset<8> x = 1<<i;
        monomial_order.push_back(x);
    }
    sbox_monomial_order = layered_combine_bin(monomial_order);
    
    //initial bit sbox
    for(int i=0;i<256;i++){
        bitset<8> tmp = SBox[i];
        bit_sbox.push_back(tmp);
    }
}

AES_ENC::AES_ENC(std::vector<uint8_t> key, int _remaining_level, long _loge, long _logn, long _logNh, long _L, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys,
    bool initialed)
    :   SEALHelper(key, _loge, _logn, _logNh, _L, _remaining_level, _final_scale, _boundary_K, _sin_cos_deg,
                   _scale_factor, _inverse_deg, _context, _keygen, _encoder, _encryptor, 
                   _decryptor, _evaluator, _relin_keys, _gal_keys, initialed)
{
    vector<bitset<8> > monomial_order;
    for (int i=0;i<8;i++){
        bitset<8> x = 1<<i;
        monomial_order.push_back(x);
    }
    sbox_monomial_order = layered_combine_bin(monomial_order);
    
    //initial bit sbox
    for(int i=0;i<256;i++){
        bitset<8> tmp = SBox[i];
        bit_sbox.push_back(tmp);
    }
}

vector<Ciphertext> AES_ENC::debug_test(std::vector<uint8_t>& ciphertexts, size_t bits){
    if (all_zero_in){
        bits = ( encoder.slot_count() ) * blocksize;
    }
    size_t num_block = ceil((double)bits / blocksize);
    
    block iv = 0;
    encrypt_key();
    encrypt_input(iv, num_block);

    vector<Ciphertext> state = input_encrypted;
    for(int i=0;i<state.size();i++){
        for(int j=0;j<total_level - remaining_level+5;j++){
            evaluator.mod_switch_to_next_inplace(state[i]);
        }
    }

    // auto  start  =  std::chrono::high_resolution_clock::now();
    // cout << "Sin Chain index after modswitch is: "<<get_chain_index(state[0]);
    vector<Ciphertext> Sout( state.begin(), state.begin()+8 );

    aes_subbyte_lut(Sout);

    for (int i=0;i<8;i++){
        debugPrint(Sout[i], "after");
    }

    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //  //  输出结果
    // std::cout  <<  "代码执行时间："  <<  duration.count()/1000  <<  "  秒 :: " << duration.count()%1000<< "毫秒"  <<  std::endl;


    // throw std::invalid_argument("finished\n\n\n\n\n");

    Ciphertext tmp = construct_number_from_bits(Sout, 10);
    debugPrint(tmp, "after construct");

    // aes_subbyte_lut(Sbox);
    // bootstrap_cipher(Sbox[0]);

    // for (int i=0;i<1;i++){
    //     debugPrint(tmp, "after");
    // }
    // state[0] = Sbox[0];
        // add cipher
    // encode_ciphertext(ciphertexts, num_block);
    // for(int i=0;i<blocksize;i++)
    // {
    //     XOR_inplace(state[i], encoded_ct[i]);
    // }

    // cout << "Sin Chain index after modswitch is: "<<get_chain_index(state[0]);
    // cout << "scale: " << input_encrypted[0].scale() << endl << endl;
    return state;
}

void AES_ENC::encrypt_key() {
    secret_key_encrypted.clear();
    secret_key_encrypted.reserve(keysize);
    for (int i=0;i<keysize;i++)
    {
        int32_t bit = (secret_key[i / 8] >> (i % 8)) & 1;
        Ciphertext sk_bit;
        vector<double> sk_duplicated;
        for (int j=0;j<encoder.slot_count();j++)
        {
            sk_duplicated.push_back(bit);
        }
        Plaintext sk_plain;
        encoder.encode(sk_duplicated, scale, sk_plain);
        Ciphertext sk_bit_encrypted;
        encryptor.encrypt(sk_plain, sk_bit_encrypted);
        secret_key_encrypted.push_back( std::move(sk_bit_encrypted) );
    }
}

void AES_ENC::encrypt_input(block iv, size_t num_block) {
    input_encrypted.clear();
    input_encrypted.reserve(blocksize);

    vector<block> input_data;
    input_data.reserve(mod_degree/2);

    for(int i=0;i<(mod_degree/2);i++){
        if(1){
            input_data[i] = 0;
        }else{
            input_data[i] = ctr(iv, i+1);
        }
    }

    for (int i=0;i<blocksize;i++)
    {
        vector<double> state_batched;
        for (int j=0;j<encoder.slot_count();j++)
        {
            state_batched.push_back( input_data[j][i] ); //encode i-th bit of j-th slots
        }
        Plaintext state_batched_plain;
        encoder.encode(state_batched, scale, state_batched_plain);
        Ciphertext state_batched_encrypted;
        encryptor.encrypt(state_batched_plain, state_batched_encrypted);
        input_encrypted.push_back( std::move(state_batched_encrypted) );
    }
}

void AES_ENC::encode_ciphertext( std::vector<uint8_t>& ciphertexts, size_t num_block ) {
    encoded_ct.clear();
    encoded_ct.reserve(blocksize);
    if(num_block < mod_degree/2) cout<<"data is not full pack, fill with 0..."<<endl;
    vector<block> encrypted_data;
    encrypted_data.reserve(num_block);
    for(size_t i=0;i<(mod_degree/2 - num_block);i++){
        block tmp = 0;
        encrypted_data.push_back( tmp );
    }

    for(int i=0;i<num_block;i++){
        if(all_zero_in){
            encrypted_data[i] = 0;
        }else{
            for (size_t k = 0; k < blocksize && i * blocksize + k < num_block*blocksize; k++) {
                size_t ind = i * blocksize + k;
                int32_t bit = (ciphertexts[ind / 8] >> ( ind % 8 )) & 0x1;
                encrypted_data[i][k] = bit;
            }
        }
    }

    for (int i=0;i<blocksize;i++)
    {
        vector<double> data_batched;
        for (int j=0;j<encoder.slot_count();j++)
        {
            data_batched.push_back( encrypted_data[j][i] ); //encode i-th bit of j-th slots
        }
        Plaintext encrypted_data_plain;
        encoder.encode(data_batched, scale, encrypted_data_plain);  
        Ciphertext encrypted_data_ctxt;
        encryptor.encrypt( encrypted_data_plain, encrypted_data_ctxt); 
        encoded_ct.push_back( std::move(encrypted_data_ctxt) );
    }
}

std::vector<Ciphertext> AES_ENC::HE_decrypt(
    std::vector<uint8_t>& ciphertexts, size_t bits) {
    if (all_zero_in){
        bits = ( encoder.slot_count() ) * blocksize;
    }
    size_t num_block = ceil((double)bits / blocksize);

    block iv = 0;

    encrypt_key();
    encrypt_input(iv, num_block);
    for(int i=0;i<input_encrypted.size();i++){
        for(int j=0;j<total_level-remaining_level+5;j++){
            evaluator.mod_switch_to_next_inplace(input_encrypted[i]);
        }
    }
    auto  start_aes  =  std::chrono::high_resolution_clock::now();
//AES encryption**********************************************
    //debug aes_round_function
    vector<Ciphertext> state = add_white_key(input_encrypted, secret_key_encrypted);

    for(int i=1;i<10;i++){
        cout << "round iterator : "<<i<<endl;
        // aes_two_round_function(state, secret_key_encrypted);
        aes_round_function(state, secret_key_encrypted);
    }
    cout << "round iterator : last round"<<endl;
    // aes_last_two_round(state, secret_key_encrypted);
    aes_last_round(state, secret_key_encrypted);
//AES encryption**********************************************
    auto  end_aes  =  std::chrono::high_resolution_clock::now();
    auto  duration_aes  =  std::chrono::duration_cast<std::chrono::milliseconds>(end_aes  -  start_aes);
     //  输出结果
    std::cout  <<  "代码执行时间："  <<  duration_aes.count()/1000  <<  "  秒 :: " << duration_aes.count()%1000<< "毫秒"  <<  std::endl;

    // add cipher
    // encode_ciphertext(ciphertexts, num_block);
    for(size_t i=0;i<blocksize;i++){
        if( (i%8 == 1)||(i%8 == 2)||(i%8 == 4)||(i%8 == 5) ){
            NOT(state[i]);
        }
    }
    for (size_t i=0;i<8;i++){
        debugPrint(state[i], "After transcipher all zero: ");
    }

    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<blocksize;i++)
    {   
        XOR_inplace(state[i], encoded_ct[i]);
    }
    return state;
}

Ciphertext AES_ENC::coefficient_mult_monomial(vector<Ciphertext> mon, const int *coeff_arr, int pos ){
    if ( mon.size() != sbox_monomial_order.size() ){
        throw std::runtime_error("monomial size must equal to sbox_monomial_order!\n");
    }
    vector<Ciphertext> ctxt_to_add;
    for(int i=0;i<mon.size();i++){
        int ind = static_cast<int>( sbox_monomial_order[i].to_ulong() ) - 1;//const item is not included, totally 255 items
        int coeff = coeff_arr[ind];
        if (coeff == 0){
            continue;
        } // when coeff = 0, do nothing!
        else if (coeff > 0){
            for(int j=0;j<coeff;j++){
                ctxt_to_add.push_back( mon[i] );
            }
        }
        else{
            int loop = abs(coeff);
            evaluator.negate_inplace(mon[i]);
            for(int j=0;j<loop;j++){
                ctxt_to_add.push_back( mon[i] );
            }
        }
    }

    for(int i=1;i<ctxt_to_add.size();i++){
        evaluator.add_inplace_reduced_error(ctxt_to_add[0], ctxt_to_add[i]);
    }
    evaluator.add_const_inplace(ctxt_to_add[0], bit_sbox[0][pos]);
    return ctxt_to_add[0];    
}

void AES_ENC::aes_subbyte_lut(vector<Ciphertext> &x){
    // construct 8-bit val of the sbox
    assert( ("The input length of the Sbox is wrong (8bit)!!", x.size() == 8 ) ); 
    // vector<Ciphertext> tmp = x;
    vector<Ciphertext> sbox_monomials = layered_combine(x);

    x[0] = coefficient_mult_monomial( sbox_monomials, sbox_0, 0 );
    x[1] = coefficient_mult_monomial( sbox_monomials, sbox_1, 1 );
    x[2] = coefficient_mult_monomial( sbox_monomials, sbox_2, 2 );
    x[3] = coefficient_mult_monomial( sbox_monomials, sbox_3, 3 );
    x[4] = coefficient_mult_monomial( sbox_monomials, sbox_4, 4 );
    x[5] = coefficient_mult_monomial( sbox_monomials, sbox_5, 5 );
    x[6] = coefficient_mult_monomial( sbox_monomials, sbox_6, 6 );
    x[7] = coefficient_mult_monomial( sbox_monomials, sbox_7, 7 );
}

// TODO...
vector<vector<Ciphertext> > AES_ENC::key_expansion(vector<int> MK){
  vector<vector<Ciphertext> > subkey;
  vector<int> tmp = MK;


  return subkey;
}

void AES_ENC::column_gf_2pow8_mult(vector<Ciphertext> &x)
{
    vector<Ciphertext> y = x;
    for(int i=0;i<4;i++){
        y[0+8*i] = x[1+8*i];
        y[1+8*i] = x[2+8*i];
        y[2+8*i] = x[3+8*i];
        y[3+8*i] = XOR( x[4+8*i], x[0+8*i]);
        y[4+8*i] = XOR( x[5+8*i], x[0+8*i]);
        y[5+8*i] = x[6+8*i];
        y[7+8*i] = x[0+8*i];
        y[6+8*i] = XOR( x[7+8*i], x[0+8*i]);
    }
    x = y;
}

// ### AES Mixcolumn Operation, apply the simplified equation:
//  $D_0$ = $x\cdot(b_0 \oplus b_1) \oplus b_1 \oplus b_2 \oplus b_3$\
//  $D_1$ = $x\cdot(b_1 \oplus b_2) \oplus b_2 \oplus b_3 \oplus b_0$\
//  $D_2$ = $x\cdot(b_2 \oplus b_3) \oplus b_3 \oplus b_0 \oplus b_1$\
//  $D_3$ = $x\cdot(b_3 \oplus b_0) \oplus b_0 \oplus b_1 \oplus b_2$,\
//  where one element $\in GF(2^8)$ multi x can be simplified to:\
//  $(a_7, a_6, a_5, a_4, a_3, a_2, a_1, a_0)$\
//  =$( a_6, a_5, a_4, a_3 \oplus a_7, a_2 \oplus a_7, a_1, a_0 \oplus a_7, a_7 )$
void AES_ENC::mixcolumn(vector<Ciphertext> &x)
{   
    vector<Ciphertext> x0, x1, x2, x3;
    for (int i = 0; i < 128; i++) {
        int mod = i % 32;
        if (mod < 8) {
            x0.push_back(x[i]);
        } else if (mod >= 8 && mod < 16) {
            x1.push_back(x[i]);
        } else if (mod >= 16 && mod < 24) {
            x2.push_back(x[i]);
        } else {
            x3.push_back(x[i]);
        }
    }

    vector<Ciphertext> y0 = x0;
    vector<Ciphertext> y1 = x1;
    vector<Ciphertext> y2 = x2;
    vector<Ciphertext> y3 = x3;

    vector<Ciphertext> z0 = x0;
    vector<Ciphertext> z1 = x1;
    vector<Ciphertext> z2 = x2;
    vector<Ciphertext> z3 = x3;

    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<32;i++){
        y0[i] = XOR(x0[i], x1[i]);
        y1[i] = XOR(x1[i], x2[i]);
        y2[i] = XOR(x2[i], x3[i]);
        y3[i] = XOR(x3[i], x0[i]);
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<32;i++){
        z0[i] = XOR(y1[i], x3[i]);
        z1[i] = XOR(y2[i], x0[i]);
        z2[i] = XOR(y3[i], x1[i]);
        z3[i] = XOR(y0[i], x2[i]);
    }

    column_gf_2pow8_mult(y0);
    column_gf_2pow8_mult(y1);
    column_gf_2pow8_mult(y2);
    column_gf_2pow8_mult(y3);

    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<32;i++){
        XOR_inplace(z0[i], y0[i]);
        XOR_inplace(z1[i], y1[i]);
        XOR_inplace(z2[i], y2[i]);
        XOR_inplace(z3[i], y3[i]);  
    }
    // link all vectors
    z0.insert(z0.end(), z1.begin(), z1.end());
    z0.insert(z0.end(), z2.begin(), z2.end());
    z0.insert(z0.end(), z3.begin(), z3.end());
    x = z0; 
}

// Operate ShiftRow operation over whole state, where the byte order
// is refferd to FIPS 197:
//         x0, x4,  x8, x12
//         x1, x5,  x9, x13
//         x2, x6, x10, x14
//         x3, x7, x11, x15
void AES_ENC::shift_row(vector<Ciphertext> &x)
{
    //shiftRow 1, which is left shift 1 byte.
    for (int i=0;i<8;i++){
        //Row one, left shift 1 byte
        x[1*8+i] = x[5*8+i];
        x[5*8+i] = x[9*8+i];
        x[9*8+i] = x[13*8+i];
        x[13*8+i] = x[1*8+i];
        //Row two, left shift 2 bytes
        x[2*8+i] = x[10*8+i];
        x[6*8+i] = x[14*8+i];
        x[10*8+i] = x[2*8+i];
        x[14*8+i] = x[6*8+i];
        //Row three, left shift 3 bytes
        x[3*8+i] = x[15*8+i];
        x[7*8+i] = x[3*8+i];
        x[11*8+i] = x[7*8+i];
        x[15*8+i] = x[11*8+i];
    }
}

// add white key before first round (confirm!)
vector<Ciphertext> AES_ENC::add_white_key(vector<Ciphertext> pt, vector<Ciphertext> key)
{
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<pt.size();i++){
        XOR_inplace(pt[i], key[i]);
    }
    return pt;
}

void AES_ENC::add_round_key(vector<Ciphertext> &state, vector<Ciphertext> key)
{
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<state.size();i++){
        XOR_inplace(state[i], key[i]);
    }
}
 
void AES_ENC::aes_round_function(vector<Ciphertext> &state, vector<Ciphertext> round_key)
{
    cout << "Chain index before sbox: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //SubByte
    for(size_t i=0;i<state.size();i++){
        int curr_level = context.get_context_data(state[i].parms_id())->chain_index();
        while(curr_level > 6){
            evaluator.mod_switch_to_next_inplace(state[i]);
            curr_level -= 1;
        };
    }

    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*8, state.begin()+(i+1)*8 );
        aes_subbyte_lut( Sout );
        for(int j=0;j<8;j++){
            state[i*8+j] = Sout[j];
        }
    }

    cout << "Chain index before BTS: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    //raise chain index and clean Ciphertext
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<state.size();i++){
        bootstrap_cipher(state[i]);
        if(i == 0) debugPrint(state[i], "BTS precise: ");
        cleanTensor(state[i]);               
    }

    auto end = std::chrono::high_resolution_clock::now();
    time_BTSreEnc += std::chrono::duration_cast<std::chrono::milliseconds>(end  -  start);

    cout << "ShiftRow Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //ShiftRow
    shift_row( state ); 
    cout << "MixColumn Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    mixcolumn( state );
    cout << "AddRoundKey Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //Add Round Key
    add_round_key( state, round_key );
}

void AES_ENC::aes_two_round_function(vector<Ciphertext> &state, vector<Ciphertext> round_key)
{
    cout << "Chain index before sbox: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //SubByte
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*8, state.begin()+(i+1)*8 );
        aes_subbyte_lut( Sout );
        for(int j=0;j<8;j++){
            state[i*8+j] = Sout[j];
        }
    }
    cout << "ShiftRow Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //ShiftRow
    shift_row( state ); 
    cout << "MixColumn Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    mixcolumn( state );
    cout << "AddRoundKey Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //Add Round Key
    add_round_key( state, round_key );
    for(size_t i=0;i<state.size();i++){
        int curr_level = context.get_context_data(state[i].parms_id())->chain_index();
        while(curr_level > 6){
            evaluator.mod_switch_to_next_inplace(state[i]);
            curr_level -= 1;
        };
    }
    cout << "Chain index before sbox: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //SubByte
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*8, state.begin()+(i+1)*8 );
        aes_subbyte_lut( Sout );
        for(int j=0;j<8;j++){
            state[i*8+j] = Sout[j];
        }
    }
    cout << "Chain index before BTS: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //raise chain index and clean Ciphertext
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<state.size();i++){
        bootstrap_cipher(state[i]);
        if(i == 0) debugPrint(state[i], "BTS precise: ");
        cleanTensor(state[i]);               
    }
    cout << "ShiftRow Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //ShiftRow
    shift_row( state ); 
    cout << "MixColumn Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    mixcolumn( state );
    cout << "AddRoundKey Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //Add Round Key
    add_round_key( state, round_key );
}

void AES_ENC::aes_last_round(vector<Ciphertext> &state, vector<Ciphertext> round_key)
{
    cout << "Chain index before sbox: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //SubByte
    for(size_t i=0;i<state.size();i++){
        int curr_level = context.get_context_data(state[i].parms_id())->chain_index();
        while(curr_level > 6){
            evaluator.mod_switch_to_next_inplace(state[i]);
            curr_level -= 1;
        };
    }
    //SubByte
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*8, state.begin()+(i+1)*8 );
        aes_subbyte_lut(Sout);
        for(int j=0;j<8;j++){
            state[i*8+j] = Sout[j];
        }
    }
    //raise chain index and clean Ciphertext
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<state.size();i++){
        bootstrap_cipher(state[i]);

        int curr_level = context.get_context_data(state[i].parms_id())->chain_index();
        while(curr_level > 12){
            evaluator.mod_switch_to_next_inplace(state[i]);
            curr_level -= 1;
        };

        cleanTensor(state[i]);               
    }
    //ShiftRow
    shift_row( state );
    //Add Round Key
    add_round_key( state, round_key );
}

void AES_ENC::aes_last_two_round(vector<Ciphertext> &state, vector<Ciphertext> round_key)
{
    cout << "Chain index before sbox: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //SubByte
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*8, state.begin()+(i+1)*8 );
        aes_subbyte_lut( Sout );
        for(int j=0;j<8;j++){
            state[i*8+j] = Sout[j];
        }
    }
    cout << "ShiftRow Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //ShiftRow
    shift_row( state ); 
    cout << "MixColumn Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    mixcolumn( state );
    cout << "AddRoundKey Chain: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //Add Round Key
    add_round_key( state, round_key );
    cout << "Chain index before sbox: "<<get_chain_index(state[0])<<", scale: "<<log2(state[0].scale())<<endl;
    //SubByte
    for(size_t i=0;i<state.size();i++){
        int curr_level = context.get_context_data(state[i].parms_id())->chain_index();
        while(curr_level > 6){
            evaluator.mod_switch_to_next_inplace(state[i]);
            curr_level -= 1;
        };
    }
    //SubByte
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*8, state.begin()+(i+1)*8 );
        aes_subbyte_lut(Sout);
        for(int j=0;j<8;j++){
            state[i*8+j] = Sout[j];
        }
    }
    //raise chain index and clean Ciphertext
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<state.size();i++){
        bootstrap_cipher(state[i]);
        cleanTensor(state[i]);               
    }
    //ShiftRow
    shift_row( state );
    //Add Round Key
    add_round_key( state, round_key );
}

}

