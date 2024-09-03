// gift64_enc.cpp
#include "gift64_enc.h"

namespace GIFT64{

std::chrono::microseconds GIFT64_ENC::time_BTSreEnc(0);

block ctr(const block& iv, uint64_t ctr) {
  block out = iv;
  for (size_t i = blocksize - 64; i < blocksize; i++) {
    out[i] = iv[i] ^ ((ctr >> (blocksize - i - 1)) & 1);
  }
  return out;
}

GIFT64_ENC::GIFT64_ENC(std::vector<uint8_t> key, int _remaining_level, long _loge, long _logn, long _logNh, long _L, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
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
    init_offline_data();
}

vector<Ciphertext> GIFT64_ENC::debug_test(std::vector<uint8_t>& ciphertexts, size_t bits){
    if (all_zero_in){
        bits = ( encoder.slot_count() ) * blocksize;
    }
    size_t num_block = ceil((double)bits / blocksize);
    
    block iv = 0;
    encrypt_key();
    encrypt_input(iv, num_block);

    vector<Ciphertext> state = input_encrypted;
    for(int i=0;i<state.size();i++){
        for(int j=0;j<total_level - remaining_level;j++){
            evaluator.mod_switch_to_next_inplace(state[i]);
        }
    }

    auto  start  =  std::chrono::high_resolution_clock::now();
    cout << "Sin Chain index after modswitch is: "<<get_chain_index(state[0]);
    cout << "scale: " << state[0].scale() << endl;
    vector<Ciphertext> Sin(state.begin(), state.begin()+4);
    for(int i=0;i<4;i++){
        debugPrint(Sin[i], " before hoist sbox");
    }
    // bootstrap_cipher(state[0]);
    subbyte_lut(Sin);
    for(int i=0;i<4;i++){
        debugPrint(Sin[i], " before hoist sbox");
    }
    
    // cout << "Sin Chain index after modswitch is: "<<get_chain_index(state[0]);
    // cout << "scale: " << input_encrypted[0].scale() << endl << endl;
    auto  end  =  std::chrono::high_resolution_clock::now();
    auto  duration  =  std::chrono::duration_cast<std::chrono::milliseconds>(end  -  start);
    std::cout  <<  "code running time"  <<  duration.count()/1000  <<  "  s :: " << duration.count()%1000<< "milisecond"  <<  std::endl;
    return state;
}

void GIFT64_ENC::encrypt_key() {
    secret_key_encrypted.clear();
    secret_key_encrypted.reserve(keysize);
    for (int i=0;i<keysize;i++)
    {
        int32_t bit = (secret_key[i / 8] >> (7 - i % 8)) & 1;
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

void GIFT64_ENC::encrypt_input(block iv, size_t num_block) {
    input_encrypted.clear();
    input_encrypted.reserve(blocksize);

    vector<block> input_data;
    input_data.reserve(num_block);
    for(int i=0;i<num_block;i++){
        if(all_zero_in){
            input_data[i] = 1;
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

void GIFT64_ENC::encode_ciphertext( std::vector<uint8_t>& ciphertexts, size_t num_block ) {
    encoded_ct.clear();
    encoded_ct.reserve(blocksize);

    vector<block> encrypted_data;
    encrypted_data.reserve(num_block);

    for(int i=0;i<num_block;i++){
        if(all_zero_in){
            encrypted_data[i] = 0;
        }else{
            for (size_t k = 0; k < blocksize && i * blocksize + k < num_block*blocksize; k++) {
                size_t ind = i * blocksize + k;
                int32_t bit = (ciphertexts[ind / 8] >> (7 - ind % 8)) & 1;
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

void GIFT64_ENC::init_offline_data(){
     //initial bit_index, from 0 to 255 in the binary form
    for (size_t i=0;i<16;i++){
        vector<int> index;
        int tmp = i;
        for(size_t j=0;j<4;j++){
            index.push_back( tmp&1 );
            tmp /= 2;
        }
        reverse( index.begin(), index.end() );
        bit_index.push_back(index);
    }

    //initial bit_sbox
    for(size_t i=0;i<16;i++){
        vector<int> tmp;
        int val = SBox[i];
        for(size_t j=0;j<4;j++){
            tmp.push_back( val & 1);
            val /= 2;
        }
        reverse(tmp.begin(), tmp.end());
        bit_sbox.push_back(tmp);
    }
}

std::vector<Ciphertext> GIFT64_ENC::HE_decrypt(
    std::vector<uint8_t>& ciphertexts, size_t bits) {
    if (all_zero_in){
        bits = ( encoder.slot_count() ) * blocksize;
    }
    size_t num_block = ceil((double)bits / blocksize);
    
    block iv = 0;

    encrypt_key();
    encrypt_input(iv, num_block);
    for(int i=0;i<input_encrypted.size();i++){
        for(int j=0;j<total_level - remaining_level+3;j++){
            evaluator.mod_switch_to_next_inplace(input_encrypted[i]);
        }
    }
    auto  start_gift  =  std::chrono::high_resolution_clock::now();
// encryption**********************************************
    //debug round_function
    vector<Ciphertext> state =  input_encrypted;
    for(int i=1;i<4;i++){
        cout << "round iterator : "<<i<<endl;
        round_function(state, secret_key_encrypted, i);
    }
    cout << "round iterator : last round"<<endl;
    last_round(state, secret_key_encrypted);
// encryption**********************************************
    
    cout << "scale: " << state[0].scale() << endl;
    for(int i=0;i<8;i++){
        debugPrint(state[i], "first sbox");
    }  
    // tester.debug_test(sk, 128);
    auto  end_gift  =  std::chrono::high_resolution_clock::now();
    auto  duration_gift  =  std::chrono::duration_cast<std::chrono::milliseconds>(end_gift  -  start_gift);

    std::cout  <<  "code running time"  <<  duration_gift.count()/1000  <<  "  s :: " << duration_gift.count()%1000<< "milisecond"  <<  std::endl;

    // add cipher
    encode_ciphertext(ciphertexts, num_block);
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<blocksize;i++)
    {
        XOR_inplace(state[i], encoded_ct[i]);
    }
    return state;
}

// //compare two 8-bit data is equal or not, if eq return "1" else return "0"(two input: 1. cipher x 2. index from 0 to 255)
Ciphertext GIFT64_ENC::IsEq(vector<Ciphertext> x, int index){
    for (int i=0;i<4;i++){
        if( bit_index[index][i] == 0 ){ 
            NOT(x[i]);
        }
    }
    //multiplied togather with multi-depth log_2(4) = 2                            
    for(int i=0;i<2;i++){
        evaluator.multiply_inplace(x[i], x[i+2]);
        evaluator.relinearize_inplace(x[i], relin_keys);
        evaluator.rescale_to_next_inplace(x[i]);
    }
    evaluator.multiply_inplace(x[0], x[1]);
    evaluator.relinearize_inplace(x[0], relin_keys);
    evaluator.rescale_to_next_inplace(x[0]);
    return x[0];
}
//compare two 8-bit data is equal or not, if eq return "1" else return "0"(two input: 1. cipher x 2. index from 0 to 255)
// Ciphertext GIFT64_ENC::IsEq(vector<Ciphertext> x, int index){
//     vector<Ciphertext> multi_queue;
//     for (int i=0;i<4;i++){
//         if( bit_index[index][i] == 1 ){ 
//             multi_queue.push_back(x[i]);
//         }
//     }
//     Ciphertext tmp;

//     cout << "before multiply_many size: "<<multi_queue.size()<<endl;
//     if ( multi_queue.size() == 1 )  return multi_queue[0];

//     for(int i=1;i<multi_queue.size();i++){
//         evaluator.multiply_inplace(multi_queue[0], multi_queue[i]);
//         evaluator.relinearize_inplace(multi_queue[0], relin_keys);
//         evaluator.rescale_to_next_inplace(multi_queue[0]);
//         for(int j=i+1;j<multi_queue.size();j++)
//             evaluator.rescale_to_next_inplace(multi_queue[j]);
//     }

//     return multi_queue[0];
// }

void GIFT64_ENC::subbyte_lut(vector<Ciphertext> &x){
    // construct 4-bit val of the sbox
    assert( ("The input length of the Sbox is wrong (4bit)!!", x.size() == 4 ) ); 
    //perform IsEq(x, i),  parallel operation
    vector<Ciphertext> flag_vec(16, x[0]);
    // #pragma omp parallel for schedule(dynamic)
    for (int i=1;i<16;i++){
        flag_vec[i] = IsEq(x, i);
    }

    vector<vector<Ciphertext> > ctxt_to_add;
    for(int i=0;i<4;i++){
        vector<Ciphertext> bit_to_add;
        for(int j=0;j<16;j++){
            if( bit_sbox[j][i] == 1 ){
                bit_to_add.push_back(flag_vec[j]);
            }
        }
        ctxt_to_add.push_back(bit_to_add);
    }
    for(int i=0;i<4;i++){
        Ciphertext result;
        evaluator.add_many(ctxt_to_add[i], result);
        // for(int j=1;j<ctxt_to_add.size();j++){
        //     evaluator.add_inplace_reduced_error(ctxt_to_add[i][0], ctxt_to_add[i][j]);
        // }   
        x[i] = result;
    }
}

// TODO...
vector<vector<Ciphertext> > GIFT64_ENC::key_expansion(vector<int> MK){
  vector<vector<Ciphertext> > subkey;
  vector<int> tmp = MK;


  return subkey;
}

void GIFT64_ENC::bit_permutation(vector<Ciphertext> &x)
{   
    for (size_t i = 0; i < 64; i++) {
        x[ BitPerm_64[i] ] = x[i]; 
    }
}

void GIFT64_ENC::add_round_key(vector<Ciphertext> &state, vector<Ciphertext> key)
{
    // #pragma omp parallel for schedule(dynamic)
    for (size_t i=0;i<16;i++){
        XOR_inplace(state[4*i+1], key[i]);
        XOR_inplace(state[4*i+2], key[i+16]);
        evaluator.mod_switch_to_next_inplace(state[4*i]);
        evaluator.mod_switch_to_next_inplace(state[4*i+3]);
    }
}
 
void GIFT64_ENC::round_function(vector<Ciphertext> &state, vector<Ciphertext> round_key, size_t rnd)
{
    cout << "Chain index before sbox is  "<<get_chain_index(state[0])<<endl;
    //SubByte
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*4, state.begin()+(i+1)*4 );
        subbyte_lut( Sout );
        for(int j=0;j<4;j++){
            state[i*4+j] = Sout[j];
        }
    }
    cout << "Chain index after sbox is  "<<get_chain_index(state[0])<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    //raise chain index and clean Ciphertext
    if(rnd%3 == 0){
        #pragma omp parallel for schedule(dynamic)
        for (int i=0;i<state.size();i++){
            bootstrap_cipher(state[i]);
            if(i == 0)
                debugPrint(state[i], "BTS precise = ");
            cleanTensor(state[i]);               
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    time_BTSreEnc += std::chrono::duration_cast<std::chrono::milliseconds>(end  -  start);
    cout << "Chain index after bootcleaning is  "<<get_chain_index(state[0])<<endl;
    //Bit Permutation
    bit_permutation( state ); 
    //Add Round Key
    add_round_key( state, round_key );
}

void GIFT64_ENC::last_round(vector<Ciphertext> &state, vector<Ciphertext> round_key)
{
    //SubByte
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<16;i++){
        vector<Ciphertext> Sout( state.begin()+i*4, state.begin()+(i+1)*4 );
        subbyte_lut(Sout);
        for(int j=0;j<4;j++){
            state[i*4+j] = Sout[j];
        }
    }
    //Bit Perm
    bit_permutation( state );
    //Add Round Key
    add_round_key( state, round_key );
}

}

