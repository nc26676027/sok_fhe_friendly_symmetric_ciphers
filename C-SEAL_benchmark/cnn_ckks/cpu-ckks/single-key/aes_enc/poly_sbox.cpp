#include "poly_sbox.h"
#include <unordered_map>
#include <bitset>
#include <assert.h>
#include "aes_enc.h"

using namespace std;

vector<bitset<8> > combine_bin(vector<bitset<8> > var1, vector<bitset<8> > var2){
    vector<bitset<8> > result;
    for(int i=0;i<var1.size();i++){
        for(int j=0;j<var2.size();j++){
            result.push_back( var1[i] ^ var2[j] );
        }
    }
    return result;
}

vector<bitset<8> > layered_combine_bin( vector<bitset<8> > variables ){
    if( (variables.size()%2) != 0 ){
        throw std::runtime_error("The input must be even");
    }
    vector<bitset<8> > result;
    if (variables.size() == 2){
        result.push_back(variables[0]);
        result.push_back(variables[1]);
        result.push_back(variables[0] ^ variables[1]);
        return result;
    }

    size_t mid = variables.size() / 2;
    vector<bitset<8> > left( variables.begin(), variables.begin()+mid );
    vector<bitset<8> > right( variables.begin()+mid, variables.end() );

    left = layered_combine_bin(left);
    right = layered_combine_bin(right);
    vector<bitset<8> > combined = combine_bin(left, right);
    result.insert( result.end(), combined.begin(), combined.end() );
    result.insert( result.end(), left.begin(), left.end() );
    result.insert( result.end(), right.begin(), right.end() );
    return result;
}



namespace AES{

//SEAL Ciphertext operations 
vector<Ciphertext> AES_ENC::combine(vector<Ciphertext> var1, vector<Ciphertext> var2){
    vector<Ciphertext> result;
    for(int i=0;i<var1.size();i++){
        for(int j=0;j<var2.size();j++){
            Ciphertext tmp = var1[i];
            evaluator.multiply_inplace_reduced_error(tmp, var2[j], relin_keys);
            evaluator.rescale_to_next_inplace(tmp);
            result.push_back( tmp );
        }
    }
    return result;
}

vector<Ciphertext> AES_ENC::layered_combine( vector<Ciphertext> variables ){
    if( (variables.size()%2) != 0 ){
        throw std::runtime_error("The input must be even");
    }
    vector<Ciphertext> result;
    if (variables.size() == 2){
        result.push_back(variables[0]);
        result.push_back(variables[1]);
        Ciphertext tmp = variables[0];
        evaluator.multiply_inplace_reduced_error(tmp, variables[1], relin_keys);
        evaluator.rescale_to_next_inplace(tmp);
        result.push_back( tmp );
        return result;
    }
    int mid = variables.size() / 2;
    vector<Ciphertext> left( variables.begin(), variables.begin()+mid );
    vector<Ciphertext> right( variables.begin()+mid, variables.end() );
    left = layered_combine(left);
    right = layered_combine(right);
    vector<Ciphertext> combined = combine(left, right);
    result.insert( result.end(), combined.begin(), combined.end() );
    result.insert( result.end(), left.begin(), left.end() );
    result.insert( result.end(), right.begin(), right.end() );

    return result;
}

}