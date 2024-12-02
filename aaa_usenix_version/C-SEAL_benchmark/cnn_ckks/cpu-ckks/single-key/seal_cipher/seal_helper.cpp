// encryption.cpp
#include "seal_helper.h"

SEALHelper::SEALHelper(std::vector<uint8_t> key, long _loge, long _logn, long _logNh, long _L, long _rL, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys)
    :   secret_key(key),
        Bootstrapper(_loge, _logn, _logNh, _L, _final_scale, _boundary_K, _sin_cos_deg,
                   _scale_factor, _inverse_deg, _context, _keygen, _encoder, _encryptor, 
                   _decryptor, _evaluator, _relin_keys, _gal_keys)
{
    scale = _final_scale;
    mod_degree = _context.first_context_data()->parms().poly_modulus_degree();
    total_level = _L;
    remaining_level = _rL;
    // scale_evaluator(context, encoder, relin_keys);

    // bootstrapping preprocessing
    cout << "Generating Optimal Minimax Polynomials..." << endl;
    prepare_mod_polynomial();

    cout << "Adding Bootstrapping Keys..." << endl;
    // addBootKeys_hoisting(gal_keys);
    addBootKeys_3(gal_keys);

    for(int i=0;i<slot_vec.size();i++){
        cout << "slot_vec: "<<slot_vec[i]<<endl;
    }
    cout << "Generating Linear Transformation Coefficients..." << endl;
    // generate_LT_coefficient();
    generate_LT_coefficient_3();
    cout << "after generate_LT_coefficient"<<endl;
}

SEALHelper::SEALHelper(std::vector<uint8_t> key, long _loge, long _logn, long _logNh, long _L, long _rL, double _final_scale, long _boundary_K, long _sin_cos_deg, long _scale_factor, long _inverse_deg,
	SEALContext &_context, 
    KeyGenerator &_keygen, 
    CKKSEncoder &_encoder, 
    Encryptor &_encryptor, 
    Decryptor &_decryptor, 
    Evaluator &_evaluator, 
    RelinKeys &_relin_keys, 
    GaloisKeys &_gal_keys,
    bool initialed)
    :   secret_key(key),
        Bootstrapper(_loge, _logn, _logNh, _L, _final_scale, _boundary_K, _sin_cos_deg,
                _scale_factor, _inverse_deg, _context, _keygen, _encoder, _encryptor, 
                _decryptor, _evaluator, _relin_keys, _gal_keys)
{
    if(initialed == false) cout <<"not initial version"<<endl;
    scale = _final_scale;
    mod_degree = _context.first_context_data()->parms().poly_modulus_degree();
    total_level = _L;
    remaining_level = _rL;
    // scale_evaluator(context, encoder, relin_keys);
}

void SEALHelper::bootstrap_print(Ciphertext &cipher)
{
    cout << "bootstrapping..." << endl;
	Ciphertext ctxt;
	chrono::high_resolution_clock::time_point time_start, time_end;
	chrono::microseconds time_diff;

	ctxt = cipher;
	time_start = chrono::high_resolution_clock::now();
	bootstrap_hoisting(cipher, ctxt);
	// bootstrapper.bootstrap_real_3(rtn, ctxt);
	time_end = chrono::high_resolution_clock::now();
	time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
	cout << "time : " << time_diff.count() / 1000 << " ms" << endl;
}

void SEALHelper::bootstrap_cipher(Ciphertext &cipher)
{
    int reEnc = 0;

    if(reEnc == 0){
        int curr_level = get_chain_index(cipher);
        if (curr_level >= 3) {
            cipher.scale() = scale;
            slim_bootstrap_inplace(cipher);
        } else {
            throw std::runtime_error("Condition failed, must to left 3 level for slim-BTS...");
        }
    }
    else if (reEnc == 1){
        Plaintext plain_result;
        decryptor.decrypt(cipher, plain_result);
        vector<double> result;
        encoder.decode(plain_result, result);
        encoder.encode(result, scale, plain_result);
        encryptor.encrypt(plain_result, cipher);
        for(int j=0;j<total_level - remaining_level;j++){
            evaluator.mod_switch_to_next_inplace(cipher);
        }     
    }
    else{
        int curr_level = context.get_context_data(cipher.parms_id())->chain_index();
        while(curr_level > 0){
            evaluator.mod_switch_to_next_inplace(cipher);
            curr_level -= 1;
        };
        bootstrap_inplace_real_3(cipher);
    }
}

void SEALHelper::print_parameters()
{
    auto &context_data = *context.key_context_data();
    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

void SEALHelper::addScalar(Ciphertext &x, double scalar){
    evaluator.add_const_inplace(x, scalar);
}

void SEALHelper::subScalar(Ciphertext &x, double scalar){
    evaluator.add_const_inplace(x, -scalar);
}

//boolean function constructed from Heaan using CKKS scheme
void SEALHelper::XOR_inplace(Ciphertext &x, Ciphertext y){
  // (x-y)^2
    evaluator.sub_inplace_reduced_error(x, y);
    evaluator.square_inplace(x);
    evaluator.relinearize_inplace(x, relin_keys);
    evaluator.rescale_to_next_inplace(x);
}
//boolean function constructed from Heaan using CKKS scheme
Ciphertext SEALHelper::XOR(Ciphertext x, Ciphertext y){
  // (x-y)^2
    evaluator.sub_inplace_reduced_error(x, y);
    evaluator.square_inplace(x);
    evaluator.relinearize_inplace(x, relin_keys);
    evaluator.rescale_to_next_inplace(x);
    return x;
}

void SEALHelper::NOT(Ciphertext &x){
  // (1-x)
    evaluator.negate_inplace(x);
    addScalar( x, 1.0 );
}

void SEALHelper::AND_inplace(Ciphertext &x, Ciphertext y){
    evaluator.multiply_inplace(x, y);
    evaluator.relinearize_inplace(x, relin_keys);
    evaluator.rescale_to_next_inplace(x);
}
Ciphertext SEALHelper::AND(Ciphertext x, Ciphertext y){
    evaluator.multiply_inplace(x, y);
    evaluator.relinearize_inplace(x, relin_keys);
    evaluator.rescale_to_next_inplace(x);
    return x;
}

void SEALHelper::OR_inplace(Ciphertext &x, Ciphertext y){
  //x + y - x\cdot y
    Ciphertext tmp = x;
    evaluator.multiply_inplace(tmp, y);
    evaluator.relinearize_inplace(tmp, relin_keys);
    evaluator.rescale_to_next_inplace(tmp);
    evaluator.add_inplace_reduced_error(x, y);    
    evaluator.sub_inplace_reduced_error(x, tmp);
}
Ciphertext SEALHelper::OR(Ciphertext x, Ciphertext y){
  //x + y - x\cdot y
    Ciphertext tmp = x;
    evaluator.multiply_inplace(tmp, y);
    evaluator.relinearize_inplace(tmp, relin_keys);
    evaluator.rescale_to_next_inplace(tmp);
    evaluator.add_inplace_reduced_error(x, y);    
    evaluator.sub_inplace_reduced_error(x, tmp);
    return x;
}

void SEALHelper::cleanTensor(Ciphertext& x){
    Ciphertext square = x;
    evaluator.square_inplace(square);
    evaluator.relinearize_inplace(square, relin_keys);
    evaluator.rescale_to_next_inplace(square);
    Ciphertext cube = square;
    evaluator.multiply_inplace_reduced_error(cube, x, relin_keys);
    evaluator.relinearize_inplace(cube, relin_keys);
    evaluator.rescale_to_next_inplace(cube);
    // computation of 3x^2-2x^3 which is a rough-but-good-enough approximation of
    // the sign fucntion 
    x = square;
    evaluator.add_inplace_reduced_error(x, square);
    evaluator.add_inplace_reduced_error(x, square);
    evaluator.sub_inplace_reduced_error(x, cube);
    evaluator.sub_inplace_reduced_error(x, cube);
}

void SEALHelper::vector_left_rotation_inplace(vector<Ciphertext> &v, int rot_num)
{
  for (int i=0;i<rot_num;i++){
    Ciphertext first = v.front();
    v.insert(v.end(), first);
    v.erase(v.begin());
  }
}

void SEALHelper::vector_right_rotation_inplace(vector<Ciphertext> &v, int rot_num)
{
  for (int i=0;i<rot_num;i++){
    Ciphertext last = v.back();
    v.insert(v.begin(), last );
    v.erase(v.end());
  }
}

std::vector<Ciphertext> SEALHelper::get_encoded_ct(){
    return encoded_ct;
}

/*
Helper function: Prints a vector of floating-point values.
*/
void SEALHelper::print_vector(std::vector<double> vec, std::size_t print_size = 4, int prec = 3)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if (slot_count <= 2 * print_size)
    {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        vec.resize(std::max(vec.size(), 2 * print_size));
        std::cout << "    [";
        for (std::size_t i = 0; i < print_size; i++)
        {
            std::cout << " " << vec[i] << ",";
        }
        if (vec.size() > 2 * print_size)
        {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}
void SEALHelper::print_vector_trunc(std::vector<double> vec, std::size_t print_size, size_t trunc, int prec)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if (slot_count <= 2 * print_size)
    {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        // vec.resize(std::max(vec.size(), 2 * print_size));

        std::cout << "    [";
        for(size_t i = 0;i<vec.size();i++){
            if ( (i%trunc) < print_size)
            {
                std::cout << " " << vec[i] << ",";
            }
            if( (i%trunc) == 0){
                std::cout << endl;
            }
        }

        if (vec.size() > 2 * print_size)
        {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}

void SEALHelper::debugPrint(Ciphertext encrypted_result, string str = ""){
    /*
    Decrypt, decode, and print the result.
    */
    cout << "Chain Index = "<<get_chain_index(encrypted_result);
    cout << "Scale: "<<encrypted_result.scale()<<endl;
    Plaintext plain_result;
    decryptor.decrypt(encrypted_result, plain_result);
    vector<double> result;
    encoder.decode(plain_result, result);
    cout << str << endl;
    // if(str.compare("image0:") == 0){
    //     print_vector_trunc(result, 100, 4096, 7); 
    // }else{
    //     print_vector(result, 7, 3);
    // }
    print_vector(result, 7, 4);
}

int SEALHelper::get_chain_index(Ciphertext encrypted_data){
	return context.get_context_data(encrypted_data.parms_id())->chain_index();
}

Ciphertext SEALHelper::BinaryTreeAdd(std::vector<Ciphertext> &vector) {

	for(size_t j = 1; j < vector.size(); j=j*2) {
		for(size_t i = 0; i<vector.size(); i = i + 2*j) {
			if ((i+j)<vector.size())
				evaluator.add_inplace_reduced_error(vector[i],vector[i+j]);
		}
	}
	return vector[0];
}

Ciphertext SEALHelper::construct_number_from_bits(vector<Ciphertext> ciphers, int error){
    Ciphertext out = ciphers[0];
    evaluator.add_inplace_reduced_error(out, ciphers[1]);
    evaluator.add_inplace_reduced_error(out, ciphers[1]);

    double error_bound = 0.5 * error;
    size_t n = ciphers.size();
    vector<Ciphertext> var2add;
    for(size_t i=2;i<n;i++){
        if (i < error_bound){
            Ciphertext b;
            double tmp = pow(2, i/2);
            evaluator.multiply_const(ciphers[i], tmp, b);
            evaluator.rescale_to_next_inplace(b);
            evaluator.square_inplace(b);
            evaluator.relinearize_inplace(b, relin_keys);
            evaluator.rescale_to_next_inplace(b);
            if( i%2 == 1 ){
                evaluator.add_inplace_reduced_error(b, b);
            }
            var2add.push_back(b);
        }
        else{
            Ciphertext b;
            double tmp = pow(2, i/4);
            evaluator.multiply_const(ciphers[i], tmp, b);
            evaluator.rescale_to_next_inplace(b);

            evaluator.square_inplace(b);
            evaluator.relinearize_inplace(b, relin_keys);
            evaluator.rescale_to_next_inplace(b);

            evaluator.square_inplace(b);
            evaluator.relinearize_inplace(b, relin_keys);
            evaluator.rescale_to_next_inplace(b);
            
            int copy_num = int( pow(2, i%4) );
            for(int j=0;j<copy_num;j++){
                var2add.push_back(b);
            }
        }  
    }
    Ciphertext sum_from_2_to_n = BinaryTreeAdd(var2add);
    evaluator.add_inplace_reduced_error(out, sum_from_2_to_n);
    return out;
}
