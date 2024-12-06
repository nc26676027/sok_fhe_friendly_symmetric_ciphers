package main

import (
	"flag"
	"fmt"
	"runtime"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/ckks_cipher"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")


func main() {

	flag.Parse()
	runtime.GOMAXPROCS(64)
	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 14

	if *flagShort {
		LogN -= 3
	}
	//==============================
	//=== 1) RESIDUAL PARAMETERS ===
	//==============================

	// First we must define the residual parameters.
	// The residual parameters are the parameters used outside of the bootstrapping circuit.
	// For this example, we have a LogN=16, logQ = 55 + 10*40 and logP = 3*61, so LogQP = 638.
	// With LogN=16, LogQP=638 and H=192, these parameters achieve well over 128-bit of security.
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,                                              // Log2 of the ring degree
		LogQ:            []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40}, // Log2 of the ciphertext prime moduli
		LogP:            []int{61, 61, 61, 61, 61},                                 // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: 40,                                                // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	//==========================================
	//=== 2) BOOTSTRAPPING PARAMETERSLITERAL ===
	//==========================================

	// The bootstrapping circuit use its own Parameters which will be automatically
	// instantiated given the residual parameters and the bootstrapping parameters.

	// !WARNING! The bootstrapping parameters are not ensure to be 128-bit secure, it is the
	// responsibility of the user to check that the meet the security requirement and tweak them if necessary.

	// Note that the default bootstrapping parameters use LogN=16 and a ternary secret with H=192 non-zero coefficients
	// which provides parameters which are at least 128-bit if their LogQP <= 1550.

	// For this first example, we do not specify any circuit specific optional field in the bootstrapping parameters literal.
	// Thus we expect the bootstrapping to give a average precision of 27.9 bits with H=192 (and 24.4 with H=N/2)
	// if the plaintext values are uniformly distributed in [-1, 1] for both the real and imaginary part.
	// See `he/float/bootstrapping/parameters_literal.go` for detailed information about the optional fields.
	btpParametersLit := bootstrapping.ParametersLiteral{
		// We specify LogN to ensure that both the residual parameters and the bootstrapping parameters
		// have the same LogN. This is not required, but we want it for this example.
		LogN: utils.Pointy(LogN),

		// In this example we need manually specify the number of auxiliary primes (i.e. #Pi) used by the
		// evaluation keys of the bootstrapping circuit, so that the size of LogQP  meets the security target.
		LogP: []int{61, 61, 61, 61, 61},

		// In this example we manually specify the bootstrapping parameters' secret distribution.
		// This is not necessary, but we ensure here that they are the same as the residual parameters.
		Xs: params.Xs(),
	}

	//===================================
	//=== 3) BOOTSTRAPPING PARAMETERS ===
	//===================================

	// Now that the residual parameters and the bootstrapping parameters literals are defined, we can instantiate
	// the bootstrapping parameters.
	// The instantiated bootstrapping parameters store their own ckks.Parameter, which are the parameters of the
	// ring used by the bootstrapping circuit.
	// The bootstrapping parameters are a wrapper of ckks.Parameters, with additional information.
	// They therefore has the same API as the ckks.Parameters and we can use this API to print some information.
	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}

	if *flagShort {
		// Corrects the message ratio Q0/|m(X)| to take into account the smaller number of slots and keep the same precision
		btpParams.Mod1ParametersLiteral.LogMessageRatio += 16 - params.LogN()
	}

	// We print some information about the residual parameters.
	fmt.Printf("Residual parameters: logN=%d, logSlots=%d, H=%d, sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.ResidualParameters.LogN(),
		btpParams.ResidualParameters.LogMaxSlots(),
		btpParams.ResidualParameters.XsHammingWeight(),
		btpParams.ResidualParameters.Xe(), params.LogQP(),
		btpParams.ResidualParameters.MaxLevel(),
		btpParams.ResidualParameters.LogDefaultScale())

	// And some information about the bootstrapping parameters.
	// We can notably check that the LogQP of the bootstrapping parameters is smaller than 1550, which ensures
	// 128-bit of security as explained above.
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	//===========================
	//=== 4) KEYGEN & ENCRYPT ===
	//===========================

	// Now that both the residual and bootstrapping parameters are instantiated, we can
	// instantiate the usual necessary object to encode, encrypt and decrypt.

	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//========================
	//=== 5) BOOTSTRAPPING ===
	//========================
	iv := make([]uint8, 16)
	symmetricKey := []byte{
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
    }
	aes, _ := ckks_cipher.NewAESCtr(symmetricKey, params, btpParams, evk, encoder, encryptor, decryptor, iv) 
	Start1 := time.Now()
	// ctxt, _ := aes.DebugTest(symmetricKey, 128 )
	ctxt := aes.HEDecrypt(symmetricKey, 128 )
	aes.DebugPrint(ctxt[0], "\n\n Finished")
	end1 := time.Since(Start1)
	fmt.Printf("The Whole Decryption took %v\n", end1)
}

