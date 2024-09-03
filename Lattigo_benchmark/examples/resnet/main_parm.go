package main

import (
	"fmt"
	"math"
	"runtime"
	"sync"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
	"github.com/ldsec/lattigo/v2/utils"
)



func main2() {
	runtime.GOMAXPROCS(64)
	var err error

	var btp *ckks_fv.Bootstrapper
	var kgen ckks_fv.KeyGenerator
	var encoder ckks_fv.CKKSEncoder
	var sk *ckks_fv.SecretKey
	var pk *ckks_fv.PublicKey
	var encryptor ckks_fv.CKKSEncryptor
	var decryptor ckks_fv.CKKSDecryptor
	var plaintext *ckks_fv.Plaintext

	// Bootstrapping parameters
	// Four sets of parameters (index 0 to 3) ensuring 128 bit of security
	// are available in github.com/ldsec/lattigo/v2/ckks/bootstrap_params
	// LogSlots is hardcoded to 15 in the parameters, but can be changed from 1 to 15.
	// When changing logSlots make sure that the number of levels allocated to CtS and StC is
	// smaller or equal to logSlots.

	btpParams := ckks_fv.DefaultBootstrapParams[0]
	
	params, err := btpParams.Params()
	if err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, h = %d, logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), btpParams.H, params.LogQP(), params.Levels(), math.Log2(params.Scale()), params.Sigma())

	// Scheme context and keys
	kgen = ckks_fv.NewKeyGenerator(params)

	sk, pk = kgen.GenKeyPairSparse(btpParams.H)

	encoder = ckks_fv.NewCKKSEncoder(params)
	decryptor = ckks_fv.NewCKKSDecryptor(params, sk)
	encryptor = ckks_fv.NewCKKSEncryptorFromPk(params, pk)

	fmt.Println()
	// fmt.Println("Generating bootstrapping keys...")
	// params.SetLogSlots(13)
	// btpParams.SetLogSlots(13)
	rotations := kgen.GenRotationIndexesForBootstrapping(params.LogSlots(), btpParams)

	rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)
	btpKey := ckks_fv.BootstrappingKey{Rlk: rlk, Rtks: rotkeys}
	

	// params.SetLogSlots(13)
	// btpParams.SetLogSlots(13)
	if btp, err = ckks_fv.NewBootstrapper(params, btpParams, btpKey); err != nil {
		panic(err)
	}
	fmt.Println("Done")
	// params.SetLogSlots(15)
	// btpParams.SetLogSlots(8)
	// Generate a random plaintext
	valuesWant := make([]complex128, params.Slots() )
	randomVec := make([]float64, params.Slots())
	for i:=0;i<params.Slots();i++{
		randomVec[i] = utils.RandFloat64(-1, 1)
	}
	for i:=0; i<params.Slots();i++ {
		valuesWant[i] = complex(randomVec[i], 0.0)
	}


	ctxtPool := make([]*ckks_fv.Ciphertext, 8)

	for i:=0;i<8;i++{
		plaintext = encoder.EncodeComplexNew(valuesWant, params.LogSlots())
		// Encrypt
		ciphertext1 := encryptor.EncryptNew(plaintext)
		ctxtPool[i] = ciphertext1
	}

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest1 := printDebug2(params, ctxtPool[0], valuesWant, decryptor, encoder)

	// Bootstrap the ciphertext (homomorphic re-encryption)
	// It takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.Scale
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.Scale) can be used at the expense of one level.
	Start := time.Now()
	fmt.Println()
	fmt.Println("Bootstrapping...")
	var wg sync.WaitGroup
	wg.Add(len(ctxtPool))
	
	for i := range ctxtPool {
		go func(i int) {
			defer wg.Done()	
			// btsCopy := btp.ShallowCopy()
			// Update the original state with the modified copy
			ctxtPool[i] = btp.MultByConstNew(ctxtPool[i], 1)
			// btsCopy.RescaleMany(ctxtPool[i], 1, ctxtPool[i])
		}(i)
	}
	wg.Wait()	
	fmt.Println("Done")

	elapsed := time.Since(Start) //
	fmt.Printf("The bootstrap operation took %v\n", elapsed)
	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrapp(ciphertext)")
	printDebug2(params, ctxtPool[0], valuesTest1, decryptor, encoder)
}

func printDebug2(params *ckks_fv.Parameters, ciphertext *ckks_fv.Ciphertext, valuesWant []complex128, decryptor ckks_fv.CKKSDecryptor, encoder ckks_fv.CKKSEncoder) (valuesTest []complex128) {

	valuesTest = encoder.DecodeComplex(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	
	fmt.Println("Full slots data: ( ")
	for i:=(1<<13)-20;i< (1<<13);i++{
		fmt.Println( i, " ,", real(valuesTest[i]) )
	}
	fmt.Println(" )")
	

	precStats := ckks_fv.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
