package main

import (
	"fmt"
	"math"
	"runtime"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_cipher"
	"github.com/ldsec/lattigo/v2/ckks_fv"
)
func main() {
	runtime.GOMAXPROCS(64)
	var err error

	var kgen ckks_fv.KeyGenerator
	var encoder ckks_fv.CKKSEncoder
	var sk *ckks_fv.SecretKey
	var pk *ckks_fv.PublicKey
	var encryptor ckks_fv.CKKSEncryptor
	var decryptor ckks_fv.CKKSDecryptor

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

	rotations := kgen.GenRotationIndexesForBootstrapping(params.LogSlots(), btpParams)

	rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)
	btpKey := ckks_fv.BootstrappingKey{Rlk: rlk, Rtks: rotkeys}
	
	symmetricKey := make([]uint8, 16)
	iv := make([]uint8, 16)
	aes, _ := ckks_cipher.NewAESCtr(symmetricKey, params, btpParams, btpKey, &encoder, &encryptor, &decryptor, iv) 
	Start1 := time.Now()
	ctxt, _ := aes.DebugTest(symmetricKey, 128 )
	aes.DebugPrint(ctxt[0], "after debug")
	end1 := time.Since(Start1)
	fmt.Printf("The bootstrap operation took %v\n", end1)
}

