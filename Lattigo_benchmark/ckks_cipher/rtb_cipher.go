package ckks_cipher

import (
	"fmt"
	"math"

	"github.com/ldsec/lattigo/v2/ckks_fv"
)

type RtBCipher struct {
    *ckks_fv.Bootstrapper
	params			*ckks_fv.Parameters
	encoder			*ckks_fv.CKKSEncoder
	encryptor		*ckks_fv.CKKSEncryptor
	decryptor		*ckks_fv.CKKSDecryptor
    totalLevel	 	int
    remainingLevel	int
    symmetricKey	[]uint8
    keyEncrypted	[]*ckks_fv.Ciphertext
    inputEncrypted	[]*ckks_fv.Ciphertext
    encodeCipher	[]*ckks_fv.Ciphertext
}

// NewRtBCipher create a new Real to boolean transcipher helper
func NewRtBCipher( key []uint8, params_ *ckks_fv.Parameters, btpParams_ *ckks_fv.BootstrappingParameters, btpKey_ ckks_fv.BootstrappingKey, encoder_ *ckks_fv.CKKSEncoder, encryptor_ *ckks_fv.CKKSEncryptor, decryptor_ *ckks_fv.CKKSDecryptor  ) (rtb *RtBCipher, err error) {
	btp, err := ckks_fv.NewBootstrapperReal(params_, btpParams_, btpKey_)
	if err!= nil {
		panic("btp in RtBCipher initial failed!")
	}
	rtb = &RtBCipher{
		Bootstrapper: 	btp,
		params: 		params_,
		encoder:  		encoder_,
		symmetricKey: 	key,
		encryptor:		encryptor_,
		decryptor:		decryptor_,
		totalLevel: 	len( params_.Qi() ) - 1,
		remainingLevel: len( btpParams_.ResidualModuli ) - 1,
		keyEncrypted: 	make([]*ckks_fv.Ciphertext, 0),
		inputEncrypted: make([]*ckks_fv.Ciphertext, 0),
		encodeCipher: 	make([]*ckks_fv.Ciphertext, 0),
	}
	return rtb, nil
}

func (rtb *RtBCipher) BootstrapCipher( ct *ckks_fv.Ciphertext ) *ckks_fv.Ciphertext{
	// for ct.Level() > 1 {
	// 	rtb.DropLevel(ct, 1)
	// }
	// ct.SetScale(math.Exp2(math.Round(math.Log2(ct.Scale())))) // rounds to the nearest power of two
	ct.SetScale(rtb.params.Scale())
	if true {
		ct = rtb.BootstrappReal(ct)
	} else {
		ct = rtb.BootReEnc(ct)		
	}
	return ct
}

func (rtb *RtBCipher) BootReEnc( ct *ckks_fv.Ciphertext ) *ckks_fv.Ciphertext {
	valuesTest := (*rtb.encoder).DecodeComplex( (*rtb.decryptor).DecryptNew(ct), rtb.params.LogSlots())
	// fmt.Println("BootReEnc debug")
	// PrintVectorTrunc(valuesTest, 7, 3)
	plainNew := ckks_fv.NewPlaintextCKKS(rtb.params, rtb.remainingLevel, rtb.params.Scale())
	(*rtb.encoder).EncodeComplex(plainNew, valuesTest, rtb.params.LogSlots())
	rtnCipher := ckks_fv.NewCiphertextCKKS(rtb.params, 1, plainNew.Level(), rtb.params.Scale())
	(*rtb.encryptor).Encrypt( plainNew, rtnCipher )
	return rtnCipher
}

func (rtb *RtBCipher) AddScalar(ct0 *ckks_fv.Ciphertext, scalar float64, Out *ckks_fv.Ciphertext) {
	rtb.AddConst(ct0, scalar, Out)
}

func (rtb *RtBCipher) SubScalar(ct0 *ckks_fv.Ciphertext, scalar float64, Out *ckks_fv.Ciphertext){
	rtb.AddConst(ct0, -scalar, Out)
}

// //boolean function constructed from Heaan using CKKS scheme
func (rtb *RtBCipher) XOR(ct0, ct1 *ckks_fv.Ciphertext, Out *ckks_fv.Ciphertext) {
  	// (x-y)^2
  	rtb.Sub(ct0, ct1, Out)
	rtb.Power(Out, 2, Out)
}

// //boolean function constructed from Heaan using CKKS scheme
func (rtb *RtBCipher) XORNew(ct0, ct1 *ckks_fv.Ciphertext) (*ckks_fv.Ciphertext) {
	// (x-y)^2
	res := rtb.SubNew(ct0, ct1)
	rtb.Power(res, 2, res)
	return res
}

func (rtb *RtBCipher) NOT( ct0 *ckks_fv.Ciphertext, Out *ckks_fv.Ciphertext ){
  	// (1-x)
	rtb.Neg(ct0, Out)
	rtb.AddConst(Out, 1.0, Out)
}

func (rtb *RtBCipher) NOTNew( ct0 *ckks_fv.Ciphertext ) (*ckks_fv.Ciphertext) {
	// (1-x)
	Out := rtb.NegNew(ct0)
	rtb.AddConst(Out, 1.0, Out)
	return Out
}

func (rtb *RtBCipher) AND( ct0, ct1 *ckks_fv.Ciphertext, Out *ckks_fv.Ciphertext){
	// xy
	rtb.MulRelin(ct0, ct1, Out)
	rtb.RescaleMany(Out, 1, Out)
}

func (rtb *RtBCipher) ANDNew( ct0, ct1 *ckks_fv.Ciphertext ) (*ckks_fv.Ciphertext) {
	// xy
	Out := rtb.MulRelinNew(ct0, ct1)
	rtb.RescaleMany(Out, 1, Out)
	return Out
}

func (rtb *RtBCipher) OR(ct0, ct1 *ckks_fv.Ciphertext, Out *ckks_fv.Ciphertext){
  	//x + y - x\cdot y
	rtb.MulRelin(ct0, ct1, Out)
	rtb.RescaleMany(Out, 1, Out)
	rtb.Sub(ct1, Out, Out)
	rtb.Add(Out, ct0, Out)
}

func (rtb *RtBCipher) ORNew(ct0, ct1 *ckks_fv.Ciphertext ) (*ckks_fv.Ciphertext) {
	//x + y - x\cdot y
  	Out := rtb.MulRelinNew(ct0, ct1)
	rtb.RescaleMany(Out, 1, Out)
	rtb.Sub(ct1, Out, Out)
	rtb.Add(Out, ct0, Out)
  return Out
}

func (rtb *RtBCipher) CleanReal(ct *ckks_fv.Ciphertext){
	squ := rtb.PowerNew(ct, 2)
	cube :=	rtb.MulRelinNew(squ, ct)
	rtb.RescaleMany(cube, 1, cube)
    // computation of 3x^2-2x^3 which is a rough-but-good-enough approximation of
    // the sign fucntion 
	rtb.Add(squ, squ, ct)
	rtb.Add(squ, ct, ct)
	rtb.Sub(ct, cube, ct)
	rtb.Sub(ct, cube, ct)
}

func (rtb *RtBCipher) BinaryTreeAdd(cts []*ckks_fv.Ciphertext) *ckks_fv.Ciphertext {
    for j := 1; j < len(cts); j *= 2 {
        for i := 0; i < len(cts); i += 2 * j {
            if i+j < len(cts) {
            	rtb.Add(cts[i], cts[i+j], cts[i])
            }
        }
    }
    return cts[0]
}

func (rtb *RtBCipher) ConstructNumberFromBits(cts []*ckks_fv.Ciphertext, error int) *ckks_fv.Ciphertext {
	fmt.Println("Construct bits...")
    Out :=	rtb.AddNew(cts[0], cts[1])
	rtb.Add(Out, cts[1], Out)

    errorBound := 0.5 * float64(error)
    n := len(cts)
    ctxPool := make([]*ckks_fv.Ciphertext, 0)
    for i := 2; i < n; i++ {
        if float64(i) < errorBound {
			tmp := math.Pow(2.0, float64(i)/2)
            b := rtb.MultByConstNew(cts[i], tmp)
			rtb.RescaleMany(b, 1, b)
        	rtb.Power(b, 2, b)
            if i%2 == 1 {
            	rtb.Add(b, b, b)
            }
            ctxPool = append(ctxPool, b)
        } else {
            tmp := math.Pow(2.0, float64(i)/4)
            b := rtb.MultByConstNew(cts[i], tmp)
			rtb.RescaleMany(b, 1, b)
        	rtb.Power(b, 4, b)
            copyNum := int( math.Pow(2.0, float64(i%4) ) )
            for j := 0; j < copyNum; j++ {
                ctxPool = append(ctxPool, b)
            }
        }
    }
    sum := rtb.BinaryTreeAdd(ctxPool)
	rtb.Add(Out, sum, Out)
    return Out
}

func (rtb *RtBCipher) DebugPrint(ct *ckks_fv.Ciphertext, descriptor string) {
	fmt.Printf("Chain Index: %d, Scale: %.2f\n", ct.Level(), math.Log2(ct.Scale()) )
	result := (*rtb.encoder).DecodeComplex( (*rtb.decryptor).DecryptNew(ct), rtb.params.LogSlots())
	fmt.Println(descriptor)
	var fltVar []float64
	for _, v := range result {
		fltVar = append(fltVar, real(v))
	}
	PrintVectorTrunc( fltVar, 7, 4 )
}

// PrintVectorTrunc prints a truncated version of the vector.
func PrintVectorTrunc(vec interface{}, printSize, prec int) {
	switch v := vec.(type) {
	case []complex128:
		printComplexVectorTrunc(v, printSize, prec)
	case []float64:
		printFloatVectorTrunc(v, printSize, prec)
	default:
		fmt.Println("Unsupported type")
	}
	fmt.Println()
}

// extractComponents extracts the real and imaginary parts of a complex number.
func extractComponents(value complex128) (float64, float64) {
	return real(value), imag(value)
}

// getDelimiter returns the appropriate delimiter for the print output.
func getDelimiter(i, maxIndex int) string {
	if i != maxIndex {
		return ", "
	}
	return " ]\n"
}

// printComplexVectorTrunc prints a truncated version of a complex vector.
func printComplexVectorTrunc(v []complex128, printSize, prec int) {
	lenV := len(v)
	if lenV <= 2*printSize {
		fmt.Printf("[")
		for i, value := range v {
			cReal, cImag := extractComponents(value)
			fmt.Printf(" (%.*f,+ %.*fi)%s", prec, cReal, prec, cImag, getDelimiter(i, lenV-1))
		}
	} else {
		fmt.Printf("[")
		for i := 0; i < printSize; i++ {
			cReal, cImag := extractComponents(v[i])
			fmt.Printf(" (%.*f+%.*fi)%s", prec, cReal, prec, cImag, getDelimiter(i, printSize-1))
		}
		fmt.Printf(" ... ")
		for i := lenV - printSize; i < lenV; i++ {
			cReal, cImag := extractComponents(v[i])
			fmt.Printf(" (%.*f+%.*fi)%s", prec, cReal, prec, cImag, getDelimiter(i, lenV-1))
		}
	}
}

// printFloatVectorTrunc prints a truncated version of a float vector.
func printFloatVectorTrunc(v []float64, printSize, prec int) {
	lenV := len(v)
	if lenV <= 2*printSize {
		fmt.Printf("[")
		for i, value := range v {
			fmt.Printf(" %.*f%s", prec, value, getDelimiter(i, lenV-1))
		}
	} else {
		fmt.Printf("[")
		for i := 0; i < printSize; i++ {
			fmt.Printf(" %.*f%s", prec, v[i], getDelimiter(i, printSize-1))
		}
		fmt.Printf(" ...")
		for i := lenV - printSize; i < lenV; i++ {
			fmt.Printf(" %.*f%s", prec, v[i], getDelimiter(i, lenV-1))
		}
	}
}