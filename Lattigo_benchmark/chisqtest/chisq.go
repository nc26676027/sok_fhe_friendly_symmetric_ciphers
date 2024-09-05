package chisqtest

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
)


const EPSILON float64 = 1.0E-08


func ReadSNPFile(headers *[]string, dataColumns *[][]float64, y *[]float64, dataFileName string, N, M int) error {
	fileName := dataFileName + ".csv"
	fmt.Fprintf(os.Stderr, "file name = %s\n", fileName)

	file, err := os.Open(fileName)
	if err != nil {
		return fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.TrimLeadingSpace = true

	// Read the header line
	line, err := reader.Read()
	if err != nil {
		return fmt.Errorf("failed to read header line: %w", err)
	}

	cols := len(line)
	for i, value := range line {
		if value != "" && i > 4 && i < M+5 {
			*headers = append(*headers, value)
		}
	}

	count := 0
	for {
		if count >= N {
			break
		}
		line, err := reader.Read()
		if err != nil {
			break
		}

		if len(line) > 2 {
			yval, err := strconv.ParseFloat(line[1], 64)
			if err != nil {
				return fmt.Errorf("failed to parse float from line[1]: %w", err)
			}
			*y = append(*y, yval)

			var row []float64
			for i := 5; i < cols; i++ {
				val, err := strconv.ParseFloat(line[i], 64)
				if err != nil {
					return fmt.Errorf("failed to parse float from line[%d]: %w", i, err)
				}
				row = append(row, val)
			}
			*dataColumns = append(*dataColumns, row)
		}
		count++
	}

	fmt.Printf("Read in data: %s\n", dataFileName)
	return nil
}

// BS 
func BS(z float64) float64 {
	y := math.Exp(-z * z / 2)
	return math.Sqrt(1-y) * (31*y/200 - 341*y*y/8000) / math.Sqrt(math.Pi)
}

// normalCFD 
func normalCFD(value float64) float64 {
	return 0.5 * math.Erfc( math.Sqrt2/(-value) )
}

// sf 
func sf(value float64) float64 {
	return 1 - normalCFD(value)
}

func Equal(a, b float64) bool {
	return EPSILON > math.Abs(a-b)
}

func Less(a, b float64) bool {
	return (a - b) < (-EPSILON)
}

func Greater(a, b float64) bool {
	return (a - b) > EPSILON
}

// IncompleteGamma 
func IncompleteGamma(val, p float64) float64 {
	if !Greater(val, 0) || !Greater(p, 0) {
		return 0
	}
	LgammaP, _ := math.Lgamma(p)
	expValue := p*math.Log(val) - val - LgammaP
	if Less(expValue, math.Log(1.0E-37)) { // underflow
		return 0
	}
	factor := math.Exp(expValue)
	if !Greater(val, 1) || Less(val, p) {
		igamma := 1.0
		term := 1.0
		for i := 1; Greater(term, EPSILON); i++ {
			term *= (val / (p + float64(i)))
			igamma += term
		}
		return (igamma * factor / p)
	}

	pn := [6]float64{1, val, val + 1, val * (2 + val - p)}
	upperIncGamma := pn[2] / pn[3]
	for j := 1; ; j++ {
		a := (float64(j)+1)*2 + val - p
		b := (1 + float64(j) - p) * float64(j)
		pn[4] = a*pn[2] - b*pn[0]
		pn[5] = a*pn[3] - b*pn[1]
		if !Equal(pn[5], 0) {
			rn := pn[4] / pn[5]
			diff := math.Abs(upperIncGamma - rn)
			if !Greater(diff, EPSILON) && !Greater(diff, (EPSILON*rn)) {
				return (1 - factor*upperIncGamma)
			}
			upperIncGamma = rn
		}
		for i := 0; i < 4; i++ {
			pn[i] = pn[i+2]
		}
		if !Greater(1.0E+37, math.Abs(pn[4])) { // overflow
			for i := 0; i < 4; i++ {
				pn[i] = pn[i] / 1.0E+37
			}
		}
	}

}

func ToComplexVec( input []float64 ) []complex128 {
	res := make([]complex128, len(input))
	for i, el := range input {
		res[i] = complex(el, 0)
	}
	return res
}

func BinaryTreeAdd(vector []*ckks_fv.Ciphertext, evaluator ckks_fv.CKKSEvaluator) *ckks_fv.Ciphertext {

	for j := 1; j < len(vector); j=j*2 {
		for i := 0; i<len(vector); i = i + 2*j {
			if (i+j)<len(vector) {
				evaluator.Add(vector[i], vector[i+j], vector[i])
			}				
		}
	}
	return vector[0]
}
 

func RunChi2(SNPDir, SNPFileName, pValue, Runtime, SampleSize, SNPs string) {
	N, err := strconv.Atoi(SampleSize)
	if err != nil {
		log.Fatalf("SampleSize turn erro: %v", err)
	}
	M, err := strconv.Atoi(SNPs)
	if err != nil {
		log.Fatalf("SNPs turn erro: %v", err)
	}
	scalingFactor := 0.1*math.Pow(float64(N), -2)

	var headersS []string
	var sData [][]float64
	var yData []float64

	startTime := time.Now()

	err = ReadSNPFile(&headersS, &sData, &yData, SNPDir+"/"+SNPFileName, N, M)
	if err != nil {
		log.Fatalf("ReadSNPFile error: %v", err)
	}
	for i:=0;i<10;i++{
		fmt.Println(headersS[i])
	}

	fmt.Println("Number of Individuals =", len(sData))
	fmt.Println("Number of SNPs =", len(sData[0]))
	fmt.Println("Number of yData =", len(yData))
	btpParams := ckks_fv.DefaultBootstrapParams[0]
	params, err := btpParams.Params()
	if err != nil {
		panic(err)
	}
	
	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, h = %d, logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), btpParams.H, params.LogQP(), params.Levels(), math.Log2(params.Scale()), params.Sigma())

	// Scheme context and keys
	
	kgen := ckks_fv.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairSparse(192)
	rlk := kgen.GenRelinearizationKey(sk)
	encoder := ckks_fv.NewCKKSEncoder(params)
	encryptor := ckks_fv.NewCKKSEncryptorFromPk(params, pk)
	decryptor := ckks_fv.NewCKKSDecryptor(params, sk)
	evaluator := ckks_fv.NewCKKSEvaluator(params, ckks_fv.EvaluationKey{rlk, nil})
	galStepsVector := []int{ M }
	rotkeys := kgen.GenRotationKeysForRotations(galStepsVector, true, sk)
	evk := ckks_fv.EvaluationKey{Rlk: rlk, Rtks: rotkeys}

	// context.params.SetLogSlots(15)
	evaluator = ckks_fv.NewCKKSEvaluator(params, evk)
	workingLevel := 12

	sCiphertexts := make([]*ckks_fv.Ciphertext, N/2)

	for i := 0; i < N/2; i++ {
		individual := make([]float64, 2*M)
		copy(individual[:M], sData[i])
		copy(individual[M:], sData[N/2 + i])

		plaintext := encoder.EncodeComplexNTTNew( ToComplexVec(individual), params.LogSlots())
		S := encryptor.EncryptNew(plaintext)
		for S.Level() > workingLevel {
			evaluator.DropLevel(S, 1)
		}
		sCiphertexts[i] = S
	}

	yCiphertexts := make([]*ckks_fv.Ciphertext, N/2)
	for i := 0; i < N/2; i++ {
		individual := make([]float64, 2*M)
		for j:=0;j<M;j++{
			individual[j] = yData[i]
			individual[j+M] = yData[N/2 + i]
		}
		plaintext := encoder.EncodeComplexNTTNew( ToComplexVec(individual), params.LogSlots() )
		yCiphertexts[i] = encryptor.EncryptNew(plaintext)
		Y := encryptor.EncryptNew(plaintext)
		for Y.Level() > workingLevel {
			evaluator.DropLevel(Y, 1)
		}
		yCiphertexts[i] = Y
	}

	dVal := make([]complex128, M)
	for i, _ := range dVal {
		tmp := float64(N)
		dVal[i] = complex(2 * tmp, 0)
	}
	// d := encoder.EncodeComplexNTTNew( dVal, params.LogSlots())

	dValScaled := make([]complex128, M)
	for i, _ := range dValScaled {
		tmp := float64(N)
		dValScaled[i] = complex(2 * tmp * scalingFactor, 0)
	}

	dScaled := encoder.EncodeComplexNTTNew(dValScaled, params.LogSlots())

	start := time.Now()

	ySum := make([]*ckks_fv.Ciphertext, len(yCiphertexts) )
	for i := range yCiphertexts {
		ySum[i] = yCiphertexts[i].CopyNew().Ciphertext()
	}
	yU := BinaryTreeAdd(ySum, evaluator)

	var chiD, chiN, orD, orN *ckks_fv.Ciphertext 
	ySMult := make([]*ckks_fv.Ciphertext, N/2)

	for i := 0; i < N/2; i++ {
		ySMult[i] = evaluator.MulRelinNew(sCiphertexts[i], yCiphertexts[i])
		evaluator.RescaleMany(ySMult[i], 1, ySMult[i])
	}


	n11 := BinaryTreeAdd(ySMult, evaluator);
	c1 := BinaryTreeAdd(sCiphertexts, evaluator);
	
	//Rotate step
	step := M
	n11Rot := evaluator.RotateNew(n11, step)
	evaluator.Add(n11, n11Rot, n11)
	
	c1Rot := evaluator.RotateNew(c1, step)
	evaluator.Add(c1, c1Rot, c1)
	
	yURot := evaluator.RotateNew(yU, step)
	evaluator.Add(yU, yURot, yU)
	
	// printDebug( "n11: ", params, n11, decryptor, encoder)
	// printDebug( "c1: ", params, c1, decryptor, encoder)
	// printDebug( "yU: ", params, yU, decryptor, encoder)

	// r1 = 2 * yU
	r1 := evaluator.AddNew(yU, yU)
	r1Scaled := evaluator.MultByConstNew(r1, complex(float64(scalingFactor), 0))
	evaluator.RescaleMany(r1Scaled, 1, r1Scaled)
	
	c1Scaled := evaluator.MultByConstNew(c1, complex(float64(scalingFactor), 0))
	evaluator.RescaleMany(c1Scaled, 1, c1Scaled)
	
	// compute Chi2
	mult1 := evaluator.MulRelinNew(n11, dScaled)
	evaluator.RescaleMany(mult1, 1, mult1)
	
	mult2 := evaluator.MulRelinNew(c1, r1Scaled)
	evaluator.RescaleMany(mult2, 1, mult2)
	
	chiN1 := evaluator.SubNew(mult1, mult2)
	chiN = evaluator.PowerNew(chiN1, 2)
	
	// denominator
	negC1Scaled := evaluator.NegNew(c1Scaled)
	chiD1 := evaluator.AddNew(negC1Scaled, dScaled)
	chiD1 = evaluator.MulRelinNew(chiD1, c1)
	evaluator.RescaleMany(chiD1, 1, chiD1)
	
	negR1Scaled := evaluator.NegNew(r1Scaled)
	chiD2 := evaluator.AddNew(negR1Scaled, dScaled)
	chiD2 = evaluator.MulRelinNew(chiD2, r1)
	evaluator.RescaleMany(chiD2, 1, chiD2)
	
	chiD = evaluator.MulRelinNew( chiD1, chiD2 )
	evaluator.RescaleMany(chiD, 1, chiD)
	
	// Odds Ratio
	n11Scaled := evaluator.MultByConstNew(n11, complex(float64(scalingFactor), 0))
	evaluator.RescaleMany(n11Scaled, 1, n11Scaled)
	
	// denominator
	or2 := evaluator.SubNew(c1, n11)
	or3 := evaluator.SubNew(r1Scaled, n11Scaled)
	orD = evaluator.MulRelinNew(or2, or3)
	evaluator.RescaleMany(orD, 1, orD)
	
	// numerator
	or1 := evaluator.SubNew(n11Scaled, r1Scaled)
	or1 = evaluator.SubNew(or1, c1Scaled)
	or1 = evaluator.AddNew(or1, dScaled)
	evaluator.RescaleMany(or1, 1, or1)
	
	orN = evaluator.MulRelinNew(n11, or1)
	evaluator.RescaleMany(orN, 1, orN)

	endToEndTime := float64(time.Since(start).Seconds())


	pN := decryptor.DecryptNew(chiN)
	pD := decryptor.DecryptNew(chiD)
	oddN := decryptor.DecryptNew(orN)
	oddD := decryptor.DecryptNew(orD)
	msg_pN := encoder.DecodeComplex(pN, params.LogSlots())
	msg_pD := encoder.DecodeComplex(pD, params.LogSlots())
	msg_oddN := encoder.DecodeComplex(oddN, params.LogSlots())
	msg_oddD := encoder.DecodeComplex(oddD, params.LogSlots())

	chival := make([]float64, len(headersS))
	pval := make([]float64, len(headersS))
	odds := make([]float64, len(headersS))

	for i := 0; i < len(headersS); i++ {
		chival[i] = real(msg_pN[i]) * 2 * float64(N) / real(msg_pD[i])
		if chival[i] < 0 {
			chival[i] = 0
		}
		pval[i] = 1 - IncompleteGamma(chival[i]/2, 0.5)
		if pval[i] < 0 {
			pval[i] = 1e-15
		} else {
			if pval[i] == 0 {
				pval[i] = BS(math.Sqrt(chival[i]))
			}
		}
		odds[i] = real(msg_oddN[i]) / real(msg_oddD[i])
	}	

	writeResults("pValue.txt", headersS, pval)
	writeResults("odds.txt", headersS, odds)
	writeResults("chi2.txt", headersS, chival)

	// keyGenTime := float64(time.Since(startKeyGen).Seconds())
	// encryptionTime := float64(time.Since(startEnc).Seconds())
	// computationTime := float64(time.Since(startComp).Seconds())
	// decryptionTime := float64(time.Since(startDec).Seconds())
	
	// fmt.Printf("\nKey Generation Time: \t\t%.3f s\n", keyGenTime)
	// fmt.Printf("Encoding and Encryption Time: \t%.3f s\n", encryptionTime)
	// fmt.Printf("Computation Time: \t\t%.3f s\n", computationTime)
	// fmt.Printf("Decryption & Decoding Time: \t%.3f s\n", decryptionTime)
	fmt.Printf("\nEnd-to-end Runtime: \t\t%.3f s\n", endToEndTime)

	// writeRuntime("runtime.txt", keyGenTime, encryptionTime, computationTime, decryptionTime, endToEndTime)

	duration := time.Since(startTime)
	fmt.Printf("End-to-end Runtime: %.2f s\n", duration.Seconds())
}

func writeResults(filename string, headersS []string, values []float64) {
	file, err := os.Create(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	for i := 0; i < len(headersS); i++ {
		fmt.Fprintf(file, "%s\t%f\n", headersS[i], values[i])
	}
}

func writeRuntime(filename string, keyGenTime, encryptionTime, computationTime, decryptionTime, endToEndTime float64) {
    file, err := os.Create(filename)
    if err != nil {
        panic(err)
    }
    defer file.Close()

    fmt.Fprintf(file, "Key Generation Time: \t\t%.3f s\n", keyGenTime)
    fmt.Fprintf(file, "Encoding and Encryption Time: \t%.3f s\n", encryptionTime)
    fmt.Fprintf(file, "Computation Time: \t\t%.3f s\n", computationTime)
    fmt.Fprintf(file, "Decryption & Decoding Time: \t%.3f s\n", decryptionTime)
    fmt.Fprintf(file, "End-to-end Runtime: \t\t%.3f s\n", endToEndTime)
}

func printDebug( str string, params *ckks_fv.Parameters, ciphertext *ckks_fv.Ciphertext, decryptor ckks_fv.CKKSDecryptor, encoder ckks_fv.CKKSEncoder) {

	fmt.Println(str)
	valuesTest := encoder.DecodeComplex(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
}


