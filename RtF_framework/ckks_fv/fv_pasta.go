package ckks_fv

import (
	"fmt"
	"strconv"

	"github.com/ldsec/lattigo/v2/ring"
	"golang.org/x/crypto/sha3"
)
type PastaParam struct {
	Blocksize    int
	PlainModulus uint64
	NumRound     int
}

const (
	PASTA3S = iota
	PASTA3M
	PASTA3L
	PASTA4S
	PASTA4M
	PASTA4L
	PASTA5S
	PASTA5M
	PASTA5L
)

var PastaParams = []PastaParam{
	{	//PASTA3S
		Blocksize:    256,
		// PlainModulus: 65537,//17bit
		PlainModulus: 0x3ee0001,//26bit
		NumRound:     3,
	},
	{	//PASTA3M
		Blocksize:    256,
		// PlainModulus: 8088322049, // 33bit
		PlainModulus: 0x1fc0001, // 25bit
		NumRound:     3,
	},
	{ 	//PASTA3L
		Blocksize:    256,
		// PlainModulus: 1096486890805657601, //60bit
		PlainModulus: 0x3ee0001, //60bit
		NumRound:     3,
	},
	{	//PASTA4S	
		Blocksize:    64,
		// PlainModulus: 65537,
		PlainModulus: 0x1fc0001,
		NumRound:     4,
	},
	{	//PASTA4M
		Blocksize:    64,
		// PlainModulus: 8088322049,
		PlainModulus: 0x3ee0001,
		NumRound:     4,
	},
	{	//PASTA4L
		Blocksize:    64,
		PlainModulus: 1096486890805657601,
		NumRound:     4,
	},
	{	//PASTA4S	
		Blocksize:    32,
		// PlainModulus: 65537,
		PlainModulus: 0x1fc0001,
		NumRound:     5,
	},
	{	//PASTA4M
		Blocksize:    32,
		// PlainModulus: 8088322049,
		PlainModulus: 0x3ee0001,
		NumRound:     5,
	},
	{	//PASTA4L
		Blocksize:    32,
		PlainModulus: 1096486890805657601,
		NumRound:     5,
	},
}

type MFVPasta interface {
	Crypt(nonce [][]byte, counter []byte, pastaModDown []int) []*Ciphertext
	CryptNoModSwitch( nonce [][]byte, counter []byte ) []*Ciphertext
	CryptAutoModSwitch(nonce [][]byte, counter []byte, noiseEstimator MFVNoiseEstimator) (res []*Ciphertext, pastaModDown []int)
	Reset(symmetricKey []*Ciphertext, nbInitModDown int)
	EncKey(key []uint64) (res []*Ciphertext)
	DebugPrint(*Ciphertext, string )
	DebugPrintPlain(*Plaintext, string )
	DebugPrintPlainMul(*PlaintextMul, string )
}

type mfvPasta struct {

	pastaParam   int
	blocksize     int
	numRound      int
	slots         int
	nbInitModDown int

	params    *Parameters
	encoder   MFVEncoder
	encryptor MFVEncryptor
	evaluator MFVEvaluator
	decryptor MFVDecryptor

	stCt []*Ciphertext
	rc   [][][]uint64    // RoundConstants[round][state][slot]
	rcMat1 [][][][]uint64	//Affine matrix [round][row][column][slot]
	rcMat2 [][][][]uint64	//Affine matrix [round][row][column][slot]
	rcPt []*Plaintext // [state] Buffer for round constants, reuse PT in each round to reduce memory cost
	rcMat1Pt [][]*PlaintextMul 	//[row][column] Buffer for round Affine matrix, reuse PT in each round to reduce memory cost
	rcMat2Pt [][]*PlaintextMul 	//[row][column] Buffer for round Affine matrix
	xof  []sha3.ShakeHash
}

func NewMFVPasta( key []uint64, pastaParam int, params *Parameters, encoder MFVEncoder, encryptor MFVEncryptor, decryptor MFVDecryptor, evaluator MFVEvaluator, nbInitModDown int) MFVPasta {
	pasta := new(mfvPasta)

	pasta.pastaParam = pastaParam
	pasta.blocksize = PastaParams[pastaParam].Blocksize
	pasta.numRound = PastaParams[pastaParam].NumRound
	pasta.slots = params.FVSlots()
	pasta.nbInitModDown = nbInitModDown

	pasta.params = params
	pasta.encoder = encoder
	pasta.encryptor = encryptor
	pasta.evaluator = evaluator
	pasta.decryptor = decryptor

	pasta.stCt = make([]*Ciphertext, pasta.blocksize)
	pasta.rcPt = make([]*Plaintext, pasta.blocksize)
	pasta.rcMat1Pt = make([][]*PlaintextMul, pasta.blocksize/2)
	for row := 0; row < len(pasta.rcMat1Pt); row++ {
		pasta.rcMat1Pt[row] = make([]*PlaintextMul, pasta.blocksize/2)
	}
	pasta.rcMat2Pt = make([][]*PlaintextMul, pasta.blocksize/2)
	for row := 0; row < len(pasta.rcMat2Pt); row++ {
		pasta.rcMat2Pt[row] = make([]*PlaintextMul, pasta.blocksize/2)
	}
	pasta.xof = make([]sha3.ShakeHash, pasta.slots)
	
	pasta.rc = make([][][]uint64, pasta.numRound+1)
	for r := 0; r <= pasta.numRound; r++ {
		pasta.rc[r] = make([][]uint64, pasta.blocksize)
		for st := 0; st < pasta.blocksize; st++ {
			pasta.rc[r][st] = make([]uint64, pasta.slots)
		}
	}

	pasta.rcMat1 = make( [][][][]uint64, pasta.numRound+1 )
	pasta.rcMat2 = make( [][][][]uint64, pasta.numRound+1 )
	for r := 0; r <= pasta.numRound; r++ {
		pasta.rcMat1[r] = make([][][]uint64, pasta.blocksize/2) // row size
		pasta.rcMat2[r] = make([][][]uint64, pasta.blocksize/2)
		for row := 0; row < pasta.blocksize/2; row++ {
			pasta.rcMat1[r][row] = make([][]uint64, pasta.blocksize/2) // column size
			pasta.rcMat2[r][row] = make([][]uint64, pasta.blocksize/2) 
			for col := 0; col < pasta.blocksize/2; col++{
				pasta.rcMat1[r][row][col] = make([]uint64, pasta.slots) 
				pasta.rcMat2[r][row][col] = make([]uint64, pasta.slots) 		
			}
		}
	}
	symmetricKey := pasta.EncKey(key)
	fmt.Println("level: ", symmetricKey[0].Level())
	// initial state Ct using input symmetricKey encrypted
	for i := 0; i < pasta.blocksize; i++ {
		pasta.stCt[i] = symmetricKey[i]
	}
	fmt.Println("level after modsw: ", pasta.stCt[0].Level())
	return pasta
}

func (pasta *mfvPasta) Reset(symmetricKey []*Ciphertext, nbInitModDown int) {
	// Precompute Initial States
	pasta.nbInitModDown = nbInitModDown

	for i := 0; i < pasta.blocksize; i++ {
		pasta.stCt[i] = symmetricKey[i]
		if nbInitModDown > 0 {
			pasta.evaluator.ModSwitchMany(pasta.stCt[i], pasta.stCt[i], nbInitModDown)
		}
	}
}

// Compute Round Constants
func (pasta *mfvPasta) init(nonce [][]byte, counter []byte) {
	slots := pasta.slots
    if slots <= 0 {
        panic("invalid number of slots")
    }
    if len(nonce) != slots {
        panic(fmt.Sprintf("nonce length does not match the number of slots: expected %d, got %d", slots, len(nonce)))
    }
    if len(pasta.xof) != slots {
        panic("pasta.xof length does not match the number of slots")
    }

	//Each slot stores the random sampled nounce.
	for i := 0; i < slots; i++ {
		pasta.xof[i] = sha3.NewShake128()
		_, err := pasta.xof[i].Write(nonce[i])
		if err != nil {
            panic(fmt.Sprintf("failed to write nonce[%d]: %v", i, err))
        }
		_, err = pasta.xof[i].Write(counter)
		if err != nil {
            panic(fmt.Sprintf("failed to write counter: %v", err))
        }
	}

	if len(pasta.rc) != pasta.numRound+1 {
        panic("pasta.rc length does not match the number of rounds")
    }
	// initial round constant used in Affine layer
	for r := 0; r <= pasta.numRound; r++ {
		for st := 0; st < pasta.blocksize; st++ {
			for slot := 0; slot < slots; slot++ {
				pasta.rc[r][st][slot] = SampleZqx(pasta.xof[slot], pasta.params.PlainModulus())
			}
		}
	}
	// initial Random matrices in each round 
	for r := 0; r<= pasta.numRound; r++{
		for slot := 0; slot < slots; slot++ {
			tmpMat := GetRandomMatrixPasta(pasta.blocksize/2, pasta.xof[slot], pasta.params.PlainModulus())
			if tmpMat == nil {
                panic("GetRandomMatrix returned nil")
            }
			for row := 0; row < len(tmpMat); row++{
				for col := 0; col < len(tmpMat[row]); col++{
					pasta.rcMat1[r][row][col][slot] = tmpMat[row][col]
				}
			}		

			tmpMat2 := GetRandomMatrixPasta(pasta.blocksize/2, pasta.xof[slot], pasta.params.PlainModulus())
			if tmpMat2 == nil {
                panic("GetRandomMatrix returned nil")
            }
			for row := 0; row < len(tmpMat2); row++{
				for col := 0; col < len(tmpMat2[row]); col++{
					pasta.rcMat2[r][row][col][slot] = tmpMat2[row][col]
				}
			}	
		}
	} 

}

func (pasta *mfvPasta) findBudgetInfo(noiseEstimator MFVNoiseEstimator) (maxInvBudget, minErrorBits int) {
	plainModulus := ring.NewUint(pasta.params.PlainModulus())
	maxInvBudget = 0
	minErrorBits = 0
	for i := 0; i < pasta.blocksize; i++ {
		invBudget := noiseEstimator.InvariantNoiseBudget(pasta.stCt[i])
		errorBits := pasta.params.LogQLvl(pasta.stCt[i].Level()) - plainModulus.BitLen() - invBudget

		if invBudget > maxInvBudget {
			maxInvBudget = invBudget
			minErrorBits = errorBits
		}
	}
	return
}

func (pasta *mfvPasta) modSwitchAuto(round int, noiseEstimator MFVNoiseEstimator, pastaModDown []int) {
	lvl := pasta.stCt[0].Level()

	QiLvl := pasta.params.Qi()[:lvl+1]
	LogQiLvl := make([]int, lvl+1)
	for i := 0; i < lvl+1; i++ {
		tmp := ring.NewUint(QiLvl[i])
		LogQiLvl[i] = tmp.BitLen()
	}

	invBudgetOld, errorBitsOld := pasta.findBudgetInfo(noiseEstimator)
	nbModSwitch, targetErrorBits := 0, errorBitsOld
	for {
		targetErrorBits -= LogQiLvl[lvl-nbModSwitch]
		if targetErrorBits > 0 {
			nbModSwitch++
		} else {
			break
		}
	}
	if nbModSwitch != 0 {
		tmp := pasta.stCt[0].CopyNew().Ciphertext()
		pasta.evaluator.ModSwitchMany(pasta.stCt[0], pasta.stCt[0], nbModSwitch)
		invBudgetNew, _ := pasta.findBudgetInfo(noiseEstimator)

		if invBudgetOld-invBudgetNew > 3 {
			nbModSwitch--
		}
		pasta.stCt[0] = tmp
	}

	if nbModSwitch > 0 {
		pastaModDown[round] = nbModSwitch
		for i := 0; i < pasta.blocksize; i++ {
			pasta.evaluator.ModSwitchMany(pasta.stCt[i], pasta.stCt[i], nbModSwitch)
		}

		invBudgetNew, errorBitsNew := pasta.findBudgetInfo(noiseEstimator)
		fmt.Printf("Hera Round %d [Budget | Error] : [%v | %v] -> [%v | %v]\n", round, invBudgetOld, errorBitsOld, invBudgetNew, errorBitsNew)
		fmt.Printf("Hera modDown : %v\n\n", pastaModDown)
	}
}

func (pasta *mfvPasta) modSwitch(nbSwitch int) {
	if nbSwitch <= 0 {
		return
	}
	for i := 0; i < pasta.blocksize; i++ {
		pasta.evaluator.ModSwitchMany(pasta.stCt[i], pasta.stCt[i], nbSwitch)
	}
}

// Compute ciphertexts without modulus switching
func (pasta *mfvPasta) CryptNoModSwitch(nonce [][]byte, counter []byte) []*Ciphertext {

	pasta.init(nonce, counter)
	pasta.linLayer(0)
	pasta.DebugPrint(pasta.stCt[10], "LinLayer 10")
	pasta.DebugPrint(pasta.stCt[pasta.blocksize/2+10], "LinLayer t+10")
	// panic("Linlayer Done")
	for r := 1; r < pasta.numRound; r++ {
		pasta.sboxFeistel() 
		str_debug := "sboxFeistel[10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[10], str_debug)
		str_debug = "sboxFeistel[t+10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[pasta.blocksize/2+10], str_debug)
		fmt.Println("level before modsw", pasta.stCt[0].Level())
		for i := 0; i < pasta.blocksize; i++ {
			pasta.evaluator.ModSwitchMany(pasta.stCt[i], pasta.stCt[i], 2)
		}
		fmt.Println("level after modsw", pasta.stCt[0].Level())
		pasta.linLayer(r)
		str_debug = "linLayer[10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[10], str_debug)
		str_debug = "linLayer[r+10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[pasta.blocksize/2+10], str_debug)
	}
	pasta.sboxCube()
	pasta.DebugPrint(pasta.stCt[0], "SboxCube Fin")
	pasta.linLayer(pasta.numRound)
	pasta.DebugPrint(pasta.stCt[0], "linLayer Fin")
	return pasta.stCt[:pasta.blocksize/2]
}

// Compute ciphertexts with automatic modulus switching
func (pasta *mfvPasta) CryptAutoModSwitch(nonce [][]byte, counter []byte, noiseEstimator MFVNoiseEstimator) ([]*Ciphertext, []int) {
	pastaModDown := make([]int, pasta.numRound+1)
	pastaModDown[0] = pasta.nbInitModDown

	pasta.init(nonce, counter)
	pasta.linLayer(0)
	pasta.DebugPrint(pasta.stCt[10], "LinLayer 10")
	pasta.DebugPrint(pasta.stCt[pasta.blocksize/2+10], "LinLayer t+10")
	for r := 1; r < pasta.numRound; r++ {
		pasta.sboxFeistel()
		str_debug := "sboxFeistel[10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[10], str_debug)
		str_debug = "sboxFeistel[t+10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[pasta.blocksize/2+10], str_debug)

		pasta.modSwitchAuto(r, noiseEstimator, pastaModDown)
		pasta.linLayer(r)
		str_debug = "linLayer[10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[10], str_debug)
		str_debug = "linLayer[r+10], Round: " + strconv.Itoa(r)
		pasta.DebugPrint(pasta.stCt[pasta.blocksize/2+10], str_debug)
	}
	pasta.sboxCube()
	pasta.modSwitchAuto(pasta.numRound, noiseEstimator, pastaModDown)
	pasta.linLayer(pasta.numRound)

	return pasta.stCt[:pasta.blocksize/2], pastaModDown
}

// Compute ciphertexts with modulus switching as given in pastaModDown
func (pasta *mfvPasta) Crypt(nonce [][]byte, counter []byte, pastaModDown []int) []*Ciphertext {
	if pastaModDown[0] != pasta.nbInitModDown {
		errorString := fmt.Sprintf("nbInitModDown expected %d but %d given", pasta.nbInitModDown, pastaModDown[0])
		panic(errorString)
	}
	pasta.init(nonce, counter)
	pasta.linLayer(0)
	for r := 1; r < pasta.numRound; r++ {
		pasta.sboxFeistel()
		pasta.modSwitch(pastaModDown[r])
		pasta.linLayer(r)
	}
	pasta.sboxCube()
	pasta.modSwitch(pastaModDown[pasta.numRound])
	pasta.linLayer(pasta.numRound)
	
	return pasta.stCt[:pasta.blocksize/2]
}

func (pasta *mfvPasta) linLayer(round int) {
	ev := pasta.evaluator

	// Encode the random matrix into 
	PastaT := pasta.blocksize/2
	
	for row := 0; row < PastaT; row++ {
		for col := 0; col < PastaT; col++ {
			pasta.rcMat1Pt[row][col] = NewPlaintextMulLvl(pasta.params, pasta.stCt[row].Level())
			pasta.encoder.EncodeUintMulSmall( pasta.rcMat1[round][row][col], pasta.rcMat1Pt[row][col] )

			pasta.rcMat2Pt[row][col] = NewPlaintextMulLvl(pasta.params, pasta.stCt[row+PastaT].Level())
			pasta.encoder.EncodeUintMulSmall( pasta.rcMat2[round][row][col], pasta.rcMat2Pt[row][col] )
		}
	}

	for i := 0; i < pasta.blocksize; i++ {
		pasta.rcPt[i] = NewPlaintextFVLvl(pasta.params, pasta.stCt[i].Level()) // may have error, can change to RingT Plaintext
		pasta.encoder.EncodeUintSmall( pasta.rc[round][i], pasta.rcPt[i] )
	}

	sumPool := make([]*Ciphertext, PastaT)
	sumPool2 := make([]*Ciphertext, PastaT)

	for row := 0; row < PastaT; row++ {
		rowElPool := make([]*Ciphertext, PastaT)
		rowElPool2 := make([]*Ciphertext, PastaT)
		
		for col := 0; col < PastaT; col++ {
			rowElPool[col] = ev.MulNew( pasta.stCt[col], pasta.rcMat1Pt[row][col] )	
			rowElPool2[col] = ev.MulNew( pasta.stCt[PastaT+col], pasta.rcMat2Pt[row][col] )	
		}

		for i := 1; i < PastaT; i++ {
			ev.AddNoMod(rowElPool[i], rowElPool[0], rowElPool[0])
			ev.AddNoMod(rowElPool2[i], rowElPool2[0], rowElPool2[0])
			if i % 6 == 0 {
				ev.Reduce(rowElPool[0], rowElPool[0])
				ev.Reduce(rowElPool2[0], rowElPool2[0])
			}
		}
		sumPool[row] = rowElPool[0]
		sumPool2[row] = rowElPool2[0]
	}
	for row := 0; row < PastaT; row++ {
		pasta.stCt[row] = ev.AddNoModNew(sumPool[row], pasta.rcPt[row])
		pasta.stCt[PastaT+row] = ev.AddNoModNew(sumPool2[row], pasta.rcPt[PastaT+row])
		ev.Reduce(pasta.stCt[row], pasta.stCt[row])
		ev.Reduce(pasta.stCt[PastaT+row], pasta.stCt[PastaT+row])
	}
	// mix tow matrix
	for i := 0; i < PastaT; i++ {
		sum := ev.AddNew(pasta.stCt[i], pasta.stCt[PastaT+i])
		pasta.stCt[i] = ev.AddNew(pasta.stCt[i], sum)
		pasta.stCt[PastaT+i] = ev.AddNew(pasta.stCt[PastaT+i], sum)
	}
	// pasta.DebugPrint(pasta.stCt[0], "after mix matrix")
}

func (pasta *mfvPasta) sboxCube() {
	ev := pasta.evaluator
	for st := 0; st < pasta.blocksize; st++ {
		x2 := ev.MulNew(pasta.stCt[st], pasta.stCt[st])
		y2 := ev.RelinearizeNew(x2)
		x3 := ev.MulNew(y2, pasta.stCt[st])
		pasta.stCt[st] = ev.RelinearizeNew(x3)
	}
}

func (pasta *mfvPasta) sboxFeistel() {
	ev := pasta.evaluator
	for i := pasta.blocksize - 1; i > 0; i-- {
		tmp := ev.MulNew(pasta.stCt[i-1], pasta.stCt[i-1])
		ev.Relinearize(tmp, tmp)
		ev.Add(pasta.stCt[i], tmp, pasta.stCt[i])
	}
}

func (pasta *mfvPasta) EncKey(key []uint64) (res []*Ciphertext) {
	slots := pasta.slots
	res = make([]*Ciphertext, pasta.blocksize)

	for i := 0; i < pasta.blocksize; i++ {
		dupKey := make([]uint64, slots)
		for j := 0; j < slots; j++ {
			dupKey[j] = key[i] % pasta.params.plainModulus
		}
		keyPt := NewPlaintextFV(pasta.params)
		pasta.encoder.EncodeUintSmall(dupKey, keyPt)
		res[i] = pasta.encryptor.EncryptNew(keyPt)
		if pasta.nbInitModDown > 0 {
			pasta.evaluator.ModSwitchMany(res[i], res[i], pasta.nbInitModDown)
		}
	}
	return res
}

func (pasta *mfvPasta) DebugPrint(ct *Ciphertext, str string) {
	fmt.Println(str, ", level: ", ct.Level())
	valuesTest := pasta.encoder.DecodeUintNew(pasta.decryptor.DecryptNew(ct))
	PrintDebugVec( valuesTest, 1, 0 )
}

func (pasta *mfvPasta) DebugPrintPlain(ct *Plaintext, str string) {
	fmt.Println(str)
	valuesTest := pasta.encoder.DecodeUintNew(ct)
	PrintDebugVec( valuesTest, 1, 0 )
}

func (pasta *mfvPasta) DebugPrintPlainMul(ct *PlaintextMul, str string) {
	fmt.Println(str)
	valuesTest := pasta.encoder.DecodeUintNew(ct)
	PrintDebugVec( valuesTest, 1, 0 )
}

func PrintDebugVec( vec []uint64, start, end int ) error {
	if start > len(vec) || end > len(vec) {
		panic("Print element size is larger than vector size!")
	}
	starIDX := len( vec ) - end
	fmt.Println( vec[:start],"  ...  ", vec[starIDX:])
	return nil
}

func PrintDebugVecC( vec []complex128, start, end int ) error {
	if start > len(vec) || end > len(vec) {
		panic("Print element size is larger than vector size!")
	}
	starIDX := len( vec ) - end
	fmt.Println( vec[:start],"  ...  ", vec[starIDX:])
	return nil
}