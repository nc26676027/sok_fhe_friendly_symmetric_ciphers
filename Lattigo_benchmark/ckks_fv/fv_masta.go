package ckks_fv

import (
	"fmt"
	"strconv"

	"github.com/ldsec/lattigo/v2/ring"
	"golang.org/x/crypto/sha3"
)
type MastaParam struct {
	Blocksize    int
	PlainModulus uint64
	NumRound     int
}

const (
	MASTA5S = iota
	MASTA5M
	MASTA5L
	MASTA6S
	MASTA6M
	MASTA6L
	MASTA7S
	MASTA7M
	MASTA7L
)

var MastaParams = []MastaParam{
	{	//MASTA5S	
		Blocksize:    64,
		// PlainModulus: 65537,
		PlainModulus: 0x1fc0001,
		NumRound:     5,
	},
	{	//MASTA5M
		Blocksize:    64,
		// PlainModulus: 8088322049,
		PlainModulus: 0x3ee0001,
		NumRound:     5,
	},
	{	//MASTA5L
		Blocksize:    64,
		PlainModulus: 1096486890805657601,
		NumRound:     5,
	},
	{	//MASTA6S	
		Blocksize:    32,
		// PlainModulus: 65537,
		PlainModulus: 0x1fc0001,
		NumRound:     6,
	},
	{	//MASTA6M
		Blocksize:    32,
		// PlainModulus: 8088322049,
		PlainModulus: 0x3ee0001,
		NumRound:     6,
	},
	{	//MASTA6L
		Blocksize:    32,
		PlainModulus: 1096486890805657601,
		NumRound:     6,
	},
	{	//MASTA7S	
		Blocksize:    16,
		// PlainModulus: 65537,
		PlainModulus: 0x1fc0001,
		NumRound:     7,
	},
	{	//MASTA7M
		Blocksize:    16,
		// PlainModulus: 8088322049,
		PlainModulus: 0x3ee0001,
		NumRound:     7,
	},
	{	//MASTA7L
		Blocksize:    16,
		PlainModulus: 1096486890805657601,
		NumRound:     7,
	},
}

type MFVMasta interface {
	Crypt(nonce [][]byte, counter []byte, pastaModDown []int) []*Ciphertext
	CryptNoModSwitch( nonce [][]byte, counter []byte ) []*Ciphertext
	CryptAutoModSwitch(nonce [][]byte, counter []byte, noiseEstimator MFVNoiseEstimator) (res []*Ciphertext, pastaModDown []int)
	Reset(symmetricKey []*Ciphertext, nbInitModDown int)
	EncKey(key []uint64) (res []*Ciphertext)
}

type mfvMasta struct {

	mastaParam   int
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
	mkCt []*Ciphertext
	rc   [][][]uint64    // RoundConstants[round][state][slot]
	rcMat [][][][]uint64	//Affine matrix [round][row][column][slot]
	rcPt []*Plaintext // [state] Buffer for round constants, reuse PT in each round to reduce memory cost
	rcMatPt [][]*PlaintextMul 	//[row][column] Buffer for round Affine matrix, reuse PT in each round to reduce memory cost
	xof  []sha3.ShakeHash
}

func NewMFVMasta( key []uint64, mastaParam int, params *Parameters, encoder MFVEncoder, encryptor MFVEncryptor, decryptor MFVDecryptor, evaluator MFVEvaluator, nbInitModDown int) MFVMasta {
	masta := new(mfvMasta)

	masta.mastaParam = mastaParam
	masta.blocksize = MastaParams[mastaParam].Blocksize
	masta.numRound = MastaParams[mastaParam].NumRound
	masta.slots = params.FVSlots()
	masta.nbInitModDown = nbInitModDown

	masta.params = params
	masta.encoder = encoder
	masta.encryptor = encryptor
	masta.evaluator = evaluator
	masta.decryptor = decryptor

	masta.stCt = make([]*Ciphertext, masta.blocksize)
	masta.mkCt = make([]*Ciphertext, masta.blocksize)
	masta.rcPt = make([]*Plaintext, masta.blocksize)
	masta.rcMatPt = make([][]*PlaintextMul, masta.blocksize)
	for row := 0; row < len(masta.rcMatPt); row++ {
		masta.rcMatPt[row] = make([]*PlaintextMul, masta.blocksize)
	}

	masta.xof = make([]sha3.ShakeHash, masta.slots)
	
	masta.rc = make([][][]uint64, masta.numRound+1)
	for r := 0; r <= masta.numRound; r++ {
		masta.rc[r] = make([][]uint64, masta.blocksize)
		for st := 0; st < masta.blocksize; st++ {
			masta.rc[r][st] = make([]uint64, masta.slots)
		}
	}

	masta.rcMat = make( [][][][]uint64, masta.numRound+1 )
	for r := 0; r <= masta.numRound; r++ {
		masta.rcMat[r] = make([][][]uint64, masta.blocksize) // row size
		for row := 0; row < masta.blocksize; row++ {
			masta.rcMat[r][row] = make([][]uint64, masta.blocksize) // column size
			for col := 0; col < masta.blocksize; col++{
				masta.rcMat[r][row][col] = make([]uint64, masta.slots) // column size
			}
		}
	}

	masta.mkCt = masta.EncKey(key)
	// initial state Ct using input symmetricKey encrypted
	for i := 0; i < masta.blocksize; i++ {
		masta.stCt[i] = masta.mkCt[i].CopyNew().Ciphertext()
	}
	return masta
}

func (masta *mfvMasta) Reset(symmetricKey []*Ciphertext, nbInitModDown int) {
	// Precompute Initial States
	masta.nbInitModDown = nbInitModDown

	for i := 0; i < masta.blocksize; i++ {
		masta.stCt[i] = symmetricKey[i]
		if nbInitModDown > 0 {
			masta.evaluator.ModSwitchMany(masta.stCt[i], masta.stCt[i], nbInitModDown)
		}
	}
}

// Compute Round Constants
func (masta *mfvMasta) init(nonce [][]byte, counter []byte) {
	slots := masta.slots
    if slots <= 0 {
        panic("invalid number of slots")
    }
    if len(nonce) != slots {
        panic(fmt.Sprintf("nonce length does not match the number of slots: expected %d, got %d", slots, len(nonce)))
    }
    if len(masta.xof) != slots {
        panic("pasta.xof length does not match the number of slots")
    }

	//Each slot stores the random sampled nounce.
	for i := 0; i < slots; i++ {
		masta.xof[i] = sha3.NewShake128()
		_, err := masta.xof[i].Write(nonce[i])
		if err != nil {
            panic(fmt.Sprintf("failed to write nonce[%d]: %v", i, err))
        }
		_, err = masta.xof[i].Write(counter)
		if err != nil {
            panic(fmt.Sprintf("failed to write counter: %v", err))
        }
	}

	if len(masta.rc) != masta.numRound+1 {
        panic("pasta.rc length does not match the number of rounds")
    }
	// initial round constant used in Affine layer
	for r := 0; r <= masta.numRound; r++ {
		for st := 0; st < masta.blocksize; st++ {
			for slot := 0; slot < slots; slot++ {
				masta.rc[r][st][slot] = SampleZqx(masta.xof[slot], masta.params.PlainModulus())
			}
		}
	}
	// initial Random matrices in each round 
	for r := 0; r<= masta.numRound; r++{
		for slot := 0; slot < slots; slot++ {
			tmpMat := GetRandomMatrixMasta(masta.blocksize, masta.xof[slot], masta.params.PlainModulus())
			if tmpMat == nil {
                panic("GetRandomMatrix returned nil")
            }
			for row := 0; row < len(tmpMat); row++{
				for col := 0; col < len(tmpMat[row]); col++{
					masta.rcMat[r][row][col][slot] = tmpMat[row][col]
				}
			}		

		}
	} 

}

func (masta *mfvMasta) findBudgetInfo(noiseEstimator MFVNoiseEstimator) (maxInvBudget, minErrorBits int) {
	plainModulus := ring.NewUint(masta.params.PlainModulus())
	maxInvBudget = 0
	minErrorBits = 0
	for i := 0; i < masta.blocksize; i++ {
		invBudget := noiseEstimator.InvariantNoiseBudget(masta.stCt[i])
		errorBits := masta.params.LogQLvl(masta.stCt[i].Level()) - plainModulus.BitLen() - invBudget

		if invBudget > maxInvBudget {
			maxInvBudget = invBudget
			minErrorBits = errorBits
		}
	}
	return
}

func (masta *mfvMasta) modSwitchAuto(round int, noiseEstimator MFVNoiseEstimator, mastaModDown []int) {
	lvl := masta.stCt[0].Level()

	QiLvl := masta.params.Qi()[:lvl+1]
	LogQiLvl := make([]int, lvl+1)
	for i := 0; i < lvl+1; i++ {
		tmp := ring.NewUint(QiLvl[i])
		LogQiLvl[i] = tmp.BitLen()
	}

	invBudgetOld, errorBitsOld := masta.findBudgetInfo(noiseEstimator)
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
		tmp := masta.stCt[0].CopyNew().Ciphertext()
		masta.evaluator.ModSwitchMany(masta.stCt[0], masta.stCt[0], nbModSwitch)
		invBudgetNew, _ := masta.findBudgetInfo(noiseEstimator)

		if invBudgetOld-invBudgetNew > 3 {
			nbModSwitch--
		}
		masta.stCt[0] = tmp
	}

	if nbModSwitch > 0 {
		mastaModDown[round] = nbModSwitch
		for i := 0; i < masta.blocksize; i++ {
			masta.evaluator.ModSwitchMany(masta.stCt[i], masta.stCt[i], nbModSwitch)
		}

		invBudgetNew, errorBitsNew := masta.findBudgetInfo(noiseEstimator)
		fmt.Printf("Masta Round %d [Budget | Error] : [%v | %v] -> [%v | %v]\n", round, invBudgetOld, errorBitsOld, invBudgetNew, errorBitsNew)
		fmt.Printf("Masta modDown : %v\n\n", mastaModDown)
	}
}

func (masta *mfvMasta) modSwitch(nbSwitch int) {
	if nbSwitch <= 0 {
		return
	}
	for i := 0; i < masta.blocksize; i++ {
		masta.evaluator.ModSwitchMany(masta.stCt[i], masta.stCt[i], nbSwitch)
	}
}

// Compute ciphertexts without modulus switching
func (masta *mfvMasta) CryptNoModSwitch(nonce [][]byte, counter []byte) []*Ciphertext {

	masta.init(nonce, counter)
	for r := 0; r < masta.numRound; r++ {
		masta.linLayer(r)
		masta.sbox()
	}
	masta.linLayer(masta.numRound)
	masta.AddKey()
	return masta.stCt
}

// Compute ciphertexts with automatic modulus switching
func (masta *mfvMasta) CryptAutoModSwitch(nonce [][]byte, counter []byte, noiseEstimator MFVNoiseEstimator) ([]*Ciphertext, []int) {
	mastaModDown := make([]int, masta.numRound+1)
	mastaModDown[0] = masta.nbInitModDown

	masta.init(nonce, counter)
	for r := 0; r < masta.numRound; r++ {
		masta.linLayer(r)
		masta.sbox()
		masta.modSwitchAuto(r, noiseEstimator, mastaModDown)
	}
	masta.linLayer(masta.numRound)
	masta.AddKey()

	return masta.stCt, mastaModDown
}

// Compute ciphertexts with modulus switching as given in pastaModDown
func (masta *mfvMasta) Crypt(nonce [][]byte, counter []byte, mastaModDown []int) []*Ciphertext {
	if mastaModDown[0] != masta.nbInitModDown {
		errorString := fmt.Sprintf("nbInitModDown expected %d but %d given", masta.nbInitModDown, mastaModDown[0])
		panic(errorString)
	}
	masta.init(nonce, counter)
	for r := 0; r < masta.numRound; r++ {
		masta.linLayer(r)
		myStr := "After linear state " + strconv.Itoa(r) + ": "
		masta.DebugPrint(masta.stCt[10], myStr)

		masta.sbox()
		myStr = "After sbox R: " + strconv.Itoa(r) + ": "
		masta.DebugPrint(masta.stCt[10], myStr)

		masta.modSwitch(mastaModDown[r+1])
		myStr = "After Mod sw R: " + strconv.Itoa(r) + ": "
		masta.DebugPrint(masta.stCt[10], myStr)
	}
	masta.linLayer(masta.numRound)
	masta.DebugPrint(masta.stCt[10], "After Linlayer Fin R: ")
	masta.AddKey()
	masta.DebugPrint(masta.stCt[10], "After Add Key Fin R: ")
	return masta.stCt
}

func (masta *mfvMasta) linLayer(round int) {
	ev := masta.evaluator
	// Encode the random matrix into 
	T := masta.blocksize
	
	for row := 0; row < T; row++ {
		for col := 0; col < T; col++ {
			masta.rcMatPt[row][col] = NewPlaintextMulLvl(masta.params, masta.stCt[row].Level())
			masta.encoder.EncodeUintMulSmall( masta.rcMat[round][row][col], masta.rcMatPt[row][col] )
		}
	}
	for i := 0; i < masta.blocksize; i++ {
		masta.rcPt[i] = NewPlaintextFVLvl(masta.params, masta.stCt[i].Level())
		masta.encoder.EncodeUintSmall( masta.rc[round][i], masta.rcPt[i] )
	}
	// masta.DebugPrint(masta.stCt[0], "before Linlayer state: ")
	// masta.DebugPrintPlainMul(masta.rcMatPt[0][0], "before linlay Mat[0][0]")
	// masta.DebugPrintPlain(masta.rcPt[0], "before linlay rc[0]")

	sumPool := make([]*Ciphertext, T)
	for row := 0; row < T; row++ {

		rowElPool := make([]*Ciphertext, T)
		for col := 0; col < T; col++ {
			rowElPool[col] = ev.MulNew( masta.stCt[col], masta.rcMatPt[row][col] )	
		}

		for i := 1; i < len(rowElPool); i++ {
			ev.Add(rowElPool[i], rowElPool[0], rowElPool[0])
			// if i % 4 == 0 {
			// 	ev.Reduce(rowElPool[0], rowElPool[0])
			// }
		}
		sumPool[row] = rowElPool[0]
	}
	
	for row := 0; row < T; row++ {
		// if round == masta.numRound {
		// 	myStr := "sumPool[" + strconv.Itoa(row) + "]: "
		// 	masta.DebugPrint(sumPool[row], myStr)
		// }
		ev.Add(sumPool[row], masta.rcPt[row], masta.stCt[row])
	}
}

func (masta *mfvMasta) sbox() {
	ev := masta.evaluator
	stCopy := make([]*Ciphertext, 2)
	stCopy[0] = masta.stCt[0].CopyNew().Ciphertext()
	stCopy[1] = masta.stCt[1].CopyNew().Ciphertext()

	for el := 0; el < masta.blocksize - 2; el++ {
		ind1 := (el + 1) % masta.blocksize
		ind2 := (ind1 + 1) % masta.blocksize
		tmp := ev.MulNew(masta.stCt[ind1], masta.stCt[ind2])
		ev.Relinearize(tmp, tmp)
		ev.AddNoMod(tmp, masta.stCt[ind2], tmp)
		ev.Add(masta.stCt[el], tmp, masta.stCt[el])
	}

	ind1 := masta.blocksize - 1
	ind2 := 0
	tmp := ev.MulNew(masta.stCt[ind1], stCopy[ind2])
	ev.Relinearize(tmp, tmp)
	ev.AddNoMod(tmp, stCopy[ind2], tmp)
	ev.Add(masta.stCt[masta.blocksize - 2], tmp, masta.stCt[masta.blocksize - 2])

	ind1 = 0
	ind2 = 1
	tmp = ev.MulNew(stCopy[ind1], stCopy[ind2])
	ev.Relinearize(tmp, tmp)
	ev.AddNoMod(tmp, stCopy[ind2], tmp)
	ev.Add(masta.stCt[masta.blocksize - 1], tmp, masta.stCt[masta.blocksize - 1])
}

func (masta *mfvMasta) AddKey() {
	ev := masta.evaluator
	for  el := 0; el < masta.blocksize; el++ {
		if el == 10 {
			fmt.Println("mkSc level: ", masta.mkCt[el].Level(), "stSc level: ", masta.stCt[el].Level())
		}
		
		nbSwitch := masta.mkCt[el].Level() - masta.stCt[el].Level()

		ev.ModSwitchMany(masta.mkCt[el], masta.mkCt[el], nbSwitch)

		ev.Add(masta.stCt[el], masta.mkCt[el], masta.stCt[el])
	}
}

func (masta *mfvMasta) EncKey(key []uint64) (res []*Ciphertext) {
	slots := masta.slots
	res = make([]*Ciphertext, masta.blocksize)

	for i := 0; i < masta.blocksize; i++ {
		dupKey := make([]uint64, slots)
		for j := 0; j < slots; j++ {
			dupKey[j] = key[i] % masta.params.plainModulus
		}
		keyPt := NewPlaintextFV(masta.params)
		masta.encoder.EncodeUintSmall(dupKey, keyPt)
		res[i] = masta.encryptor.EncryptNew(keyPt)
		if masta.nbInitModDown > 0 {
			masta.evaluator.ModSwitchMany(res[i], res[i], masta.nbInitModDown)
		}
	}
	masta.DebugPrint(res[0], "Gen EncKey: ")
	return
}

func (masta *mfvMasta) DebugPrint(ct *Ciphertext, str string) {
	fmt.Println(str, ", level: ", ct.Level())
	valuesTest := masta.encoder.DecodeUintNew(masta.decryptor.DecryptNew(ct))
	PrintDebugVec( valuesTest, 1, 0 )
}

func (masta *mfvMasta) DebugPrintPlain(ct *Plaintext, str string) {
	fmt.Println(str)
	valuesTest := masta.encoder.DecodeUintNew(ct)
	PrintDebugVec( valuesTest, 1, 0 )
}

func (masta *mfvMasta) DebugPrintPlainMul(ct *PlaintextMul, str string) {
	fmt.Println(str)
	valuesTest := masta.encoder.DecodeUintNew(ct)
	PrintDebugVec( valuesTest, 1, 0 )
}