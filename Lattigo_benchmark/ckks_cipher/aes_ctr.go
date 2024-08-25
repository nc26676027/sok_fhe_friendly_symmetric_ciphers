package ckks_cipher

import (
	"fmt"
	"math"
	"strconv"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
)

type AESCtr struct {
    *RtBCipher
	iv         	[]byte
    blockSize	int
	keySize	    int
	rounds		int
	// bitIndex   	[]*BitSet
    bitSbox	   	[]*BitSet
    sboxMonomialOrder	[]*BitSet
	allZeroIn	bool
}

func NewAESCtr(key_ []uint8, params_ *ckks_fv.Parameters, btpParams_ *ckks_fv.BootstrappingParameters, btpKey_ ckks_fv.BootstrappingKey, encoder_ *ckks_fv.CKKSEncoder, encryptor_ *ckks_fv.CKKSEncryptor, decryptor_ *ckks_fv.CKKSDecryptor, iv_ []byte) (*AESCtr, error) {
    rtb, err := NewRtBCipher(key_, params_, btpParams_, btpKey_, encoder_, encryptor_, decryptor_)
    if err!= nil {
        return nil, err
    }
    aes := &AESCtr{
        RtBCipher: 	rtb,
        blockSize:	128,
		keySize: 	128,
		rounds:		10,
        iv:         iv_,
		allZeroIn: 	true, //debug mode
    }
	var monomialOrder []*BitSet
    for i := 0; i < 8; i++ {
        x := NewBitSet(8)
		x.Set(1<<i)
        monomialOrder = append(monomialOrder, x)
    }

    aes.sboxMonomialOrder, _ = LayeredCombineBin(monomialOrder)

    for i := 0; i < 256; i++ {
		tmp := NewBitSet(8)
		tmp.Set( int(AESSbox[i]) )
        aes.bitSbox = append(aes.bitSbox, tmp)
    }

	return aes, err
}

func (aes *AESCtr) DebugTest(ciphertexts []byte, bits int) ([]*ckks_fv.Ciphertext, error) {
	if aes.allZeroIn {
		bits = aes.params.Slots() * aes.blockSize
	}
	numBlocks := int(math.Ceil(float64(bits) / float64(aes.blockSize)))

	iv := NewBitSet(aes.blockSize)
	iv.Set(0) // set iv all zero
	aes.EncryptKey()
	aes.EncryptInput(iv, numBlocks)
	aes.DebugPrint(aes.inputEncrypted[0], "after Enc in")
	state := aes.inputEncrypted
	for i := 0; i < len(state); i++ {
		if state[i] == nil { 
			strMsg := "state[" + strconv.Itoa(i) + "] is nil, exit!"	
			panic(strMsg) 
		}
		for j := 0; j < aes.totalLevel-aes.remainingLevel+10; j++ {
			(*aes).DropLevel(state[i], 1)
		}
	}

	for i:=0;i<8;i++{
		aes.DebugPrint(state[i], "before sbox: \n")
	}	

	aes.aesSubbyteLUT(state[:8])

	// valuesTest := (*aes.encoder).DecodeComplex( (*aes.decryptor).DecryptNew(state[0]), aes.params.LogSlots())
	// fmt.Println("BootReEnc debug")
	// PrintVectorTrunc(valuesTest, 7, 3)

	for i:=0;i<8;i++{
		aes.DebugPrint(state[i], "after sbox: \n")
	}
	return state, nil
}

func (aes *AESCtr) HEDecrypt(ciphertexts []uint8, bits int) []*ckks_fv.Ciphertext {
	if aes.allZeroIn {
		bits = aes.params.Slots() * aes.blockSize
	}
	numBlock := int(math.Ceil(float64(bits) / float64(aes.blockSize)))
	iv := NewBitSet(aes.blockSize)
	iv.Set(0) // set iv all zero
	aes.EncryptKey()
	aes.EncryptInput(iv, numBlock)

	for i := 0; i < len(aes.inputEncrypted); i++ {
		for j := 0; j < aes.totalLevel-aes.remainingLevel+5; j++ {
			(*aes).DropLevel(aes.inputEncrypted[i], 1)
		}
	}
	startAES := time.Now()

	// AES encryption **********************************************
	state := aes.AddWhiteKey(aes.inputEncrypted, aes.keyEncrypted)
	for i := 1; i < 10; i++ {
		fmt.Printf("round iterator : %d\n", i)
		aes.aesRoundFunction(state, aes.keyEncrypted)
	}
	fmt.Println("round iterator : last round")
	// aes.AesRoundFunction(state, aes.keyEncrypted)
	// AES encryption **********************************************

	endAES := time.Now()
	durationAES := endAES.Sub(startAES)
	fmt.Printf("代码执行时间：%d 秒 :: %d 毫秒\n", int(durationAES.Seconds()), int(durationAES.Milliseconds())%1000)

	// Add cipher
	// encode_ciphertext(ciphertexts, num_block);
	for i := 0; i < len(state); i++ {
		if (i%8 == 1) || (i%8 == 2) || (i%8 == 4) || (i%8 == 5) {
			aes.NOT(state[i], state[i])
		}
	}

	ch := make(chan int, len(state) )
	for i := 0; i < len(state); i++ {
		go func(i int) {
			aes.XOR(state[i], aes.encodeCipher[i], state[i])
			ch <- i
		}(i)
	}
	for i := 0; i < len(state); i++ {
		<-ch
	}

	return state
}

func (aes *AESCtr) EncryptKey() {
    if aes.encoder == nil || aes.encryptor == nil {
        panic("encoder or encryptor is not initialized")
    }

    // fmt.Println("Starting encryption")
    aes.keyEncrypted = make([]*ckks_fv.Ciphertext, 0, aes.keySize)

    for i := 0; i < aes.keySize; i++ {
        if i >= len(aes.symmetricKey)*8 {
            panic("input symmetric key size is not match!")
        }

        bit := (aes.symmetricKey[i/8] >> uint(i%8)) & 1

        skDuplicated := make([]complex128, aes.params.Slots())
        for j := 0; j < aes.params.Slots(); j++ {
            skDuplicated[j] = complex(float64(bit), 0.0)
        }

        // fmt.Printf("Encoding bit %d which is %d\n", i, bit)
        skPlain := (*aes.encoder).EncodeComplexAtLvlNew(aes.remainingLevel, skDuplicated, aes.params.LogSlots())
        if skPlain == nil {
            panic("skPlain is nil after encoding")
        }

        skBitEncrypted := new(ckks_fv.Ciphertext)
        // fmt.Println("Before encryption")
        skBitEncrypted = (*aes.encryptor).EncryptNew(skPlain)
        // fmt.Println("After encryption")

        if skBitEncrypted == nil {
            panic("skBitEncrypted is nil after encryption")
        }

        aes.keyEncrypted = append(aes.keyEncrypted, skBitEncrypted)
    }

    // fmt.Println("Encryption completed")
}

func (aes *AESCtr) EncryptInput(iv *BitSet, numBlock int) {
	aes.inputEncrypted = make([]*ckks_fv.Ciphertext, 0, aes.blockSize)
	inputData := make([]*BitSet, aes.params.Slots() )
	for i, block := range inputData{
		block = NewBitSet(aes.blockSize)
		inputData[i] = block
	}

	for i := range inputData {
		if aes.allZeroIn {
			inputData[i].Set(0)
		} else {
			inputData[i].Set( int(ctr(iv, uint64(i+1) ).ToULong()) )
		}
	}
	for i := 0; i < aes.blockSize; i++ {
		stateBatched := make([]complex128, aes.params.Slots())
		for j := 0; j < aes.params.Slots(); j++ {
			stateBatched[j] = complex( float64(inputData[j].bits[i]), 0 )
		}
		stateBatchedPlain := (*aes.encoder).EncodeComplexNTTAtLvlNew( aes.remainingLevel, stateBatched, aes.params.LogSlots() )
		stateBatchedEncrypted := (*aes.encryptor).EncryptNew(stateBatchedPlain)
		aes.inputEncrypted = append(aes.inputEncrypted, stateBatchedEncrypted)
	}
	if aes.inputEncrypted[0] == nil {panic("input is not stored in aesStruct")}
}

func (aes *AESCtr) EncodeCiphertext(ciphertexts []uint8, numBlock int) {
	aes.encodeCipher = make([]*ckks_fv.Ciphertext, aes.blockSize)
	if numBlock < aes.params.Slots() {
		fmt.Println("data is not full pack, fill with 0...")
	}
	encryptedData := make([]*BitSet, numBlock)
	for i, bit := range encryptedData {
		bit = NewBitSet(aes.blockSize)
		encryptedData[i] = bit 
	}

	for i := range encryptedData {
		if i < numBlock {
			if aes.allZeroIn {
				encryptedData[i].Set(0)
			} else {
				for k := 0; k < aes.blockSize && i*aes.blockSize+k < numBlock*aes.blockSize; k++ {
					ind := i*aes.blockSize + k
					bit := (ciphertexts[ind/8] >> uint(ind%8)) & 1
					encryptedData[i].bits[k] = uint8(bit)
				}
			}
		} else {
			encryptedData[i].Set(0)
		}
	}
	for i := 0; i < aes.blockSize; i++ {
		var dataBatched []complex128
		for j := 0; j < aes.params.Slots(); j++ {
			dataBatched = append(dataBatched, complex( float64( encryptedData[j].bits[i] ), 0.0 ) )
		}
		encryptedDataPlain := (*aes.encoder).EncodeComplexNTTAtLvlNew(aes.remainingLevel, dataBatched, aes.params.LogSlots())
		encryptedDataCtxt := (*aes.encryptor).EncryptNew(encryptedDataPlain)
		aes.encodeCipher = append(aes.encodeCipher, encryptedDataCtxt)
	}
}

func (aes *AESCtr) aesRoundFunction(state []*ckks_fv.Ciphertext, roundKey []*ckks_fv.Ciphertext) {
	fmt.Printf("Chain index before sbox: %d, scale: %f\n", state[0].Level(), math.Log2( state[0].Scale() ) )

	// SubByte
	for i := range state {
		currLevel := state[i].Level()
		for currLevel > 3 {
			aes.DropLevel(state[i], 1)
			currLevel--
		}
	}
	// // Parallel processing using goroutines
	// ch := make(chan bool, 16)
	// for i := 0; i < 16; i++ {
	// 	go func(i int) {
	// 		aes.aesSubbyteLUT(state[i*8 : (i+1)*8])
	// 		ch <- true
	// 	}(i)
	// }
	// for i := 0; i < 16; i++ {
	// 	<-ch
	// }
	for i := 0; i < 16; i++ {
		aes.aesSubbyteLUT(state[i*8 : (i+1)*8])
	}

	for i:=0;i<8;i++{
		aes.DebugPrint(state[i], "after Sbox")
	}

	start := time.Now()

	// // Parallel processing for bootstrapping and cleaning tensor
	// ch := make(chan bool, len(state))
	// for i := range state {
	// 	go func(i int) {
	// 		state[i] = aes.BootstrapCipher(state[i])
	// 		if i == 0 {
	// 			aes.DebugPrint(state[i], "BTS precise: ")
	// 		}
	// 		aes.CleanReal(state[i])
	// 		ch <- true
	// 	}(i)
	// }
	// for i := 0; i < len(state); i++ {
	// 	<-ch
	// }
	for i := range state {
		state[i] = aes.BootstrapCipher(state[i])
		if i == 0 {
			aes.DebugPrint(state[i], "BTS precise: ")
		}
		aes.CleanReal(state[i])
	}
	aes.DebugPrint(state[0], "after btp Clean: ")
	end := time.Now()
	timeBTSreEnc := end.Sub(start)
	// Debug output time
	fmt.Printf("Bootstrap and re-encryption time: %v\n", timeBTSreEnc)

	fmt.Printf("ShiftRow Chain: %d, scale: %f\n", state[0].Level(), math.Log2( state[0].Scale() ) )
	// ShiftRow
	aes.ShiftRow(state)
	fmt.Printf("MixColumn Chain: %d, scale: %f\n", state[0].Level(), math.Log2( state[0].Scale() ) )
	// MixColumn
	aes.MixColumn(state)
	fmt.Printf("AddRoundKey Chain: %d, scale: %f\n", state[0].Level(), math.Log2( state[0].Scale() ) )
	// AddRoundKey
	aes.AddRoundKey(state, roundKey)
}

func (aes *AESCtr) coefficientMultMonomial(mon []*ckks_fv.Ciphertext, coeffArr []int, pos int) (*ckks_fv.Ciphertext) {
    if len(mon)!= len(aes.sboxMonomialOrder) {
		panic("monomial size must equal to sbox_monomial_order!")
    }
    ctxtPool := []*ckks_fv.Ciphertext{}
    for i, m := range mon {
        ind := int(aes.sboxMonomialOrder[i].ToULong()) - 1
        coeff := coeffArr[ind]
        if coeff == 0 {
            continue
        } else if coeff > 0 {
            for j := 0; j < coeff; j++ {
                ctxtPool = append(ctxtPool, m)
            }
        } else {
            tmp := aes.NegNew(m)
            for j := 0; j < (-coeff); j++ {
                ctxtPool = append(ctxtPool, tmp)
            }
        }
    }
	for i:=1;i<len(ctxtPool);i++{
		aes.Add(ctxtPool[0], ctxtPool[i], ctxtPool[0])
	}
	// sum := aes.BinaryTreeAddNew(ctxtPool)
    return aes.AddConstNew(ctxtPool[0], int( aes.bitSbox[0].bits[pos] ) )
}

func (aes *AESCtr) aesSubbyteLUT(SBoxIn []*ckks_fv.Ciphertext) {
    // construct 8-bit val of the sbox
    if len(SBoxIn)!= 8 {
		panic("The input length of the Sbox is wrong (8bit)!!")
    }
    sboxMonomials, _ := aes.LayeredCombine(SBoxIn)
    SBoxIn[0] = aes.coefficientMultMonomial(sboxMonomials, Sbox0[:], 0)
    SBoxIn[1] = aes.coefficientMultMonomial(sboxMonomials, Sbox1[:], 1)
    SBoxIn[2] = aes.coefficientMultMonomial(sboxMonomials, Sbox2[:], 2)
    SBoxIn[3] = aes.coefficientMultMonomial(sboxMonomials, Sbox3[:], 3)
    SBoxIn[4] = aes.coefficientMultMonomial(sboxMonomials, Sbox4[:], 4)
    SBoxIn[5] = aes.coefficientMultMonomial(sboxMonomials, Sbox5[:], 5)
    SBoxIn[6] = aes.coefficientMultMonomial(sboxMonomials, Sbox6[:], 6)
    SBoxIn[7] = aes.coefficientMultMonomial(sboxMonomials, Sbox7[:], 7)
}

func (aes *AESCtr) GF2FieldMul(x []*ckks_fv.Ciphertext) {
    y := make([]*ckks_fv.Ciphertext, len(x))
    for i := 0; i < 4; i++ {
        y[0+8*i] = x[1+8*i]
        y[1+8*i] = x[2+8*i]
        y[2+8*i] = x[3+8*i]
        y[3+8*i] = aes.XORNew(x[4+8*i], x[0+8*i])
        y[4+8*i] = aes.XORNew(x[5+8*i], x[0+8*i])
        y[5+8*i] = x[6+8*i]
        y[7+8*i] = x[0+8*i]
		y[6+8*i] = aes.XORNew(x[7+8*i], x[0+8*i])
    }
	copy(x, y)
}

func (aes *AESCtr) MixColumn(x []*ckks_fv.Ciphertext) {
	x0, x1, x2, x3 := []*ckks_fv.Ciphertext{}, []*ckks_fv.Ciphertext{}, []*ckks_fv.Ciphertext{}, []*ckks_fv.Ciphertext{}
	
	for i := 0; i < 128; i++ {
		mod := i % 32
		if mod < 8 {
			x0 = append(x0, x[i])
		} else if mod >= 8 && mod < 16 {
			x1 = append(x1, x[i])
		} else if mod >= 16 && mod < 24 {
			x2 = append(x2, x[i])
		} else {
			x3 = append(x3, x[i])
		}
	}

	y0, y1, y2, y3 := make([]*ckks_fv.Ciphertext, len(x0)), make([]*ckks_fv.Ciphertext, len(x1)), make([]*ckks_fv.Ciphertext, len(x2)), make([]*ckks_fv.Ciphertext, len(x3))
	z0, z1, z2, z3 := make([]*ckks_fv.Ciphertext, len(x0)), make([]*ckks_fv.Ciphertext, len(x1)), make([]*ckks_fv.Ciphertext, len(x2)), make([]*ckks_fv.Ciphertext, len(x3))

	for i := 0; i < 32; i++ {
		y0[i] = aes.XORNew(x0[i], x1[i])
		y1[i] = aes.XORNew(x1[i], x2[i])
		y2[i] = aes.XORNew(x2[i], x3[i])
		y3[i] = aes.XORNew(x3[i], x0[i])
	}

	for i := 0; i < 32; i++ {
		z0[i] = aes.XORNew(y1[i], x3[i])
		z1[i] = aes.XORNew(y2[i], x0[i])
		z2[i] = aes.XORNew(y3[i], x1[i])
		z3[i] = aes.XORNew(y0[i], x2[i])
	}

	aes.GF2FieldMul(y0)
	aes.GF2FieldMul(y1)
	aes.GF2FieldMul(y2)
	aes.GF2FieldMul(y3)

	for i := 0; i < 32; i++ {
		z0[i] = aes.XORNew(z0[i], y0[i])
		z1[i] = aes.XORNew(z1[i], y1[i])
		z2[i] = aes.XORNew(z2[i], y2[i])
		z3[i] = aes.XORNew(z3[i], y3[i])
	}

	z0 = append(z0, z1...)
	z0 = append(z0, z2...)
	z0 = append(z0, z3...)
	copy(x, z0)
}
func (aes *AESCtr) ShiftRow(x []*ckks_fv.Ciphertext) {
	for i := 0; i < 8; i++ {
		x[1*8+i], x[5*8+i], x[9*8+i], x[13*8+i] = x[5*8+i], x[9*8+i], x[13*8+i], x[1*8+i]
		x[2*8+i], x[10*8+i], x[14*8+i], x[6*8+i] = x[10*8+i], x[14*8+i], x[6*8+i], x[2*8+i]
		x[3*8+i], x[15*8+i], x[11*8+i], x[7*8+i] = x[15*8+i], x[11*8+i], x[7*8+i], x[3*8+i]
	}
}

func (aes *AESCtr) AddWhiteKey(pt, key []*ckks_fv.Ciphertext) []*ckks_fv.Ciphertext {
	for i := 0; i < len(pt); i++ {
		aes.XOR(pt[i], key[i], pt[i])
	}
	return pt
}

func (aes *AESCtr) AddRoundKey(state, key []*ckks_fv.Ciphertext) {
	for i := 0; i < len(state); i++ {
		aes.XOR(state[i], key[i], state[i])
	}
}