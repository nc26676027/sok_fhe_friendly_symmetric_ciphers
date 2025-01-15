package cnn_infer

import "github.com/ldsec/lattigo/v2/ckks_fv"

func ProcessTranscipherResult( cipher *ckks_fv.Ciphertext, init_p, dataLenth, cur_pos int, context *testParams) (*ckks_fv.Ciphertext) {
	
	//Rotate and mask zero to get data toped
	zero_mask_top := make([]complex128, context.params.Slots())
	
	for i:=0; i<(context.params.Slots()/init_p);i++ {
		if( i < dataLenth ) {
			zero_mask_top[i] = complex(1.0, 0.0)
		} else {
			zero_mask_top[i] = complex(0.0, 0.0)
		} 
	}

	var rotate2top *ckks_fv.Ciphertext
	if (cur_pos != 0){
		rotate2top = memorySaveRotateCipher(cipher, context, (dataLenth*cur_pos), context.params.Slots())
	} else {
		rotate2top = cipher
	}
	
	zero_mask := ckks_fv.NewPlaintextCKKS(context.params, cipher.Level(), cipher.Scale() )
	context.encoder.EncodeComplexNTT( zero_mask, zero_mask_top, context.params.LogSlots() )
	context.evaluator.Mul(rotate2top, zero_mask, rotate2top)
	context.evaluator.RescaleMany(rotate2top, 1, rotate2top)
	
	var rotatedCipher []*ckks_fv.Ciphertext
	rotatedCipher = append(rotatedCipher, rotate2top)
	for i:=1; i<init_p; i++{
		tmp := memorySaveRotateCipher(rotate2top, context, (context.params.Slots()/init_p)*i, context.params.Slots())
		rotatedCipher = append(rotatedCipher, tmp)
	}
	for i:=1; i<len(rotatedCipher); i++{
		context.evaluator.Add(rotatedCipher[0], rotatedCipher[i], rotatedCipher[0])
	}
	return rotatedCipher[0]
}
