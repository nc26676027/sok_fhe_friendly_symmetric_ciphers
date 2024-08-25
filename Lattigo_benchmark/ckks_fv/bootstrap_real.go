package ckks_fv

import (
	"fmt"
	"math"
)

// NewBootstrapper creates a new Bootstrapper.
func NewBootstrapperReal(params *Parameters, btpParams *BootstrappingParameters, btpKey BootstrappingKey) (btp *Bootstrapper, err error) {

	if btpParams.SinType == SinType(Sin) && btpParams.SinRescal != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	btp = NewbootstrapperReal(params, btpParams)

	btp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = btp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	btp.ckksEvaluator = btp.ckksEvaluator.WithKey(EvaluationKey{btpKey.Rlk, btpKey.Rtks}).(*ckksEvaluator)

	return btp, nil
}

// newBootstrapper is a constructor of "dummy" bootstrapper to enable the generation of bootstrapping-related constants
// without providing a bootstrapping key. To be replaced by a propper factorization of the bootstrapping pre-computations.
func NewbootstrapperReal(params *Parameters, btpParams *BootstrappingParameters) (btp *Bootstrapper) {
	btp = new(Bootstrapper)

	btp.params = params.Copy()
	btp.BootstrappingParameters = *btpParams.Copy()

	// btp.dslots = params.Slots()
	// btp.logdslots = params.LogSlots()
	btp.dslots = 1<<btpParams.LogSlots
	btp.logdslots = btpParams.LogSlots
	// fmt.Println("logdslots: ", btp.logdslots)
	// if params.logSlots < params.MaxLogSlots() {
	if btpParams.LogSlots < params.MaxLogSlots() {
		btp.dslots <<= 1
		btp.logdslots++
	}

	btp.prescale = math.Exp2(math.Round(math.Log2(float64(params.qi[0]) / btp.MessageRatio)))
	btp.sinescale = math.Exp2(math.Round(math.Log2(btp.SineEvalModuli.ScalingFactor)))
	btp.postscale = btp.sinescale / btp.MessageRatio

	btp.encoder = NewCKKSEncoder(params)
	btp.ckksEvaluator = NewCKKSEvaluator(params, EvaluationKey{}).(*ckksEvaluator) // creates an evaluator without keys for genDFTMatrices

	btp.genSinePoly()
	btp.genDFTMatricesHalf()

	btp.ctxpool = NewCiphertextCKKS(params, 1, params.MaxLevel(), 0)

	return btp
}


// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller or equal to round(Q0/2^{10}).
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp *Bootstrapper) SlimBootstrappReal(ct *Ciphertext) *Ciphertext {

	//var t time.Time
	var ct0, ct1 *Ciphertext

	// Drops the level to 1
	for ct.Level() > 1 {
		btp.ckksEvaluator.DropLevel(ct, 1)
	}

	// Brings the ciphertext scale to Q0/2^{10}
	if ct.Level() == 1 {

		// if one level is available, then uses it to match the scale
		btp.ckksEvaluator.SetScale(ct, btp.prescale)

		// then drops to level 0
		for ct.Level() != 0 {
			btp.ckksEvaluator.DropLevel(ct, 1)
		}

	} else {
		// else drop to level 0
		for ct.Level() != 0 {
			btp.ckksEvaluator.DropLevel(ct, 1)
		}

		// and does an integer constant mult by round((Q0/Delta_m)/ctscle)

		if btp.prescale < ct.Scale() {
			panic("ciphetext scale > Q[0]/(Q[0]/Delta_m)")
		}
		btp.ckksEvaluator.ScaleUp(ct, math.Round(btp.prescale/ct.Scale()), ct)
	}
	ct0 = SlotsToCoeffsHalf(ct0, ct1, btp.pDFT, btp.ckksEvaluator)

	// ModUp ct_{Q_0} -> ct_{Q_L}
	//t = time.Now()
	ct = btp.modUp(ct)
	//log.Println("After ModUp  :", time.Now().Sub(t), ct.Level(), ct.Scale())

	// Brings the ciphertext scale to sineQi/(Q0/scale) if its under
	btp.ckksEvaluator.ScaleUp(ct, math.Round(btp.postscale/ct.Scale()), ct)

	//SubSum X -> (N/dslots) * Y^dslots
	//t = time.Now()
	ct = btp.subSum(ct)
	//log.Println("After SubSum :", time.Now().Sub(t), ct.Level(), ct.Scale())
	// Part 1 : Coeffs to slots

	//t = time.Now()
	ct0, ct1 = CoeffsToSlots(ct, btp.pDFTInv, btp.ckksEvaluator)
	//log.Println("After CtS    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 2 : SineEval
	//t = time.Now()
	ct0, ct1 = btp.evaluateSine(ct0, ct1)
	//log.Println("After Sine   :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 3 : Slots to coeffs
	//t = time.Now()

	

	ct0.SetScale( math.Exp2(math.Round(math.Log2(ct0.Scale()))) ) // rounds to the nearest power of two and scale*0.5 to compensate 2 factor introduced by conjct add
	//log.Println("After StC    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())
    // var conjct *Ciphertext
	conjct := btp.ckksEvaluator.ConjugateNew(ct0)
	btp.ckksEvaluator.Add(ct0, conjct, ct0)
	return ct0
}

// Bootstrapp re-encrypt a ciphertext at lvl Q0 to a ciphertext at MaxLevel-k where k is the depth of the bootstrapping circuit.
// If the input ciphertext level is zero, the input scale must be an exact power of two smaller or equal to round(Q0/2^{10}).
// If the input ciphertext is at level one or more, the input scale does not need to be an exact power of two as one level
// can be used to do a scale matching.
func (btp *Bootstrapper) BootstrappReal(ct *Ciphertext) *Ciphertext {

	//var t time.Time
	var ct0, ct1 *Ciphertext

	// Drops the level to 1
	for ct.Level() > 1 {
		btp.ckksEvaluator.DropLevel(ct, 1)
	}

	// Brings the ciphertext scale to Q0/2^{10}
	if ct.Level() == 1 {

		// if one level is available, then uses it to match the scale
		btp.ckksEvaluator.SetScale(ct, btp.prescale)

		// then drops to level 0
		for ct.Level() != 0 {
			btp.ckksEvaluator.DropLevel(ct, 1)
		}

	} else {

		// else drop to level 0
		for ct.Level() != 0 {
			btp.ckksEvaluator.DropLevel(ct, 1)
		}

		// and does an integer constant mult by round((Q0/Delta_m)/ctscle)

		if btp.prescale < ct.Scale() {
			panic("ciphetext scale > Q[0]/(Q[0]/Delta_m)")
		}
		btp.ckksEvaluator.ScaleUp(ct, math.Round(btp.prescale/ct.Scale()), ct)
	}

	// ModUp ct_{Q_0} -> ct_{Q_L}
	//t = time.Now()
	ct = btp.modUp(ct)
	//log.Println("After ModUp  :", time.Now().Sub(t), ct.Level(), ct.Scale())

	// Brings the ciphertext scale to sineQi/(Q0/scale) if its under
	btp.ckksEvaluator.ScaleUp(ct, math.Round(btp.postscale/ct.Scale()), ct)

	//SubSum X -> (N/dslots) * Y^dslots
	//t = time.Now()
	ct = btp.subSum(ct)
	//log.Println("After SubSum :", time.Now().Sub(t), ct.Level(), ct.Scale())
	// Part 1 : Coeffs to slots

	//t = time.Now()
	ct0, ct1 = CoeffsToSlots(ct, btp.pDFTInv, btp.ckksEvaluator)
	//log.Println("After CtS    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 2 : SineEval
	//t = time.Now()
	ct0, ct1 = btp.evaluateSine(ct0, ct1)
	//log.Println("After Sine   :", time.Now().Sub(t), ct0.Level(), ct0.Scale())

	// Part 3 : Slots to coeffs
	//t = time.Now()

	ct0 = SlotsToCoeffsHalf(ct0, ct1, btp.pDFT, btp.ckksEvaluator)

	ct0.SetScale( math.Exp2(math.Round(math.Log2(ct0.Scale()))) ) // rounds to the nearest power of two and scale*0.5 to compensate 2 factor introduced by conjct add
	//log.Println("After StC    :", time.Now().Sub(t), ct0.Level(), ct0.Scale())
    // var conjct *Ciphertext
	conjct := btp.ckksEvaluator.ConjugateNew(ct0)
	btp.ckksEvaluator.Add(ct0, conjct, ct0)
	return ct0
}

func SlotsToCoeffsHalf(ct0, ct1 *Ciphertext, pDFT []*PtDiagMatrix, eval CKKSEvaluator) (ct *Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ct1 != nil {
		eval.MultByi(ct1, ct1)
		eval.Add(ct0, ct1, ct0)
	}

	ct1 = nil

	return dftHalf(ct0, pDFT, false, eval)
}

func SlimSlotsToCoeffs(ct0, ct1 *Ciphertext, pDFT []*PtDiagMatrix, eval CKKSEvaluator) (ct *Ciphertext) {

	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ct1 != nil {
		eval.MultByi(ct1, ct1)
		eval.Add(ct0, ct1, ct0)
	}

	ct1 = nil

	return dftHalf(ct0, pDFT, false, eval)
}

func dftHalf(vec *Ciphertext, plainVectors []*PtDiagMatrix, forward bool, eval CKKSEvaluator) *Ciphertext {

	// Sequentially multiplies w with the provided dft matrices.
	for _, plainVector := range plainVectors {
		scale := vec.Scale()
		vec = eval.LinearTransform(vec, plainVector)[0]
		if err := eval.Rescale(vec, scale, vec); err != nil {
			panic(err)
		}
	}

	return vec
}

