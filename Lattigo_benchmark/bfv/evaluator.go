package bfv

import (
	"fmt"
	"math/big"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/rlwe"
	"github.com/ldsec/lattigo/v2/utils"

	"unsafe"
)

// Evaluator is an interface implementing the public methodes of the eval.
type Evaluator interface {
	Add(op0, op1 Operand, ctOut *Ciphertext)
	AddNew(op0, op1 Operand) (ctOut *Ciphertext)
	AddNoMod(op0, op1 Operand, ctOut *Ciphertext)
	AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext)
	Sub(op0, op1 Operand, ctOut *Ciphertext)
	SubNew(op0, op1 Operand) (ctOut *Ciphertext)
	SubNoMod(op0, op1 Operand, ctOut *Ciphertext)
	SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext)
	Neg(op Operand, ctOut *Ciphertext)
	NegNew(op Operand) (ctOut *Ciphertext)
	Reduce(op Operand, ctOut *Ciphertext)
	ReduceNew(op Operand) (ctOut *Ciphertext)
	MulScalar(op Operand, scalar uint64, ctOut *Ciphertext)
	MulScalarNew(op Operand, scalar uint64) (ctOut *Ciphertext)
	Mul(op0 *Ciphertext, op1 Operand, ctOut *Ciphertext)
	MulNew(op0 *Ciphertext, op1 Operand) (ctOut *Ciphertext)
	Relinearize(ct0 *Ciphertext, ctOut *Ciphertext)
	RelinearizeNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	SwitchKeys(ct0 *Ciphertext, switchKey *SwitchingKey, ctOut *Ciphertext)
	SwitchKeysNew(ct0 *Ciphertext, switchkey *SwitchingKey) (ctOut *Ciphertext)
	RotateColumnsNew(ct0 *Ciphertext, k int) (ctOut *Ciphertext)
	RotateColumns(ct0 *Ciphertext, k int, ctOut *Ciphertext)
	RotateRows(ct0 *Ciphertext, ctOut *Ciphertext)
	RotateRowsNew(ct0 *Ciphertext) (ctOut *Ciphertext)
	InnerSum(ct0 *Ciphertext, ctOut *Ciphertext)
	ShallowCopy() Evaluator
	WithKey(EvaluationKey) Evaluator
}

// evaluator is a struct that holds the necessary elements to perform the homomorphic operations between ciphertexts and/or plaintexts.
// It also holds a small memory pool used to store intermediate computations.
type evaluator struct {
	*evaluatorBase
	*evaluatorBuffers

	rlk  *RelinearizationKey
	rtks *RotationKeySet

	baseconverterQ1Q2 *ring.FastBasisExtender
	baseconverterQ1P  *ring.FastBasisExtender
}

type evaluatorBase struct {
	params   *Parameters
	ringQ    *ring.Ring
	ringP    *ring.Ring
	ringQMul *ring.Ring

	decomposer *ring.Decomposer

	t     uint64
	pHalf *big.Int

	deltaMont []uint64
}

func newEvaluatorPrecomp(params *Parameters) *evaluatorBase {
	var err error
	ev := new(evaluatorBase)

	ev.params = params.Copy()

	ev.t = params.t

	if ev.ringQ, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	// Generates #QiMul primes such that Q * QMul > Q*Q*N
	logQTimesN := ev.ringQ.ModulusBigint.BitLen() + params.LogN()
	var nbQiMul, logQMul int
	for logQMul < logQTimesN {
		nbQiMul++
		logQMul += 61
	}

	qiMul := ring.GenerateNTTPrimesP(61, 2*params.N(), nbQiMul)

	if ev.ringQMul, err = ring.NewRing(params.N(), qiMul); err != nil {
		panic(err)
	}

	ev.pHalf = new(big.Int).Rsh(ev.ringQMul.ModulusBigint, 1)
	ev.deltaMont = GenLiftParams(ev.ringQ, params.t)

	if len(params.pi) != 0 {

		if ev.ringP, err = ring.NewRing(params.N(), params.pi); err != nil {
			panic(err)
		}
		ev.decomposer = ring.NewDecomposer(ev.ringQ.Modulus, ev.ringP.Modulus)
	}
	return ev
}

type evaluatorBuffers struct {
	poolQ    [][]*ring.Poly
	poolQmul [][]*ring.Poly

	poolQKS [4]*ring.Poly
	poolPKS [3]*ring.Poly

	tmpPt *Plaintext
}

func newEvaluatorBuffer(eval *evaluatorBase) *evaluatorBuffers {
	evb := new(evaluatorBuffers)
	evb.poolQ = make([][]*ring.Poly, 4)
	evb.poolQmul = make([][]*ring.Poly, 4)
	for i := 0; i < 4; i++ {
		evb.poolQ[i] = make([]*ring.Poly, 6)
		evb.poolQmul[i] = make([]*ring.Poly, 6)
		for j := 0; j < 6; j++ {
			evb.poolQ[i][j] = eval.ringQ.NewPoly()
			evb.poolQmul[i][j] = eval.ringQMul.NewPoly()
		}
	}
	if eval.ringP != nil {
		evb.poolQKS = [4]*ring.Poly{eval.ringQ.NewPoly(), eval.ringQ.NewPoly(), eval.ringQ.NewPoly(), eval.ringQ.NewPoly()}
		evb.poolPKS = [3]*ring.Poly{eval.ringP.NewPoly(), eval.ringP.NewPoly(), eval.ringP.NewPoly()}
	}

	evb.tmpPt = NewPlaintext(eval.params)

	return evb
}

// NewEvaluator creates a new Evaluator, that can be used to do homomorphic
// operations on ciphertexts and/or plaintexts. It stores a small pool of polynomials
// and ciphertexts that will be used for intermediate values.
func NewEvaluator(params *Parameters, evaluationKey EvaluationKey) Evaluator {
	ev := new(evaluator)
	ev.evaluatorBase = newEvaluatorPrecomp(params)
	ev.evaluatorBuffers = newEvaluatorBuffer(ev.evaluatorBase)
	ev.baseconverterQ1Q2 = ring.NewFastBasisExtender(ev.ringQ, ev.ringQMul)
	if len(params.pi) != 0 {
		ev.baseconverterQ1P = ring.NewFastBasisExtender(ev.ringQ, ev.ringP)
	}
	ev.rlk = evaluationKey.Rlk
	ev.rtks = evaluationKey.Rtks
	return ev
}

// NewEvaluators creates n evaluators sharing the same read-only data-structures.
func NewEvaluators(params *Parameters, evaluationKey EvaluationKey, n int) []Evaluator {
	if n <= 0 {
		return []Evaluator{}
	}
	evas := make([]Evaluator, n, n)
	for i := range evas {
		if i == 0 {
			evas[0] = NewEvaluator(params, evaluationKey)
		} else {
			evas[i] = evas[i-1].ShallowCopy()
		}
	}
	return evas
}

// ShallowCopy creates a shallow copy of this evaluator in which the read-only data-structures are
// shared with the receiver.
func (eval *evaluator) ShallowCopy() Evaluator {
	return &evaluator{
		evaluatorBase:     eval.evaluatorBase,
		evaluatorBuffers:  newEvaluatorBuffer(eval.evaluatorBase),
		baseconverterQ1Q2: eval.baseconverterQ1Q2.ShallowCopy(),
		baseconverterQ1P:  eval.baseconverterQ1P.ShallowCopy(),
		rlk:               eval.rlk,
		rtks:              eval.rtks,
	}
}

// ShallowCopyWithKey creates a shallow copy of this evaluator in which the read-only data-structures are
// shared with the receiver but the EvaluationKey is evaluationKey.
func (eval *evaluator) WithKey(evaluationKey EvaluationKey) Evaluator {
	return &evaluator{
		evaluatorBase:     eval.evaluatorBase,
		evaluatorBuffers:  eval.evaluatorBuffers,
		baseconverterQ1Q2: eval.baseconverterQ1Q2,
		baseconverterQ1P:  eval.baseconverterQ1P,
		rlk:               evaluationKey.Rlk,
		rtks:              evaluationKey.Rtks,
	}
}

// Add adds op0 to op1 and returns the result in ctOut.
func (eval *evaluator) Add(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxInt(op0.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.Add)
}

// AddNew adds op0 to op1 and creates a new element ctOut to store the result.
func (eval *evaluator) AddNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, utils.MaxInt(op0.Degree(), op1.Degree()))
	eval.Add(op0, op1, ctOut)
	return
}

// AddNoMod adds op0 to op1 without modular reduction, and returns the result in cOut.
func (eval *evaluator) AddNoMod(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxInt(op0.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.AddNoMod)
}

// AddNoModNew adds op0 to op1 without modular reduction and creates a new element ctOut to store the result.
func (eval *evaluator) AddNoModNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, utils.MaxInt(op0.Degree(), op1.Degree()))
	eval.AddNoMod(op0, op1, ctOut)
	return
}

// Sub subtracts op1 from op0 and returns the result in cOut.
func (eval *evaluator) Sub(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxInt(op0.Degree(), op1.Degree()), true)
	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.Sub)

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.Neg(ctOut.Value()[i], ctOut.Value()[i])
		}
	}
}

// SubNew subtracts op1 from op0 and creates a new element ctOut to store the result.
func (eval *evaluator) SubNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, utils.MaxInt(op0.Degree(), op1.Degree()))
	eval.Sub(op0, op1, ctOut)
	return
}

// SubNoMod subtracts op1 from op0 without modular reduction and returns the result on ctOut.
func (eval *evaluator) SubNoMod(op0, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, utils.MaxInt(op0.Degree(), op1.Degree()), true)

	eval.evaluateInPlaceBinary(el0, el1, elOut, eval.ringQ.SubNoMod)

	if el0.Degree() < el1.Degree() {
		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			eval.ringQ.Neg(ctOut.Value()[i], ctOut.Value()[i])
		}
	}
}

// SubNoModNew subtracts op1 from op0 without modular reduction and creates a new element ctOut to store the result.
func (eval *evaluator) SubNoModNew(op0, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, utils.MaxInt(op0.Degree(), op1.Degree()))
	eval.SubNoMod(op0, op1, ctOut)
	return
}

// Neg negates op and returns the result in ctOut.
func (eval *evaluator) Neg(op Operand, ctOut *Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(op, ctOut, op.Degree())
	evaluateInPlaceUnary(el0, elOut, eval.ringQ.Neg)
}

// NegNew negates op and creates a new element to store the result.
func (eval *evaluator) NegNew(op Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, op.Degree())
	eval.Neg(op, ctOut)
	return ctOut
}

// Reduce applies a modular reduction to op and returns the result in ctOut.
func (eval *evaluator) Reduce(op Operand, ctOut *Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(op, ctOut, op.Degree())
	evaluateInPlaceUnary(el0, elOut, eval.ringQ.Reduce)
}

// ReduceNew applies a modular reduction to op and creates a new element ctOut to store the result.
func (eval *evaluator) ReduceNew(op Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, op.Degree())
	eval.Reduce(op, ctOut)
	return ctOut
}

// MulScalar multiplies op by a uint64 scalar and returns the result in ctOut.
func (eval *evaluator) MulScalar(op Operand, scalar uint64, ctOut *Ciphertext) {
	el0, elOut := eval.getElemAndCheckUnary(op, ctOut, op.Degree())
	fun := func(el, elOut *ring.Poly) { eval.ringQ.MulScalar(el, scalar, elOut) }
	evaluateInPlaceUnary(el0, elOut, fun)
}

// MulScalarNew multiplies op by a uint64 scalar and creates a new element ctOut to store the result.
func (eval *evaluator) MulScalarNew(op Operand, scalar uint64) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, op.Degree())
	eval.MulScalar(op, scalar, ctOut)
	return
}

// tensorAndRescale computes (ct0 x ct1) * (t/Q) and stores the result in ctOut.
func (eval *evaluator) tensorAndRescale(ct0, ct1, ctOut *Element) {

	c0Q1 := eval.poolQ[0]
	c0Q2 := eval.poolQmul[0]

	c1Q1 := eval.poolQ[1]
	c1Q2 := eval.poolQmul[1]

	// Prepares the ciphertexts for the Tensoring by extending their
	// basis from Q to QP and transforming them to NTT form
	eval.modUpAndNTT(ct0, c0Q1, c0Q2)

	if ct0 != ct1 {
		eval.modUpAndNTT(ct1, c1Q1, c1Q2)
	}

	// Tensoring: multiplies each elements of the ciphertexts together
	// and adds them to their corresponding position in the new ciphertext
	// based on their respective degree

	// Case where both Elements are of degree 1
	if ct0.Degree() == 1 && ct1.Degree() == 1 {
		eval.tensoreLowDeg(ct0, ct1)
		// Case where at least one element is not of degree 1
	} else {
		eval.tensortLargeDeg(ct0, ct1)
	}

	eval.quantize(ctOut)
}

func (eval *evaluator) modUpAndNTT(ct *Element, cQ, cQMul []*ring.Poly) {
	levelQ := len(eval.ringQ.Modulus) - 1
	for i := range ct.value {
		eval.baseconverterQ1Q2.ModUpSplitQP(levelQ, ct.value[i], cQMul[i])
		eval.ringQ.NTTLazy(ct.value[i], cQ[i])
		eval.ringQMul.NTTLazy(cQMul[i], cQMul[i])
	}
}

func (eval *evaluator) tensoreLowDeg(ct0, ct1 *Element) {

	c0Q1 := eval.poolQ[0]
	c0Q2 := eval.poolQmul[0]

	c1Q1 := eval.poolQ[1]
	c1Q2 := eval.poolQmul[1]

	c2Q1 := eval.poolQ[2]
	c2Q2 := eval.poolQmul[2]

	c00Q := eval.poolQ[3][0]
	c00Q2 := eval.poolQmul[3][0]
	c01Q := eval.poolQ[3][1]
	c01P := eval.poolQmul[3][1]

	eval.ringQ.MForm(c0Q1[0], c00Q)
	eval.ringQMul.MForm(c0Q2[0], c00Q2)

	eval.ringQ.MForm(c0Q1[1], c01Q)
	eval.ringQMul.MForm(c0Q2[1], c01P)

	// Squaring case
	if ct0 == ct1 {

		// c0 = c0[0]*c0[0]
		eval.ringQ.MulCoeffsMontgomery(c00Q, c0Q1[0], c2Q1[0])
		eval.ringQMul.MulCoeffsMontgomery(c00Q2, c0Q2[0], c2Q2[0])

		// c1 = 2*c0[0]*c0[1]
		eval.ringQ.MulCoeffsMontgomery(c00Q, c0Q1[1], c2Q1[1])
		eval.ringQMul.MulCoeffsMontgomery(c00Q2, c0Q2[1], c2Q2[1])

		eval.ringQ.AddNoMod(c2Q1[1], c2Q1[1], c2Q1[1])
		eval.ringQMul.AddNoMod(c2Q2[1], c2Q2[1], c2Q2[1])

		// c2 = c0[1]*c0[1]
		eval.ringQ.MulCoeffsMontgomery(c01Q, c0Q1[1], c2Q1[2])
		eval.ringQMul.MulCoeffsMontgomery(c01P, c0Q2[1], c2Q2[2])

		// Normal case
	} else {

		// c0 = c0[0]*c1[0]
		eval.ringQ.MulCoeffsMontgomery(c00Q, c1Q1[0], c2Q1[0])
		eval.ringQMul.MulCoeffsMontgomery(c00Q2, c1Q2[0], c2Q2[0])

		// c1 = c0[0]*c1[1] + c0[1]*c1[0]
		eval.ringQ.MulCoeffsMontgomery(c00Q, c1Q1[1], c2Q1[1])
		eval.ringQMul.MulCoeffsMontgomery(c00Q2, c1Q2[1], c2Q2[1])

		eval.ringQ.MulCoeffsMontgomeryAndAddNoMod(c01Q, c1Q1[0], c2Q1[1])
		eval.ringQMul.MulCoeffsMontgomeryAndAddNoMod(c01P, c1Q2[0], c2Q2[1])

		// c2 = c0[1]*c1[1]
		eval.ringQ.MulCoeffsMontgomery(c01Q, c1Q1[1], c2Q1[2])
		eval.ringQMul.MulCoeffsMontgomery(c01P, c1Q2[1], c2Q2[2])
	}
}

func (eval *evaluator) tensortLargeDeg(ct0, ct1 *Element) {

	c0Q1 := eval.poolQ[0]
	c0Q2 := eval.poolQmul[0]

	c1Q1 := eval.poolQ[1]
	c1Q2 := eval.poolQmul[1]

	c2Q1 := eval.poolQ[2]
	c2Q2 := eval.poolQmul[2]

	for i := 0; i < ct0.Degree()+ct1.Degree()+1; i++ {
		c2Q1[i].Zero()
		c2Q2[i].Zero()
	}

	// Squaring case
	if ct0 == ct1 {

		c00Q1 := eval.poolQ[3]
		c00Q2 := eval.poolQmul[3]

		for i := range ct0.value {
			eval.ringQ.MForm(c0Q1[i], c00Q1[i])
			eval.ringQMul.MForm(c0Q2[i], c00Q2[i])
		}

		for i := 0; i < ct0.Degree()+1; i++ {
			for j := i + 1; j < ct0.Degree()+1; j++ {
				eval.ringQ.MulCoeffsMontgomery(c00Q1[i], c0Q1[j], c2Q1[i+j])
				eval.ringQMul.MulCoeffsMontgomery(c00Q2[i], c0Q2[j], c2Q2[i+j])

				eval.ringQ.Add(c2Q1[i+j], c2Q1[i+j], c2Q1[i+j])
				eval.ringQMul.Add(c2Q2[i+j], c2Q2[i+j], c2Q2[i+j])
			}
		}

		for i := 0; i < ct0.Degree()+1; i++ {
			eval.ringQ.MulCoeffsMontgomeryAndAdd(c00Q1[i], c0Q1[i], c2Q1[i<<1])
			eval.ringQMul.MulCoeffsMontgomeryAndAdd(c00Q2[i], c0Q2[i], c2Q2[i<<1])
		}

		// Normal case
	} else {
		for i := range ct0.value {
			eval.ringQ.MForm(c0Q1[i], c0Q1[i])
			eval.ringQMul.MForm(c0Q2[i], c0Q2[i])
			for j := range ct1.value {
				eval.ringQ.MulCoeffsMontgomeryAndAdd(c0Q1[i], c1Q1[j], c2Q1[i+j])
				eval.ringQMul.MulCoeffsMontgomeryAndAdd(c0Q2[i], c1Q2[j], c2Q2[i+j])
			}
		}
	}
}

func (eval *evaluator) quantize(ctOut *Element) {

	levelQ := len(eval.ringQ.Modulus) - 1
	levelQMul := len(eval.ringQMul.Modulus) - 1

	c2Q1 := eval.poolQ[2]
	c2Q2 := eval.poolQmul[2]

	// Applies the inverse NTT to the ciphertext, scales down the ciphertext
	// by t/q and reduces its basis from QP to Q
	for i := range ctOut.value {
		eval.ringQ.InvNTTLazy(c2Q1[i], c2Q1[i])
		eval.ringQMul.InvNTTLazy(c2Q2[i], c2Q2[i])

		// Extends the basis Q of ct(x) to the basis P and Divides (ct(x)Q -> P) by Q
		eval.baseconverterQ1Q2.ModDownSplitQP(levelQ, levelQMul, c2Q1[i], c2Q2[i], c2Q2[i])

		// Centers (ct(x)Q -> P)/Q by (P-1)/2 and extends ((ct(x)Q -> P)/Q) to the basis Q
		eval.ringQMul.AddScalarBigint(c2Q2[i], eval.pHalf, c2Q2[i])
		eval.baseconverterQ1Q2.ModUpSplitPQ(levelQMul, c2Q2[i], ctOut.value[i])
		eval.ringQ.SubScalarBigint(ctOut.value[i], eval.pHalf, ctOut.value[i])

		// Option (2) (ct(x)/Q)*T, doing so only requires that Q*P > Q*Q, faster but adds error ~|T|
		eval.ringQ.MulScalar(ctOut.value[i], eval.t, ctOut.value[i])
	}
}

// Mul multiplies op0 by op1 and returns the result in ctOut.
func (eval *evaluator) Mul(op0 *Ciphertext, op1 Operand, ctOut *Ciphertext) {
	el0, el1, elOut := eval.getElemAndCheckBinary(op0, op1, ctOut, op0.Degree()+op1.Degree(), false)
	switch op1 := op1.(type) {
	case *PlaintextMul:
		eval.mulPlaintextMul(op0, op1, ctOut)
	case *PlaintextRingT:
		eval.mulPlaintextRingT(op0, op1, ctOut)
	case *Plaintext, *Ciphertext:
		eval.tensorAndRescale(el0, el1, elOut)
	default:
		panic(fmt.Errorf("invalid operand type for Mul: %T", op1))
	}

}

func (eval *evaluator) mulPlaintextMul(ct0 *Ciphertext, ptRt *PlaintextMul, ctOut *Ciphertext) {
	for i := range ct0.value {
		eval.ringQ.NTTLazy(ct0.value[i], ctOut.value[i])
		eval.ringQ.MulCoeffsMontgomeryConstant(ctOut.value[i], ptRt.value, ctOut.value[i])
		eval.ringQ.InvNTT(ctOut.value[i], ctOut.value[i])
	}
}

func (eval *evaluator) mulPlaintextRingT(ct0 *Ciphertext, ptRt *PlaintextRingT, ctOut *Ciphertext) {
	ringQ := eval.ringQ

	coeffs := ptRt.value.Coeffs[0]
	coeffsNTT := eval.poolQ[0][0].Coeffs[0]

	for i := range ct0.value {

		// Copies the inputCT on the outputCT and switches to the NTT domain
		eval.ringQ.NTTLazy(ct0.value[i], ctOut.value[i])

		// Switches the outputCT in the Montgomery domain
		eval.ringQ.MForm(ctOut.value[i], ctOut.value[i])

		// For each qi in Q
		for j := range ringQ.Modulus {

			tmp := ctOut.value[i].Coeffs[j]
			qi := ringQ.Modulus[j]
			nttPsi := ringQ.NttPsi[j]
			bredParams := ringQ.BredParams[j]
			mredParams := ringQ.MredParams[j]

			// Transforms the plaintext in the NTT domain of that qi
			ring.NTTLazy(coeffs, coeffsNTT, ringQ.N, nttPsi, qi, mredParams, bredParams)

			// Multiplies NTT_qi(pt) * NTT_qi(ct)
			for k := 0; k < eval.ringQ.N; k = k + 8 {

				x := (*[8]uint64)(unsafe.Pointer(&coeffsNTT[k]))
				z := (*[8]uint64)(unsafe.Pointer(&tmp[k]))

				z[0] = ring.MRed(z[0], x[0], qi, mredParams)
				z[1] = ring.MRed(z[1], x[1], qi, mredParams)
				z[2] = ring.MRed(z[2], x[2], qi, mredParams)
				z[3] = ring.MRed(z[3], x[3], qi, mredParams)
				z[4] = ring.MRed(z[4], x[4], qi, mredParams)
				z[5] = ring.MRed(z[5], x[5], qi, mredParams)
				z[6] = ring.MRed(z[6], x[6], qi, mredParams)
				z[7] = ring.MRed(z[7], x[7], qi, mredParams)
			}
		}

		// Switches the ciphertext out of the NTT domain
		eval.ringQ.InvNTT(ctOut.value[i], ctOut.value[i])
	}
}

// MulNew multiplies op0 by op1 and creates a new element ctOut to store the result.
func (eval *evaluator) MulNew(op0 *Ciphertext, op1 Operand) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, op0.Degree()+op1.Degree())
	eval.Mul(op0, op1, ctOut)
	return
}

// relinearize is a method common to Relinearize and RelinearizeNew. It switches ct0 to the NTT domain, applies the keyswitch, and returns the result out of the NTT domain.
func (eval *evaluator) relinearize(ct0 *Ciphertext, ctOut *Ciphertext) {

	if ctOut != ct0 {
		eval.ringQ.Copy(ct0.value[0], ctOut.value[0])
		eval.ringQ.Copy(ct0.value[1], ctOut.value[1])
	}

	for deg := uint64(ct0.Degree()); deg > 1; deg-- {
		eval.switchKeysInPlace(ct0.value[deg], eval.rlk.Keys[deg-2], eval.poolQKS[1], eval.poolQKS[2])
		eval.ringQ.Add(ctOut.value[0], eval.poolQKS[1], ctOut.value[0])
		eval.ringQ.Add(ctOut.value[1], eval.poolQKS[2], ctOut.value[1])
	}

	ctOut.SetValue(ctOut.value[:2])
}

// Relinearize relinearizes the ciphertext ct0 of degree > 1 until it is of degree 1, and returns the result in cOut.
//
// It requires a correct evaluation key as additional input:
//
// - it must match the secret-key that was used to create the public key under which the current ct0 is encrypted.
//
// - it must be of degree high enough to relinearize the input ciphertext to degree 1 (e.g., a ciphertext
// of degree 3 will require that the evaluation key stores the keys for both degree 3 and degree 2 ciphertexts).
func (eval *evaluator) Relinearize(ct0 *Ciphertext, ctOut *Ciphertext) {

	if eval.rlk == nil {
		panic("evaluator has no relinearization key")
	}

	if ct0.Degree()-1 > len(eval.rlk.Keys) {
		panic("input ciphertext degree is too large to allow relinearization with the evluator's relinearization key")
	}

	if ct0.Degree() < 2 {
		if ct0 != ctOut {
			ctOut.Copy(ct0.El())
		}
	} else {
		eval.relinearize(ct0, ctOut)
	}
}

// RelinearizeNew relinearizes the ciphertext ct0 of degree > 1 until it is of degree 1, and creates a new ciphertext to store the result.
//
// Requires a correct evaluation key as additional input:
//
// - it must match the secret-key that was used to create the public key under which the current ct0 is encrypted
//
// - it must be of degree high enough to relinearize the input ciphertext to degree 1 (e.g., a ciphertext
// of degree 3 will require that the evaluation key stores the keys for both degree 3 and degree 2 ciphertexts).
func (eval *evaluator) RelinearizeNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1)
	eval.Relinearize(ct0, ctOut)
	return
}

// SwitchKeys applies the key-switching procedure to the ciphertext ct0 and returns the result in ctOut. It requires as an additional input a valid switching-key:
// it must encrypt the target key under the public key under which ct0 is currently encrypted.
func (eval *evaluator) SwitchKeys(ct0 *Ciphertext, switchKey *SwitchingKey, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot SwitchKeys: input and output must be of degree 1 to allow key switching")
	}

	eval.switchKeysInPlace(ct0.value[1], &switchKey.SwitchingKey, eval.poolQKS[1], eval.poolQKS[2])

	eval.ringQ.Add(ct0.value[0], eval.poolQKS[1], ctOut.value[0])
	eval.ringQ.Copy(eval.poolQKS[2], ctOut.value[1])
}

// SwitchKeysNew applies the key-switching procedure to the ciphertext ct0 and creates a new ciphertext to store the result. It requires as an additional input a valid switching-key:
// it must encrypt the target key under the public key under which ct0 is currently encrypted.
func (eval *evaluator) SwitchKeysNew(ct0 *Ciphertext, switchkey *SwitchingKey) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1)
	eval.SwitchKeys(ct0, switchkey, ctOut)
	return
}

// RotateColumns rotates the columns of ct0 by k positions to the left and returns the result in ctOut. As an additional input it requires a RotationKeys struct:
//
// - it must either store all the left and right power-of-2 rotations or the specific rotation that is requested.
//
// If only the power-of-two rotations are stored, the numbers k and n/2-k will be decomposed in base-2 and the rotation with the lowest
// hamming weight will be chosen; then the specific rotation will be computed as a sum of powers of two rotations.
func (eval *evaluator) RotateColumns(ct0 *Ciphertext, k int, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot RotateColumns: input and or output must be of degree 1")
	}

	if k == 0 {

		ctOut.Copy(ct0.El())

	} else {

		galElL := eval.params.GaloisElementForColumnRotationBy(k)
		// Looks in the rotation key if the corresponding rotation has been generated or if the input is a plaintext
		if swk, inSet := eval.rtks.GetRotationKey(galElL); inSet {

			eval.permute(ct0, galElL, swk, ctOut)

		} else {
			panic(fmt.Errorf("evaluator has no rotation key for rotation by %d", k))
		}
	}
}

// RotateColumnsNew applies RotateColumns and returns the result in a new Ciphertext.
func (eval *evaluator) RotateColumnsNew(ct0 *Ciphertext, k int) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1)
	eval.RotateColumns(ct0, k, ctOut)
	return
}

// RotateRows rotates the rows of ct0 and returns the result in ctOut.
func (eval *evaluator) RotateRows(ct0 *Ciphertext, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot RotateRows: input and/or output must be of degree 1")
	}

	galEl := eval.params.GaloisElementForRowRotation()

	if key, inSet := eval.rtks.GetRotationKey(galEl); inSet {
		eval.permute(ct0, galEl, key, ctOut)
	} else {
		panic("evaluator has no rotation key for row rotation")
	}
}

// RotateRowsNew rotates the rows of ct0 and returns the result a new Ciphertext.
func (eval *evaluator) RotateRowsNew(ct0 *Ciphertext) (ctOut *Ciphertext) {
	ctOut = NewCiphertext(eval.params, 1)
	eval.RotateRows(ct0, ctOut)
	return
}

// InnerSum computes the inner sum of ct0 and returns the result in ctOut. It requires a rotation key storing all the left powers of two rotations.
// The resulting vector will be of the form [sum, sum, .., sum, sum].
func (eval *evaluator) InnerSum(ct0 *Ciphertext, ctOut *Ciphertext) {

	if ct0.Degree() != 1 || ctOut.Degree() != 1 {
		panic("cannot InnerSum: input and output must be of degree 1")
	}

	cTmp := NewCiphertext(eval.params, 1)

	ctOut.Copy(ct0.El())

	for i := 1; i < int(eval.ringQ.N>>1); i <<= 1 {
		eval.RotateColumns(ctOut, i, cTmp)
		eval.Add(cTmp, ctOut, ctOut.Ciphertext())
	}

	eval.RotateRows(ctOut, cTmp)
	eval.Add(ctOut, cTmp, ctOut)
}

// permute performs a column rotation on ct0 and returns the result in ctOut
func (eval *evaluator) permute(ct0 *Ciphertext, generator uint64, switchKey *rlwe.SwitchingKey, ctOut *Ciphertext) {

	eval.switchKeysInPlace(ct0.value[1], switchKey, eval.poolQKS[1], eval.poolQKS[2])

	eval.ringQ.Add(eval.poolQKS[1], ct0.value[0], eval.poolQKS[1])

	eval.ringQ.Permute(eval.poolQKS[1], generator, ctOut.value[0])
	eval.ringQ.Permute(eval.poolQKS[2], generator, ctOut.value[1])
}

// switchKeys applies the general key-switching procedure of the form [c0 + cx*evakey[0], c1 + cx*evakey[1]]
func (eval *evaluator) switchKeysInPlace(cx *ring.Poly, evakey *rlwe.SwitchingKey, pool2Q, pool3Q *ring.Poly) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	pool2P := eval.poolPKS[1]
	pool3P := eval.poolPKS[2]

	level := len(ringQ.Modulus) - 1

	c2QiQ := eval.poolQKS[0]
	c2QiP := eval.poolPKS[0]
	c2 := eval.poolQKS[3]

	evakey0Q := new(ring.Poly)
	evakey1Q := new(ring.Poly)
	evakey0P := new(ring.Poly)
	evakey1P := new(ring.Poly)

	// We switch the element on which the key-switching operation will be conducted out of the NTT domain
	ringQ.NTTLazy(cx, c2)

	var reduce int

	// Key switching with CRT decomposition for the Qi
	for i := 0; i < eval.params.Beta(); i++ {

		eval.decomposeAndSplitNTT(level, i, c2, cx, c2QiQ, c2QiP)

		evakey0Q.Coeffs = evakey.Value[i][0].Coeffs[:level+1]
		evakey1Q.Coeffs = evakey.Value[i][1].Coeffs[:level+1]
		evakey0P.Coeffs = evakey.Value[i][0].Coeffs[level+1:]
		evakey1P.Coeffs = evakey.Value[i][1].Coeffs[level+1:]

		if i == 0 {
			ringQ.MulCoeffsMontgomeryLvl(level, evakey0Q, c2QiQ, pool2Q)
			ringQ.MulCoeffsMontgomeryLvl(level, evakey1Q, c2QiQ, pool3Q)
			ringP.MulCoeffsMontgomery(evakey0P, c2QiP, pool2P)
			ringP.MulCoeffsMontgomery(evakey1P, c2QiP, pool3P)
		} else {
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey0Q, c2QiQ, pool2Q)
			ringQ.MulCoeffsMontgomeryAndAddNoModLvl(level, evakey1Q, c2QiQ, pool3Q)
			ringP.MulCoeffsMontgomeryAndAddNoMod(evakey0P, c2QiP, pool2P)
			ringP.MulCoeffsMontgomeryAndAddNoMod(evakey1P, c2QiP, pool3P)
		}

		if reduce&3 == 3 {
			ringQ.ReduceLvl(level, pool2Q, pool2Q)
			ringQ.ReduceLvl(level, pool3Q, pool3Q)
			ringP.Reduce(pool2P, pool2P)
			ringP.Reduce(pool3P, pool3P)
		}

		reduce++
	}

	if (reduce-1)&3 != 3 {
		ringQ.ReduceLvl(level, pool2Q, pool2Q)
		ringQ.ReduceLvl(level, pool3Q, pool3Q)
		ringP.Reduce(pool2P, pool2P)
		ringP.Reduce(pool3P, pool3P)
	}

	ringQ.InvNTTLazy(pool2Q, pool2Q)
	ringQ.InvNTTLazy(pool3Q, pool3Q)
	ringP.InvNTTLazy(pool2P, pool2P)
	ringP.InvNTTLazy(pool3P, pool3P)

	eval.baseconverterQ1P.ModDownSplitPQ(level, pool2Q, pool2P, pool2Q)
	eval.baseconverterQ1P.ModDownSplitPQ(level, pool3Q, pool3P, pool3Q)
}
func (eval *evaluator) getRingQElem(op Operand) *Element {
	switch o := op.(type) {
	case *Ciphertext, *Plaintext:
		return o.El()
	case *PlaintextRingT:
		scaleUp(eval.ringQ, eval.deltaMont, o.value, eval.tmpPt.value)
		return eval.tmpPt.Element
	default:
		panic(fmt.Errorf("invalid operand type for operation: %T", o))
	}
}

// getElemAndCheckBinary unwraps the elements from the operands and checks that the receiver has sufficiently large degree.
func (eval *evaluator) getElemAndCheckBinary(op0, op1, opOut Operand, opOutMinDegree int, ensureRingQ bool) (el0, el1, elOut *Element) {
	if op0 == nil || op1 == nil || opOut == nil {
		panic("operands cannot be nil")
	}

	if op0.Degree()+op1.Degree() == 0 {
		panic("operands cannot be both plaintexts")
	}

	if opOut.Degree() < opOutMinDegree {
		panic("receiver operand degree is too small")
	}

	if ensureRingQ {
		return eval.getRingQElem(op0), eval.getRingQElem(op1), opOut.El() // lifts from Rt to Rq if necessary
	}

	return op0.El(), op1.El(), opOut.El()
}

func (eval *evaluator) getElemAndCheckUnary(op0, opOut Operand, opOutMinDegree int) (el0, elOut *Element) {
	if op0 == nil || opOut == nil {
		panic("operand cannot be nil")
	}

	if op0.Degree() == 0 {
		panic("operand cannot be plaintext")
	}

	if opOut.Degree() < opOutMinDegree {
		panic("receiver operand degree is too small")
	}
	el0, elOut = op0.El(), opOut.El()
	return
}

// evaluateInPlaceBinary applies the provided function in place on el0 and el1 and returns the result in elOut.
func (eval *evaluator) evaluateInPlaceBinary(el0, el1, elOut *Element, evaluate func(*ring.Poly, *ring.Poly, *ring.Poly)) {

	smallest, largest, _ := getSmallestLargest(el0, el1)

	for i := 0; i < smallest.Degree()+1; i++ {
		evaluate(el0.value[i], el1.value[i], elOut.value[i])
	}

	// If the inputs degrees differ, it copies the remaining degree on the receiver.
	if largest != nil && largest != elOut { // checks to avoid unnecessary work.
		for i := smallest.Degree() + 1; i < largest.Degree()+1; i++ {
			elOut.value[i].Copy(largest.value[i])
		}
	}
}

// evaluateInPlaceUnary applies the provided function in place on el0 and returns the result in elOut.
func evaluateInPlaceUnary(el0, elOut *Element, evaluate func(*ring.Poly, *ring.Poly)) {
	for i := range el0.value {
		evaluate(el0.value[i], elOut.value[i])
	}
}

// decomposeAndSplitNTT decomposes the input polynomial into the target CRT basis.
func (eval *evaluator) decomposeAndSplitNTT(level, beta int, c2NTT, c2InvNTT, c2QiQ, c2QiP *ring.Poly) {

	ringQ := eval.ringQ
	ringP := eval.ringP

	eval.decomposer.DecomposeAndSplit(level, beta, c2InvNTT, c2QiQ, c2QiP)

	p0idxst := beta * eval.params.Alpha()
	p0idxed := p0idxst + eval.decomposer.Xalpha()[beta]

	// c2_qi = cx mod qi mod qi
	for x := 0; x < level+1; x++ {

		qi := ringQ.Modulus[x]
		nttPsi := ringQ.NttPsi[x]
		bredParams := ringQ.BredParams[x]
		mredParams := ringQ.MredParams[x]

		if p0idxst <= x && x < p0idxed {
			p0tmp := c2NTT.Coeffs[x]
			p1tmp := c2QiQ.Coeffs[x]
			for j := 0; j < ringQ.N; j++ {
				p1tmp[j] = p0tmp[j]
			}
		} else {
			ring.NTTLazy(c2QiQ.Coeffs[x], c2QiQ.Coeffs[x], ringQ.N, nttPsi, qi, mredParams, bredParams)
		}
	}
	// c2QiP = c2 mod qi mod pj
	ringP.NTTLazy(c2QiP, c2QiP)
}
