package ckks_fv

import (
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
)

// CKKSEncryptor in an interface for encryptors
//
// encrypt with pk : ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk : ciphertext = [-a*sk + m + e, a]
type CKKSEncryptor interface {
	// EncryptNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	// The level of the output ciphertext is plaintext.Level().
	EncryptNew(plaintext *Plaintext) *Ciphertext

	// Encrypt encrypts the input plaintext using the stored key, and returns
	// the result on the receiver ciphertext. The encryption is done by first
	// encrypting zero in QP, dividing by P and then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
	Encrypt(plaintext *Plaintext, ciphertext *Ciphertext)

	// EncryptFastNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	// The level of the output ciphertext is plaintext.Level().
	EncryptFastNew(plaintext *Plaintext) *Ciphertext

	// EncryptFast encrypts the input plaintext using the stored-key, and returns
	// the result on the receiver ciphertext. The encryption is done by first
	// encrypting zero in Q and then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level()).
	EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext)

	// EncryptFromCRPNew encrypts the input plaintext using the stored key and returns
	// the result on a newly created ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), len(CRP.Coeffs)-1).
	EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext

	// EncryptFromCRP encrypts the input plaintext using the stored key and returns
	// the result tge receiver ciphertext. The encryption is done by first encrypting
	// zero in QP, using the provided polynomial as the uniform polynomial, dividing by P and
	// then adding the plaintext.
	// The level of the output ciphertext is min(plaintext.Level(), ciphertext.Level(), len(CRP.Coeffs)-1).
	EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly)
}

// encryptor is a struct used to encrypt Plaintexts. It stores the public-key and/or secret-key.
type ckksEncryptor struct {
	params *Parameters

	ringQ *ring.Ring
	ringP *ring.Ring

	poolQ [3]*ring.Poly
	poolP [3]*ring.Poly

	baseconverter   *ring.FastBasisExtender
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	uniformSampler  *ring.UniformSampler
}

type pkCKKSEncryptor struct {
	ckksEncryptor
	pk *PublicKey
}

type skCKKSEncryptor struct {
	ckksEncryptor
	sk *SecretKey
}

// NewCKKSEncryptorFromPk creates a new Encryptor with the provided public-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewCKKSEncryptorFromPk(params *Parameters, pk *PublicKey) CKKSEncryptor {
	enc := newCKKSEncryptor(params)

	if pk.Value[0].Degree() != params.N() || pk.Value[1].Degree() != params.N() {
		panic("cannot newEncryptor: pk ring degree does not match params ring degree")
	}

	return &pkCKKSEncryptor{enc, pk}
}

// NewCKKSEncryptorFromSk creates a new Encryptor with the provided secret-key.
// This Encryptor can be used to encrypt Plaintexts, using the stored key.
func NewCKKSEncryptorFromSk(params *Parameters, sk *SecretKey) CKKSEncryptor {
	enc := newCKKSEncryptor(params)

	if sk.Value.Degree() != params.N() {
		panic("cannot newEncryptor: sk ring degree does not match params ring degree")
	}

	return &skCKKSEncryptor{enc, sk}
}

func newCKKSEncryptor(params *Parameters) ckksEncryptor {

	var q, p *ring.Ring
	var err error
	if q, err = ring.NewRing(params.N(), params.qi); err != nil {
		panic(err)
	}

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var baseconverter *ring.FastBasisExtender
	var poolP [3]*ring.Poly
	if params.PiCount() != 0 {

		if p, err = ring.NewRing(params.N(), params.pi); err != nil {
			panic(err)
		}

		baseconverter = ring.NewFastBasisExtender(q, p)

		poolP = [3]*ring.Poly{p.NewPoly(), p.NewPoly(), p.NewPoly()}
	}

	return ckksEncryptor{
		params:          params.Copy(),
		ringQ:           q,
		ringP:           p,
		poolQ:           [3]*ring.Poly{q.NewPoly(), q.NewPoly(), q.NewPoly()},
		poolP:           poolP,
		baseconverter:   baseconverter,
		gaussianSampler: ring.NewGaussianSampler(prng),
		ternarySampler:  ring.NewTernarySampler(prng, q, 0.5, false),
		uniformSampler:  ring.NewUniformSampler(prng, q),
	}
}

// EncryptNew encrypts the input Plaintext using the stored key and returns
// the result on a newly created Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkCKKSEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {

	if encryptor.baseconverter == nil {
		panic("Cannot EncryptNew : modulus P is empty -> use instead EncryptFastNew")
	}

	ciphertext := NewCiphertextCKKS(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.encrypt(plaintext, ciphertext, false)

	return ciphertext
}

func (encryptor *pkCKKSEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {

	if encryptor.baseconverter == nil {
		panic("Cannot Encrypt : modulus P is empty -> use instead EncryptFast")
	}

	encryptor.encrypt(plaintext, ciphertext, false)
}

func (encryptor *pkCKKSEncryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	ciphertext := NewCiphertextCKKS(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.encrypt(plaintext, ciphertext, true)

	return ciphertext
}

func (encryptor *pkCKKSEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.encrypt(plaintext, ciphertext, true)
}

func (encryptor *pkCKKSEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

func (encryptor *pkCKKSEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	panic("Cannot encrypt with CRP using an encryptor created with the public-key")
}

// Encrypt encrypts the input Plaintext using the stored key, and returns the result
// on the receiver Ciphertext.
//
// encrypt with pk: ciphertext = [pk[0]*u + m + e_0, pk[1]*u + e_1]
// encrypt with sk: ciphertext = [-a*sk + m + e, a]
func (encryptor *pkCKKSEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, fast bool) {

	lvl := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := encryptor.poolQ[0]
	poolQ1 := encryptor.poolQ[1]
	poolQ2 := encryptor.poolQ[2]
	poolP0 := encryptor.poolP[0]
	poolP1 := encryptor.poolP[1]
	poolP2 := encryptor.poolP[2]

	// We sample a R-WLE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)

	ringQ := encryptor.ringQ

	if fast {

		encryptor.ternarySampler.ReadLvl(lvl, poolQ2)
		ringQ.NTTLvl(lvl, poolQ2, poolQ2)
		ringQ.MFormLvl(lvl, poolQ2, poolQ2)

		// ct0 = u*pk0
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ2, encryptor.pk.Value[0], ciphertext.value[0])
		// ct1 = u*pk1
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ2, encryptor.pk.Value[1], ciphertext.value[1])

		// ct1 = u*pk1 + e1
		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0, ringQ, encryptor.params.sigma, int(6*encryptor.params.sigma))
		ringQ.NTTLvl(lvl, poolQ0, poolQ0)
		ringQ.AddLvl(lvl, ciphertext.value[1], poolQ0, ciphertext.value[1])

		if !plaintext.isNTT {

			// ct0 = u*pk0 + e0
			encryptor.gaussianSampler.ReadLvl(lvl, poolQ0, ringQ, encryptor.params.sigma, int(6*encryptor.params.sigma))
			// ct0 = (u*pk0 + e0)/P + m
			ringQ.AddLvl(lvl, poolQ0, plaintext.value, poolQ0)
			ringQ.NTTLvl(lvl, poolQ0, poolQ0)
			ringQ.AddLvl(lvl, ciphertext.value[0], poolQ0, ciphertext.value[0])

		} else {
			// ct0 = u*pk0 + e0
			encryptor.gaussianSampler.ReadLvl(lvl, poolQ0, ringQ, encryptor.params.sigma, int(6*encryptor.params.sigma))
			ringQ.NTTLvl(lvl, poolQ0, poolQ0)
			ringQ.AddLvl(lvl, ciphertext.value[0], poolQ0, ciphertext.value[0])
			ringQ.AddLvl(lvl, ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}

	} else {

		ringP := encryptor.ringP

		encryptor.ternarySampler.ReadLvl(lvl, poolQ2)

		extendBasisSmallNormAndCenter(ringQ, ringP, poolQ2, poolP2)

		// (#Q + #P) NTT
		ringQ.NTTLvl(lvl, poolQ2, poolQ2)
		ringP.NTT(poolP2, poolP2)

		ringQ.MFormLvl(lvl, poolQ2, poolQ2)
		ringP.MForm(poolP2, poolP2)

		pk0P := new(ring.Poly)
		pk1P := new(ring.Poly)
		pk0P.Coeffs = encryptor.pk.Value[0].Coeffs[len(ringQ.Modulus):]
		pk1P.Coeffs = encryptor.pk.Value[1].Coeffs[len(ringQ.Modulus):]

		// ct0 = u*pk0
		// ct1 = u*pk1
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ2, encryptor.pk.Value[0], poolQ0)
		ringQ.MulCoeffsMontgomeryLvl(lvl, poolQ2, encryptor.pk.Value[1], poolQ1)
		ringP.MulCoeffsMontgomery(poolP2, pk0P, poolP0)
		ringP.MulCoeffsMontgomery(poolP2, pk1P, poolP1)

		// 2*(#Q + #P) NTT
		ringQ.InvNTTLvl(lvl, poolQ0, poolQ0)
		ringQ.InvNTTLvl(lvl, poolQ1, poolQ1)
		ringP.InvNTT(poolP0, poolP0)
		ringP.InvNTT(poolP1, poolP1)

		// ct0 = u*pk0 + e0
		encryptor.gaussianSampler.ReadLvl(lvl, poolQ2, ringQ, encryptor.params.sigma, int(6*encryptor.params.sigma))
		extendBasisSmallNormAndCenter(ringQ, ringP, poolQ2, poolP2)
		ringQ.AddLvl(lvl, poolQ0, poolQ2, poolQ0)
		ringP.Add(poolP0, poolP2, poolP0)

		// ct1 = u*pk1 + e1
		encryptor.gaussianSampler.ReadLvl(lvl, poolQ2, ringQ, encryptor.params.sigma, int(6*encryptor.params.sigma))
		extendBasisSmallNormAndCenter(ringQ, ringP, poolQ2, poolP2)
		ringQ.AddLvl(lvl, poolQ1, poolQ2, poolQ1)
		ringP.Add(poolP1, poolP2, poolP1)

		// ct0 = (u*pk0 + e0)/P
		encryptor.baseconverter.ModDownSplitPQ(lvl, poolQ0, poolP0, ciphertext.value[0])

		// ct1 = (u*pk1 + e1)/P
		encryptor.baseconverter.ModDownSplitPQ(lvl, poolQ1, poolP1, ciphertext.value[1])

		if !plaintext.isNTT {
			ringQ.AddLvl(lvl, ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}

		// 2*#Q NTT
		ringQ.NTTLvl(lvl, ciphertext.value[0], ciphertext.value[0])
		ringQ.NTTLvl(lvl, ciphertext.value[1], ciphertext.value[1])

		if plaintext.isNTT {
			// ct0 = (u*pk0 + e0)/P + m
			ringQ.AddLvl(lvl, ciphertext.value[0], plaintext.value, ciphertext.value[0])
		}
	}

	ciphertext.value[0].Coeffs = ciphertext.value[0].Coeffs[:lvl+1]
	ciphertext.value[1].Coeffs = ciphertext.value[1].Coeffs[:lvl+1]

	ciphertext.isNTT = true
}

func (encryptor *skCKKSEncryptor) EncryptNew(plaintext *Plaintext) *Ciphertext {
	ciphertext := NewCiphertextCKKS(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.Encrypt(plaintext, ciphertext)
	return ciphertext
}

func (encryptor *skCKKSEncryptor) Encrypt(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.encryptSample(plaintext, ciphertext)
}

func (encryptor *skCKKSEncryptor) EncryptFastNew(plaintext *Plaintext) *Ciphertext {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFastNew() -> use instead EncryptNew()")
}

func (encryptor *skCKKSEncryptor) EncryptFast(plaintext *Plaintext, ciphertext *Ciphertext) {
	panic("Cannot Encrypt : SkEncryptor doesn't support EncryptFast() -> use instead Encrypt()")
}

func (encryptor *skCKKSEncryptor) EncryptFromCRPNew(plaintext *Plaintext, crp *ring.Poly) *Ciphertext {
	ciphertext := NewCiphertextCKKS(encryptor.params, 1, plaintext.Level(), plaintext.Scale())
	encryptor.EncryptFromCRP(plaintext, ciphertext, crp)
	return ciphertext
}

func (encryptor *skCKKSEncryptor) EncryptFromCRP(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {
	encryptor.ringQ.Copy(crp, ciphertext.value[1])
	ciphertext.value[0].Coeffs = ciphertext.value[0].Coeffs[:len(crp.Coeffs)]
	ciphertext.value[1].Coeffs = ciphertext.value[1].Coeffs[:len(crp.Coeffs)]
	encryptor.encrypt(plaintext, ciphertext, ciphertext.value[1])
}

func (encryptor *skCKKSEncryptor) encryptSample(plaintext *Plaintext, ciphertext *Ciphertext) {
	encryptor.uniformSampler.Readlvl(utils.MinInt(plaintext.Level(), ciphertext.Level()), ciphertext.value[1])
	encryptor.encrypt(plaintext, ciphertext, ciphertext.value[1])
}

func (encryptor *skCKKSEncryptor) encrypt(plaintext *Plaintext, ciphertext *Ciphertext, crp *ring.Poly) {

	ringQ := encryptor.ringQ

	lvl := utils.MinInt(plaintext.Level(), ciphertext.Level())

	poolQ0 := encryptor.poolQ[0]

	ringQ.MulCoeffsMontgomeryLvl(lvl, ciphertext.value[1], encryptor.sk.Value, ciphertext.value[0])
	ringQ.NegLvl(lvl, ciphertext.value[0], ciphertext.value[0])

	if plaintext.isNTT {
		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0, ringQ, encryptor.params.sigma, int(6*encryptor.params.sigma))
		ringQ.NTTLvl(lvl, poolQ0, poolQ0)
		ringQ.AddLvl(lvl, ciphertext.value[0], poolQ0, ciphertext.value[0])
		ringQ.AddLvl(lvl, ciphertext.value[0], plaintext.value, ciphertext.value[0])
	} else {
		encryptor.gaussianSampler.ReadLvl(lvl, poolQ0, ringQ, encryptor.params.sigma, int(6*encryptor.params.sigma))
		ringQ.AddLvl(lvl, poolQ0, plaintext.value, poolQ0)
		ringQ.NTTLvl(lvl, poolQ0, poolQ0)
		ringQ.AddLvl(lvl, ciphertext.value[0], poolQ0, ciphertext.value[0])
	}

	ciphertext.value[0].Coeffs = ciphertext.value[0].Coeffs[:lvl+1]
	ciphertext.value[1].Coeffs = ciphertext.value[1].Coeffs[:lvl+1]

	ciphertext.isNTT = true
}

func extendBasisSmallNormAndCenter(ringQ, ringP *ring.Ring, polQ, polP *ring.Poly) {
	var coeff, Q, QHalf, sign uint64
	Q = ringQ.Modulus[0]
	QHalf = Q >> 1

	for j := 0; j < ringQ.N; j++ {

		coeff = polQ.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i, pi := range ringP.Modulus {
			polP.Coeffs[i][j] = (coeff * sign) | (pi-coeff)*(sign^1)
		}

	}
}
