Residual parameters: logN=14, logSlots=13, H=192, sigma={3.200000 19.200000}, logQP=720.000008, levels=9, scale=2^40
Bootstrapping parameters: logN=14, logSlots=13, H(192; 32), sigma={3.200000 19.200000}, logQP=1541.000010, levels=25, scale=2^40

Generating bootstrapping evaluation keys...
Done
round iterator : 1
Chain index before sbox: 3, scale: 39.999996
Chain Index: 0, Scale: 40.00
after Sbox
[ 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]
 ... 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]

Chain Index: 0, Scale: 40.00
after Sbox
[ 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]
 ... 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]

Chain Index: 0, Scale: 40.00
after Sbox
[ 0.0000,  0.0000,  0.0000,  0.0000,  -0.0000,  -0.0000,  -0.0000 ]
 ... 0.0000,  0.0000,  -0.0000,  0.0000,  -0.0000,  0.0000,  0.0000 ]

Chain Index: 0, Scale: 40.00
after Sbox
[ -0.0000,  -0.0000,  0.0000,  -0.0000,  0.0000,  -0.0000,  0.0000 ]
 ... -0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  -0.0000,  0.0000 ]

Chain Index: 0, Scale: 40.00
after Sbox
[ 0.0000,  -0.0000,  0.0000,  -0.0000,  0.0000,  0.0000,  0.0000 ]
 ... 0.0000,  0.0000,  0.0000,  -0.0000,  0.0000,  -0.0000,  -0.0000 ]

Chain Index: 0, Scale: 40.00
after Sbox
[ 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]
 ... 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]

Chain Index: 0, Scale: 40.00
after Sbox
[ 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]
 ... 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000 ]

Chain Index: 0, Scale: 40.00
after Sbox
[ -0.0000,  0.0000,  0.0000,  0.0000,  -0.0000,  0.0000,  -0.0000 ]
 ... -0.0000,  -0.0000,  0.0000,  -0.0000,  0.0000,  -0.0000,  -0.0000 ]

panic: cannot RealToComplex: opOut ring degree must be twice ctIn ring degree

goroutine 228 [running]:
github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping.Evaluator.RealToComplexNew({{{{0xe, {...}, {...}, {...}, {...}, 0xc0000e22d0, 0xc0000e2050, 0x0, {...}, 0x1}}, ...}, ...}, ...)
	/home/niuchao/ckks-rns-seal/lattigo/circuits/ckks/bootstrapping/evaluator.go:861 +0x105
github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping.Evaluator.EvaluateConjugateInvariant({{{{0xe, {...}, {...}, {...}, {...}, 0xc0000e22d0, 0xc0000e2050, 0x0, {...}, 0x1}}, ...}, ...}, ...)
	/home/niuchao/ckks-rns-seal/lattigo/circuits/ckks/bootstrapping/evaluator.go:498 +0x49
github.com/tuneinsight/lattigo/v6/ckks_cipher.BootstrapConjInv(...)
	/home/niuchao/ckks-rns-seal/lattigo/ckks_cipher/rtb_cipher.go:68
github.com/tuneinsight/lattigo/v6/ckks_cipher.(*AESCtr).RoundFunction.func2(0x2)
	/home/niuchao/ckks-rns-seal/lattigo/ckks_cipher/aes_ctr.go:299 +0xca
created by github.com/tuneinsight/lattigo/v6/ckks_cipher.(*AESCtr).RoundFunction in goroutine 1
	/home/niuchao/ckks-rns-seal/lattigo/ckks_cipher/aes_ctr.go:298 +0x426
panic: cannot RealToComplex: opOut ring degree must be twice ctIn ring degree

goroutine 226 [running]:
github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping.Evaluator.RealToComplexNew({{{{0xe, {...}, {...}, {...}, {...}, 0xc0000e22d0, 0xc0000e2050, 0x0, {...}, 0x1}}, ...}, ...}, ...)
	/home/niuchao/ckks-rns-seal/lattigo/circuits/ckks/bootstrapping/evaluator.go:861 +0x105
github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping.Evaluator.EvaluateConjugateInvariant({{{{0xe, {...}, {...}, {...}, {...}, 0xc0000e22d0, 0xc0000e2050, 0x0, {...}, 0x1}}, ...}, ...}, ...)
	/home/niuchao/ckks-rns-seal/lattigo/circuits/ckks/bootstrapping/evaluator.go:498 +0x49
github.com/tuneinsight/lattigo/v6/ckks_cipher.BootstrapConjInv(...)
	/home/niuchao/ckks-rns-seal/lattigo/ckks_cipher/rtb_cipher.go:68
github.com/tuneinsight/lattigo/v6/ckks_cipher.(*AESCtr).RoundFunction.func2(0x0)
	/home/niuchao/ckks-rns-seal/lattigo/ckks_cipher/aes_ctr.go:299 +0xca
created by github.com/tuneinsight/lattigo/v6/ckks_cipher.(*AESCtr).RoundFunction in goroutine 1
	/home/niuchao/ckks-rns-seal/lattigo/ckks_cipher/aes_ctr.go:298 +0x426
exit status 2
