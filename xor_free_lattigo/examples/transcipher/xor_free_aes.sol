Residual parameters: logN=15, logSlots=14, H=192, sigma={3.200000 19.200000}, logQP=840.000013, levels=15, scale=2^41
Bootstrapping parameters: logN=15, logSlots=14, H(192; 32), sigma={3.200000 19.200000}, logQP=840.000013, levels=16, scale=2^41

Generating bootstrapping evaluation keys...
Done
round iterator : 1
Chain index before sbox: 14, scale: 39.999997
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 1.0013,  1.0000,  1.0068,  1.0025,  0.9999,  1.0014,  1.0013 ]
 ... 1.0000,  1.0000,  1.0012,  0.9995,  1.0002,  0.9997,  1.0001 ]

round iterator : 2
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 1.0000,  1.0000,  0.9967,  1.0000,  0.9999,  1.0000,  0.9991 ]
 ... 0.9978,  1.0000,  0.9994,  0.9998,  0.9997,  1.0004,  1.0015 ]

round iterator : 3
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.2013,  0.3665,  0.4448,  0.3705,  0.3470,  0.4197,  0.2428 ]
 ... 0.3592,  0.3541,  0.1975,  0.3515,  0.3678,  0.3461,  0.3742 ]

round iterator : 4
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.4782,  0.6551,  0.6489,  0.6564,  0.4744,  0.7813,  0.5293 ]
 ... 0.6373,  0.6368,  0.2563,  0.5603,  0.7283,  0.5825,  0.6884 ]

round iterator : 5
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.6958,  0.1487,  0.6740,  0.1239,  0.1497,  0.2869,  0.0017 ]
 ... 0.0207,  0.0168,  0.2857,  0.0870,  0.0112,  0.1130,  0.1232 ]

round iterator : 6
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.9412,  0.9521,  0.0709,  0.1466,  0.9953,  0.9194,  0.1245 ]
 ... 0.3787,  0.9998,  0.6779,  0.3054,  0.1139,  0.7913,  0.9923 ]

round iterator : 7
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.5469,  0.5389,  0.0043,  0.6988,  0.3189,  0.7581,  0.5902 ]
 ... 0.5549,  0.9946,  0.0030,  0.0206,  0.4623,  0.7679,  0.8308 ]

round iterator : 8
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.7878,  0.1309,  0.1488,  0.9923,  0.9853,  0.9923,  0.7181 ]
 ... 0.9888,  0.4840,  0.5225,  0.0002,  0.5140,  0.0571,  0.1667 ]

round iterator : 9
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.1119,  0.2615,  0.8441,  0.1265,  0.5329,  0.0027,  0.2224 ]
 ... 0.9171,  0.0088,  0.0331,  0.9004,  0.9998,  0.5916,  0.5450 ]

round iterator : last round
Chain index before sbox: 5, scale: 41.000000
Chain Index: 5, Scale: 41.00
BTS precise: 
[ 0.5273,  0.9615,  0.0501,  0.5965,  0.3090,  0.9240,  0.0071 ]
 ... 0.1191,  0.0020,  0.0465,  0.0004,  0.8528,  0.9889,  0.9070 ]

代码执行时间：334 秒 :: 504 毫秒
panic: finished

goroutine 1 [running]:
github.com/tuneinsight/lattigo/v6/ckks_cipher.(*AESCtr).HEDecrypt(0xc0002dc000, {0xc000020f80?, 0x10?, 0x10?}, 0xc00001c460?)
	/home/niuchao/sok_fhe_friendly_symmetric_ciphers/xor_free_lattigo/ckks_cipher/aes_ctr.go:133 +0x725
main.main()
	/home/niuchao/sok_fhe_friendly_symmetric_ciphers/xor_free_lattigo/examples/transcipher/main.go:160 +0x1705
exit status 2
