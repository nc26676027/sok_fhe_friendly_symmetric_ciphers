#! /usr/bin/env sage

log2 = lambda x: log(1.0*x)/log(2.0)

time_halfBTP = 20000
compare_length = 32
MUX_time = 26

instances = {
    # (block size, output size, slots num, plaintext modulu, dec time, )
    'TFHE LowMC': (128, 128, 1, 2, [ 5368000 ], 64),
    # 'helib LowMC': (128, 128, 1, 2, [ 1037285 ], 1),
    
    'TFHE Kreyvium12': (288, 64, 1, 2, [ 248633 ], 64),
    'TFHE Kreyvium13': (288, 125, 1, 2, [ 266377 ], 64),

    'TFHE FiLIP1216': (1216, 64, 1, 2, [ 1324213 ], 64),
    'TFHE FiLIP1280': (1280, 64, 1, 2, [ 1358175 ], 64),

    'TFHE Rasta5': (525, 525, 1, 2, [ 4921574 ], 64),
    'TFHE Rasta6': (351, 351, 1, 2, [ 2858397 ], 64),

    'TFHE Dasta5': (525, 525, 1, 2, [ 5031231 ], 64),
    'TFHE Dasta6': (351, 351, 1, 2, [ 2933646 ], 64),

    'TFHE AES': (128, 128, 1, 2, [ 252000 ], 1),

    ##########################  Fp ciphers

    # 'Masta4S_helib': (128, 128, 32768, 65537, [ 28214 ], 1),
    # 'Masta4M_helib': (128, 128, 65536, 8088322049, [87698.0 ], 1),
    'Masta4M_seal': (128, 128, 65536, 8088322049, [ 219811 ], 1),

    # 'Masta5S_helib': (64, 64, 32768, 65537, [ 70089  ], 1),
    # 'Masta5M_helib': (64, 64, 65536, 8088322049, [ 69508 ], 1),
    'Masta5M_seal': (64, 64, 65536, 8088322049, [ 164377 ], 1),

    # 'Masta6S': (32, 32, 65536, 0x1fc0001, [933.0 ]),
    # 'Masta6M': (32, 32, 65536, 0x3e00001, [938.0 ]),
    # 'Masta7S': (16, 16, 65536, 0x1fc0001, [349.0 ]),
    # 'Masta7M': (16, 16, 65536, 0x3e00001, [350.0 ]),

    # 'Pasta3S': (64, 32, 65536, 0x1fc0001, [1483.0 ]),
    'Pasta3M': (256, 128, 65536, 8088322049, [ 172891.0 ], 1),
    # 'Pasta4S': (32, 16, 65536, 0x1fc0001, [495.8 ]),
    'Pasta4M': (64, 32, 65536, 8088322049, [ 97043 ], 1),

    'Pasta2_3M': (256, 128, 65536, 8088322049, [180590.0 ], 1),
    'Pasta2_4M': (64, 32, 65536, 8088322049, [92386.0 ], 1),

    # 'HeraS': (16, 16, 65536, 0x1fc0001, 137),
    'Hera_seal': (16, 16, 65536, 0x1e21a0001, [ 76742.0 ], 1),

    'Rubato5S': (16, 12, 65536, 0x3e00001, [79446 ], 1), # 5 rounds masta7 
    'Rubato3M': (36, 32, 65536, 0x1fc0001, [77647 ], 1), # 3 rounds masta6
    'Rubato2L': (64, 60, 65536, 0x1fc0001, [78722 ], 1), # 2 rounds 

    
    # 'GGIFT': (64, 64, 32768, 1, [1100.0 ]),
}

for inst in instances.keys():

    print ("## "+inst)
    blockSize, outSize, slotsNum, t, decTimeAll, threads = instances[inst]
    decTime = 0
    for item in decTimeAll :
        decTime += item
    decTime /= len(decTimeAll)
    decTime /= float(threads)

    if t == 2:
        print ("Parameters: blocksize={} output size={} slots num={} Decrypt Time={:.2f}ms" \
           .format(blockSize, outSize, slotsNum, decTime) )
        
        tfhe_compare_time = 7 * compare_length * MUX_time # 7n MUX, where MUX 26ms
        totalTime = decTime + tfhe_compare_time
        print ( "TFHE two {}-bits number compare time           : {:.2f}ms\n".format( compare_length, totalTime ) )

    else:
        print ("Parameters: blocksize={} output size={} slots num={} plain modulus bits={:.2f} Decrypt Time={:.2f}ms" \
            .format(blockSize, outSize, slotsNum, log2(t), decTime) )
        
        FVtoCKKStime = time_halfBTP
        Max_time = 6630 # Min/Max evaluating time in CKKS
        totalTime = decTime + FVtoCKKStime + Max_time

        print ( "BFV->CKKS two {}-bits number compare time           : {:.2f}ms\n".format( compare_length, totalTime ) )

