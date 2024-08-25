#! /usr/bin/env sage

log2 = lambda x: log(1.0*x)/log(2.0)

time_halfBTP = 20

instances = {

    'Kreyvium12': (12, 46, 1, 2, [ 440.092 ]), # phi(m)=14112 // L:250 k=107
    'Kreyvium13': (13, 125, 1, 2, [ 515.639 ]), # phi(m)=15004 // L:270 k=113

    'FiLIP1216': (1216, 64, 1, 2, [ 1192.232 ]), # phi(m)=5292 // L:100 k=113
    'FiLIP1280': (1280, 64, 1, 2, [ 1889.830 ]), # phi(m)=5760 // L:110 k=119

    'Rasta5': (525, 525, 1, 2, [ 255.304 ]),  #phi(m)=10752 // L:200 k=110
    'Rasta5_packed': (525, 525, 1, 2, [ 33.837 ]), # phi(m)=15004 // L:270 k=113
    'Rasta6': (351, 351, 1, 2, [ 182.912 ]), # phi(m)=10752 // L:200 k=110
    'Rasta6_packed': (351, 351, 1, 2, [ 45.986 ]), # phi(m)=15004 // L:300 k=100

    'Dasta5': (525, 525, 1, 2, [ 254.351 ]),  #phi(m)=10752 // L:200 k=110
    'Dasta5_packed': (525, 525, 1, 2, [ 35.130 ]), # phi(m)=15004 // L:270 k=113
    'Dasta6': (351, 351, 1, 2, [ 184.411 ]), # phi(m)=10752 // L:200 k=110
    'Dasta6_packed': (351, 351, 1, 2, [ 47.086 ]), # phi(m)=15004 // L:300 k=100

    # (block size, output size, slots num, plaintext modulu, dec time, )
    'Masta5S': (64, 64, 65536, 0x1fc0001, [2400.0 ]),

    # 'Masta6S': (32, 32, 65536, 0x1fc0001, [933.0 ]),
    'Masta6M': (32, 32, 65536, 0x3e00001, [938.0 ]),
    # 'Masta7S': (16, 16, 65536, 0x1fc0001, [349.0 ]),
    'Masta7M': (16, 16, 65536, 0x3e00001, [350.0 ]),

    # 'Pasta4S': (64, 32, 65536, 0x1fc0001, [1483.0 ]),
    'Pasta4M': (64, 32, 65536, 0x3e00001, [1484.0 ]),
    # 'Pasta5S': (32, 16, 65536, 0x1fc0001, [495.8 ]),
    'Pasta5M': (32, 16, 65536, 0x3e00001, [496.0 ]),

    # 'HeraS': (16, 16, 65536, 0x1fc0001, 137),
    'Hera': (16, 16, 65536, 0xffa0001, [138.0 ]),

    'Rubato5S': (16, 12, 65536, 0x3e00001, [86.0 ]),
    'Rubato3M': (36, 32, 65536, 0x1fc0001, [114.0 ]),
    'Rubato2L': (64, 60, 65536, 0x1fc0001, [148.0 ]),

    'AES': (128, 128, 32768, 2, [2300.0 ]),
    'GGIFT': (64, 64, 32768, 2, [1100.0 ]),
    'FiLIP-1216': (128, 1, 1, 2, [24.37 ]),
    'FiLIP-1280': (128, 1, 1, 2, [26.59 ]),
}

for inst in instances.keys():

    print ("## "+inst)
    blockSize, outSize, slotsNum, t, decTimeAll = instances[inst]
    decTime = 0
    for item in decTimeAll :
        decTime += item
    decTime /= len(decTimeAll)

    if t == 2:
        print ("Parameters: blocksize={} output size={} slots num={} plain modulus bits={:.2f} Decrypt Time={:.2f}s" \
           .format(blockSize, outSize, slotsNum, log2(t), decTime) )
        totalTime = decTime + 20
        totalData = outSize * slotsNum * log2(t)
        throughput = totalData / ( totalTime * 8 * 1024 )
        chi2_time = totalTime + 20 # chi2 time 20s
        resnet_time = totalTime + 1091 # resnet infer time 1091s

        print ( "Troughput           : {:.9f}".format( throughput ) )
        print ( "Chi^2 E2E time           : {:.3f}".format( chi2_time ) )
        print ( "NN infer E2E time           : {:.3f}\n".format( resnet_time ) )

    else:
        print ("Parameters: blocksize={} output size={} slots num={} plain modulus bits={:.2f} Decrypt Time={:.2f}s" \
            .format(blockSize, outSize, slotsNum, log2(t), decTime) )
        
        FVtoCKKStime = (time_halfBTP * outSize * 2) / 64
        totalTime = decTime + FVtoCKKStime
        latency = decTime + time_halfBTP

        totalData = outSize * slotsNum * log2(t)
        throughput = totalData / ( totalTime * 8 * 1024 )
        
        chi2_time = latency + 20 # chi2 time 20s
        resnet_time = latency + 1064 # resnet infer time 1063s

        print ( "Troughput           : {:.3f}".format( throughput ) )
        print ( "Chi^2 E2E time           : {:.3f}".format( chi2_time ) )
        print ( "NN infer E2E time           : {:.3f}\n".format( resnet_time ) )


print( '\n\n\n\n\n') 
print("程序即将终止")
exit()

Client_instances = {
    # (block size, output size, plaintext modulu, enc time )

    'LowMCv3_14_196_63': (196, 196, 1, [63997, 63897] ),
    'LowMCv3_14_256_63': (256, 256, 1, [83576, 83778, 83351, 84539, 83575] ),

    'Kreyvium12': (288, 46, 1, [19, 18, 19, 20]),
    'Kreyvium13': (288, 125, 1, [19, 20, 20, 19]),

    'Rasta_5': (525, 525, 1, [67118, 66844, 67059, 66973, 66853, 66864]),
    'Rasta_6': (351, 351, 1, [16757, 16835, 16802, 16765, 16813, 16820]),

    'Dasta_5': (525, 525, 1, [637, 632, 632, 636, 626, 618]),
    'Dasta_6': (351, 351, 1, [288, 288, 286, 285, 287, 290]),

    'FiLIP1216': (1216, 46, 1, [17357, 17403, 17319, 17375, 17369]),
    'FiLIP1280': (1280, 46, 1, [4721, 4784, 4722, 4751, 4742]),

    'Masta4S': (128, 128, 65537, [895, 841, 877, 823, 826]),
    'Masta4M': (128, 128, 8088322049, [812, 847, 812, 815, 933]),
    'Masta4L': (128, 128, 1096486890805657601, [819, 846, 896, 910, 846]),

    'Masta5S': (64, 64, 65537, [259, 279, 268, 275, 277]),
    'Masta5M': (64, 64, 8088322049, [265, 266, 294, 262, 256]),
    'Masta5L': (64, 64, 1096486890805657601, [275, 260, 255, 262, 267]),

    'Masta6S': (32, 32, 0x1fc0001, [83, 84, 82, 83, 83, 80]), 
    'Masta6M': (32, 32, 0x3e00001, [82, 88, 82, 81, 83] ), 

    'Masta7S': (16, 16, 0x1fc0001, [33, 33, 33, 32, 33]), 
    'Masta7M': (16, 16, 0x3e00001, [33, 34, 33, 32, 33] ), 

    'Pasta3S': (128, 64, 65537, [1801, 1761, 1768, 1852, 1843]), 
    'Pasta3M': (128, 64, 8088322049, [1665, 1723, 1685, 1719, 1649]),

    'Pasta4S': (64, 32, 0x1fc0001, [174, 174, 181, 174, 175]),
    'Pasta4M': (64, 32, 0x3e00001, [175, 175, 173, 174, 174]),

    'Pasta5S': (32, 16, 0x1fc0001, [60, 56, 55, 65, 56, 56]),
    'Pasta5M': (32, 16, 0x3e00001, [58, 56, 63, 75, 56, 56]),

    'Pasta2_3': (128, 64, 0x1fc0001, [797, 798, 798, 797, 807]),
    'Pasta2_4': (64, 32, 0x3e00001, [71, 71, 72, 71, 71]),

    'Hera5': (16, 16, 0xffa0001, [11, 10, 11, 10, 11]),

    'Rubato5S': (16, 12, 0x3e00001, [4531, 4527.27, 4543.58, 4536.75, 4528.58]),
    'Rubato3M': (36, 32, 0x1fc0001, [5054.39, 5054.65, 5055.35, 5055.2, 5054.06]),
    'Rubato2L': (64, 60, 0x1fc0001, [5471.57, 5477.69, 5469.04, 5474.36, 5472]),

}

for inst in Client_instances.keys():

    blockSize, outSize, t, decCycleAll = Client_instances[inst]

    decCycle = 0
    for item in decCycleAll:
        decCycle += item
    decCycle /= len(decCycleAll)
    decCycle *= 2700
    print ("## "+inst)
    print ("Parameters: blocksize={} output size={} plain modulus bits={:.2f} Decrypt Cycles={}" \
    .format(blockSize, outSize, log2(t), decCycle) )
    
    if t == 1:
        totalBytes = ( log2(2) * outSize ) / 8
    else:
        totalBytes = ( log2(t) * outSize ) / 8

    CperB = decCycle / totalBytes

    print ( "Cycle per Byte is           : {:.3f}".format( CperB ) )
    print()


AES_instance = {
    'AES-NI': (128, 128, 1, [38.55, 38.55, 38.55, 38.67]),
    'AES-SW': (128, 128, 1, [233.12]),
    'GGIFT': (64, 64, 1, [26])
}

for inst in AES_instance.keys():

    blockSize, outSize, t, decCycleAll = AES_instance[inst]

    decCycle = 0
    for item in decCycleAll:
        decCycle += item
    decCycle /= len(decCycleAll)
    print ("\n\n\n## "+inst)
    print ("Parameters: blocksize={} output size={} plain modulus bits={:.2f} Decrypt Cycles={}" \
    .format(blockSize, outSize, log2(t), decCycle) )
    
    if t == 1:
        totalBytes = ( log2(2) * outSize ) / 8
    else:
        totalBytes = ( log2(t) * outSize ) / 8

    CperB = decCycle / totalBytes

    print ( "Cycle per Byte is           : {:.3f}".format( CperB ) )
    print()
