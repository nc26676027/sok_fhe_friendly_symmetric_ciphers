# Alongside code for paper : "SoK: FHE-Friendly Symmetric Ciphers and Transciphering"

Please see the [Wiki](https://github.com/nc26676027/sok_fhe_friendly_symmetric_ciphers/wiki) for more information.

## This source code consists of three parts. For each part, you should enter the respective subdirectory and run the code.

# 1. RtF framework part

## New Package
We implement the hybrid framework in [./RtF_framework/ckks_fv/](./RtF_framework/ckks_fv/), which contains the following functionalities.
- Evaluation of the HERA/Rubato cipher in the FV scheme
- Evaluation of the Masta/Pasta cipher in the FV scheme

We implement the ResNet-20 cnn_infer in [./RtF_framework/cnn_infer/](./RtF_framework/cnn_infer/), which contains the following functionalities.
- ResNet20 inference scheme (located at (./cnn_infer))
- Masta/Pasta cipher in the FV scheme subquent the cnn_infer
- HERA/Rubato cipher in the FV scheme subquent the cnn_infer

We implement the GAWS chisqtest in [./RtF_framework/chisqtest/](./RtF_framework/chisqtest/), which contains the following functionalities.
- chisq test scheme (located at (./chisq.go))
- Masta/Pasta cipher in the FV scheme subquent the chisqtest
- HERA/Rubato cipher in the FV scheme subquent the chisqtest

We implement the MiniMax compare in [./RtF_framework/minimax_comp/](./RtF_framework/minimax_comp/), which contains the following functionalities.
- compare scheme (located at (./minimax_comp.go))
- Masta/Pasta cipher in the FV scheme subquent the minimax_comp
- HERA/Rubato cipher in the FV scheme subquent the minimax_comp

An example of Running ResNet20 inference in the RtF framework is given in [examples/resnet](./RtF_framework/examples/resnet/main.go).

An example of Running compare in the RtF framework is given in [examples/compare](./RtF_framework/examples/compare/main.go).

An example of Running GWAS chisqtest in the RtF framework is given in [examples/chisq](./RtF_framework/examples/chisq/main.go).

## Run RtF experiment

cd into corresponding dir
```PowerShell
cd ./RtF_framework/examples/chisq/main.go
go run main.go
```

# 2. BtR framework part

We implement the hybrid framework in [./BtR_framework/ckks_cipher/](./BtR_framework/ckks_cipher/), which contains the following functionalities.
- Evaluation of the AES-CTR in the CKKS scheme


We implement the ResNet-20 cnn_infer in [./BtR_framework/cnn_infer/](./BtR_framework/cnn_infer/), which contains the following functionalities.
- ResNet20 inference scheme (located at (./cnn_infer))

We implement the GAWS chisqtest in [./BtR_framework/chisqtest/](./BtR_framework/chisqtest/), which contains the following functionalities.
- chisq test scheme (located at (./chisq.go))

An example of Running ResNet20 inference in the BtR framework is given in [examples/resnet](./BtR_framework/examples/resnet/main.go).

An example of Running GWAS chisqtest in the BtR framework is given in [BtR_framework/chisqtest/](./BtR_framework/chisqtest/chisq.go).

## Run BtR experiment

cd into corresponding dir
```PowerShell
cd ./BtR_framework/examples/btr_cipher/main.go
go run main.go
```

# 3. CGGI-based transciphering

This folder contains all the experiments we reevaluated using the open-sourced code obtained from corresponding papers.
