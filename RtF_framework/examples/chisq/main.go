package main

import (
	"runtime"

	"github.com/ldsec/lattigo/v2/chisqtest"
)

func main() {
	runtime.GOMAXPROCS(1)
	// ./chi2 --SNPdir "../data" --SNPfilename "random_sample" --pvalue "pvalue.txt" --runtime "result.txt" --samplesize="200" --snps="16384"
	
	SNPDir := "../../chisqtest/data"
	SNPFileName := "random_sample"
	pValue := "pvalue.txt"
	Runtime := "result.txt"
	SampleSize := "200"
	SNPs := "16384"
	chisqtest.RunChi2(SNPDir, SNPFileName, pValue, Runtime, SampleSize, SNPs);

}