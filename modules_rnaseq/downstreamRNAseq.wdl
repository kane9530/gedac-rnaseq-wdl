version 1.0

task basicRNAseqAnalysis {
	meta{
        description: "Runs the EdgeR-limma-voom pipeline on the unfiltered count matrix for downstream analysis."
        author: "Kane Toh"
        email: "kanetoh@nus.edu.sg"
    }
	input {
        String docker
		Int dockerMemoryGB = 10
		Int numberCpuThreads = 4
		Int numberMaxRetries = 3
		Int awsBatchRetryAttempts = 3
        File countMatrix
		Array[String] samplesName
		Array[String] samplesType
		Array[String] library
		Array[String] lane
		Array[Int] timepoint
		String species
	}

	runtime {
        docker: "~{docker}" + "downstreamrnaseq"
		memory: "~{dockerMemoryGB}G"
		cpu: numberCpuThreads
		maxRetries: numberMaxRetries
		awsBatchRetryAttempts: awsBatchRetryAttempts
    }

	command <<<
		set -e
		outputDir="results"
		mkdir -p $outputDir
		Rscript --vanilla /scripts/main.R ~{species} ~{sep="," samplesName} ~{sep="," samplesType} ~{sep="," lane} ~{sep="," library} ~{sep="," timepoint} ~{countMatrix} $outputDir
		>>>

	output {
		Array[File] resultsDir = glob("results/*")
	}

}