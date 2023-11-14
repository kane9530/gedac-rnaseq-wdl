version 1.0

import "dataStructures.wdl"

workflow trimQC {
    meta {
        description: "Prepare fastq files for processing, and report QC metrics. Especifically, remove adapters, process UMI (if necessary), and split reads to parallelize mapping"
        author: "Juan Felipe Ortiz"
        email: "jf.ortiz@nus.edu.sg"
    }
    parameter_meta {
        designMatrix: {
            description: "DesignMatrix object containing list of read pairs in an experimental group. Each element of the list is a pair of fastq files correspoing to the experimental group",
            value: "The DesignMatrix type contains three fields: sampleName: name of sample sampleType: can be either 'case' or 'control'fastqs: list of read pairs. Each element of the list is a sample, and the read pair is a list with two elements, containing paths two paired-end fastq files"
        }
        umiLoc: {
            description: "Location of UMI in the read",
            value: "index1: the first index is used as UMI. In paired-end reads, this UMI is used for both reads. index2: similar to index 1. read1: the head of read 1 is used as UMI. In paired-end reads, this UMI is used for both reads. per_index: index1_index2 is used as UMI for read1 and read2. per_read: head of read1 is UMI of read 1, and head of read2 is UMI of read2. empty string (''): do not process UMI"
        }
        umiLen: {
            description: "Length of UMI to be processed. This is relevant when umiLoc is either read1, read2, or per_read",
            value: "integer"
        }
        chopReads: {
            description: "Specifies if reads should be split for parallel mapping",
            value: "boolean"
        }
        numberCpuThreads: {
            description: "Number of threads to use on the fastp task",
            value: "Int"
        }
        dockerMemoryGB: {
            description: "Memory to request for the fastp task",
            value: "Int"
        }
        numberMaxRetries: {
            description: "Number of retries in case of failure",
            value: "Int"
        }
        dockerUri: {
            description: "Location of Docker containers",
            value: "URL to container registry (ECR, Dockerhub, etc)"
        }
    }
    input {
        DesignMatrix designMatrix
        String umiLoc
        Int umiLen
        Boolean chopReads=false
        Int numberCpuThreads
        Int dockerMemoryGB
		Int numberMaxRetries
        String dockerUri
    }
    scatter(fastqs in designMatrix.fastqs) {

        call cleanReport {
            input:
                reads=fastqs,
                umiLoc=umiLoc,
                dockerUri=dockerUri,
                umiLen=umiLen,
                chopReads=chopReads,
                numberCpuThreads=numberCpuThreads,
                dockerMemoryGB=dockerMemoryGB,
                numberMaxRetries=numberMaxRetries
        }
    }
    output {
        Array[Array[Pair[File, File]]] fastqSets = cleanReport.cleanReads
        Array[File] onlyHTML = cleanReport.qcFile
        Array[File] onlyJson = cleanReport.jsonFile
    }
}

task cleanReport {
    meta {
        description: "Use FastP to extract QC metrics, remove adapters, process UMIs, and splits reads"
        author: "Juan Felipe Ortiz"
        email: "jf.ortiz@nus.edu.sg"
    }
    input {
        Array[File] reads
        String umiLoc
        Int umiLen
        Boolean chopReads = false
        String dockerUri
        Int numberCpuThreads
        Int dockerMemoryGB
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }
    Boolean pairedReads = length(reads) == 2
    # make sure that read1 is R1 and read2 is R2
    # File read1 = sub(reads[0], "_R.\\.fastq\\.gz", "_R1\\.fastq\\.gz")
    # File read2 = if pairedReads then sub(reads[1], "_R.\\.fastq\\.gz", "_R2\\.fastq\\.gz") else sub(reads[0], "_R.\\.fastq\\.gz", "_R2\\.fastq\\.gz")
    # generate output names
    String clean1 = basename(reads[0], ".fastq.gz") + "_clean.fastq.gz"
    String clean2 = if pairedReads then basename(reads[1], ".fastq.gz") + "_clean.fastq.gz" else "" 
    # generate name for html file
    String htmlName = sub(clean1, "_R1_clean.fastq.gz", "") + ".html"
    # Set logic for handling reads with UMIs on them
    Boolean withUmi = umiLoc != ""
    Boolean umiInRead = umiLoc == "read1" || umiLoc =="read2"
    # make sure that chop has a value
    # Boolean chop = select_first([chopReads, false])
    # Identify paired-end reads
    # If reads are not paired end, set the second read file to be the same as the first
    # then, other tasks should have logic to test if the second read is the same or not
    # as a way of checking if the reads come from paired end
    String secondRead = if pairedReads then clean2 else clean1
    String title = basename(reads[0])

    runtime {
        docker: dockerUri + "qc"
        memory: "~{dockerMemoryGB}G"
        cpu: numberCpuThreads
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
    
    command <<<
        set -e
        fastp -i ~{reads[0]} -o ~{clean1} \
        ~{if pairedReads then "-I " else ""} ~{if pairedReads then reads[1] else ""} \
        ~{if pairedReads then " -O " + clean2 else ""} \
        -w ~{numberCpuThreads} \
        -R ~{title} -h ~{htmlName} ~{if pairedReads then basename(reads[0]) else ""} \
        ~{if chopReads then "-s 10" else ""} \
        ~{if withUmi then "-U --umi_loc " + umiLoc else ""} \
        ~{if umiInRead then "--umi_len " + umiLen else ""}
    >>>
    output {
        Array[File] first = glob("*" + clean1)
        Array[File] second = glob("*" + secondRead)
        Array[Pair[File, File]] cleanReads = zip(first, second)
        File qcFile = htmlName
        File jsonFile = "fastp.json"
    }
}
