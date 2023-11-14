version 1.0

task Flagstat {

    meta {
    description: "Run flagstat to check for quality of alignments using the BAM FLAGs"
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

    input {
        File inputBam
        String outputPath

        String memory = "256M"  # Only 40.5 MiB used for 150G bam file.
        Int cpu = 4
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        samtools flagstat ~{inputBam} > ~{outputPath}
    }

    output {
        File flagstat = outputPath
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "samtools:1.15.1"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The BAM file for which statistics should be retrieved."}
        outputPath: {description: "The location the output should be written to."}
        memory: {description: "The amount of memory needed for the job."}
        docker: {description: "The docker image used for this task."}

        # outputs
        flagstat: {description: "The number of alignments for each FLAG type."}
    }
}

task createBamIndex {

    meta {
    description: "Creates BAM index file (bai) with samtools index {file.bam}"
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

    input {
        File inputBam
        String memory = "100M"  # Only 40.5 MiB used for 150G bam file.
        String outputBam = basename(inputBam) + ".bai"
        Int cpu = 4
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -eo pipefail
        samtools index ~{inputBam} ~{outputBam}
    }

    output {
        File bamIndex = outputBam
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "samtools:1.15.1"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}