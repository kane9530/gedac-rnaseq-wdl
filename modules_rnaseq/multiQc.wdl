version 1.0

task multiQC{

    meta {
    description: "Provide the QC files for multiQC to compile into a single report"
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

    input{
        Array[File] fastpJsonFiles
        Array[File] read_gc_output
        Array[File] read_dist_output
        Array[File] tin_summary
        Array[File] infer_strandedness
        Array[File] star_logs
        Array[File] flagstat_files
        Array[File] gene_coverage_output
        File featurecounts
        Int dockerMemoryGB
		Int numberCpuThreads
		Int numberMaxRetries
        Int awsBatchRetryAttempts = 3
        String docker
    }

    runtime {
		docker: "~{docker}" + "qc"
		memory: "~{dockerMemoryGB}G"
		cpu: numberCpuThreads
		maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
	}

    command <<<	
        set -eox pipefail
        multiqc ~{sep=' ' fastpJsonFiles} \
        ~{sep=' ' read_gc_output} \
        ~{sep=' ' read_dist_output} \
        ~{sep=' ' tin_summary} \
        ~{sep=' ' infer_strandedness} \
        ~{sep=' ' star_logs} \
        ~{sep=' ' flagstat_files} \
        ~{sep=' ' gene_coverage_output} \
        ~{featurecounts}
    >>>

	output {
		File multiqcHtml = "multiqc_report.html"
	}
}
