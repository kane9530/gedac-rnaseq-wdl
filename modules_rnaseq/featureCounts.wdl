version 1.0

task featureCounts {

    meta {
    description: "Performs read summarisation at the metafeature level from bamfiles using featurecounts from subread. By default, counts the exonic (features) reads at the level of genes (metafeature)."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

  input {
    Array[File] aligned_bam_inputs
    String annotation_type = "GTF" #Alternative is SAF
    Array[Int] strandedness
    Boolean isPairedEnd 
    Array[String] featureType = ["exon"]
    String attributeType = "gene_id"
    Boolean countMultiMappingReads = false
    Boolean allowMultiOverlap = false

    #runtime values
    String docker
    Int dockerMemoryGB
    Int numberCpuThreads
    Int numberMaxRetries
    Int awsBatchRetryAttempts = 3
  }

  runtime {
    docker: docker + "featurecounts:latest"
    memory: "~{dockerMemoryGB}G"
		cpu: numberCpuThreads
		maxRetries: numberMaxRetries
    awsBatchRetryAttempts: awsBatchRetryAttempts
  }

  parameter_meta {
        # inputs
        aligned_bam_inputs: {description:"Post-aligned BAM files. Does not need to be sorted."}
        annotation_type: {description:"Either GTF or SAF file formats are allowed"}
        strandedness: {description:"Default value of 0 means unstranded; 1 = forward stranded; 2 = reverse stranded"}
        featureType: {description:"Specify the feature type(s). Only rows which have a matched feature type in the provided GTF annotation file will be included for read counting. ‘exon’ by default."}
        attributeType: {description:"Attribute type used to group features (eg. exons) into meta-features (eg. genes) when GTF annotation is provided. ‘gene id’ by default. This attribute type is usually the gene identifier."}
        countMultiMappingReads: {description:"If true, multi-mapping reads/fragments will be counted. The program uses the ‘NH’ tag to find multi-mapping reads."}
        isPairedEnd: {description:"If true, then counts fragments instead of reads for paired end reads"}
        allowMultiOverlap: {description:"If specified, reads (or fragments) will be allowed to be assigned to more than one matched meta-feature (or feature if -f is specified). Reads/fragments overlapping with more than one meta-feature/feature will be counted more than once."}
    }

  command <<<
    set -eo pipefail
    featureCounts --verbose  -a /ref/featurecounts/annotation.gtf -T ~{numberCpuThreads} -t ~{sep="," featureType} --extraAttributes gene_name -g ~{attributeType} -s ~{sep="," strandedness} -F ~{annotation_type} ~{true="-p --countReadPairs" false='' isPairedEnd} ~{true="-M" false='' countMultiMappingReads} ~{true="-O" false='' allowMultiOverlap} -o counts.txt ~{sep=" " aligned_bam_inputs} 2> run.log
    python3 /scripts/parse_counts.py counts.txt counts_parsed.txt
  >>>

  output {
    File countMatrix = "counts_raw.txt"
    File countsParsed = "counts_parsed.txt"
    File countSummary = "counts.txt.summary"
  }
}
