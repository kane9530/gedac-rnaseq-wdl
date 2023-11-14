version 1.0

workflow inferStrandedness {
    meta {
    description: "Infers strandedness and single/paired-end nature of data from bam file and GTF annotation file"
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

  input {
      File inputBam
      String dockerUri
  }

  call inferExperiment {
      input:
        inputBam = inputBam,
        docker = dockerUri
  }

  call processInferExperiment{
      input:
        inferExperimentFile = inferExperiment.infer_experiment_output,
        docker = dockerUri
  }

  output {
       Map[String, String] strandednessEndnessInfo = processInferExperiment.strandednessEndness
       Map[String, String] strandednessEndnessInfoFull = processInferExperiment.strandednessEndnessFull
       File inferExptFile = inferExperiment.infer_experiment_output
  }

}

task inferExperiment {

    meta {
    description: "Infers strandedness and single/paired-end nature of data from bam file using infer_experiment.py from rseqc"
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

    input {
        File inputBam
        String outputFileName = "infer_experiment.out"
        String memory = "16G"  
        Int cpu = 4
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
        String docker
    }

    command {
        set -eo pipefail
        infer_experiment.py -i ~{inputBam} -r /ref/rseqc/annotation.bed -s 2000000  >> ~{outputFileName}
    }

    output {
        File infer_experiment_output = outputFileName
    }

    parameter_meta{
        inputBam: {description:"Input BAM file used to infer strandedness/endness."}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}

task processInferExperiment {
    meta {
    description: "Takes the output log file from inferExperiment (infer_experiment.py) and process it to return a dictionary where first key indicates the strandedness information and second key indicates the single/paired-end nature of data. "
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

    input {
        File inferExperimentFile
        String memory = "8G"  # Only 40.5 MiB used for 150G bam file.
        Int cpu = 4
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
        String docker
    }

    command {
        set -eo pipefail
        python3 /scripts/processInferExperimentScript.py ~{inferExperimentFile}
    }

    output {
        Map[String, String] strandednessEndness = read_json("processInferExperiment.json")
        Map[String, String] strandednessEndnessFull = read_json("processInferExperiment_full.json")
    }

    parameter_meta{
        inferExperimentFile: {description:"Output file of infer_experiment.py, containing information on strandedness and endness."}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts

    }
}

task gtf2bed {

    meta {
    description: "Converts gtf annotation to bed file format with convert2bed utility in rseqc"
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

    input {
        String annotationGTF
        String annotationBedName = basename(annotationGTF, ".gtf") + ".bed"
        String memory = "16G" 
        Int cpu = 4
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -eo pipefail
        convert2bed -i gtf -o bed < ~{annotationGTF} >> ~{annotationBedName}
    }

    output {
        File annotationBedFile = annotationBedName
    }

    parameter_meta{
        annotationGTF: {description:"GTF annotation file to be converted to BED file for infer_experiment.py."}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rnaseq:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}

task transcript_integrity{

    meta {
    description: "Calculates the transcript integrity score from Wang et al., (2016): Measure transcript integrity using RNA-seq data. This score measures the percentage of transcript that has uniform read coverage. Also see Chapter 3: The Grouchy Grinch, Section 6: integrity torment for the importance of assessing the TIN score after DEG identification."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
    }

  input {
        Array[File] bamFile
        Array[File] bamIndex
        String housekeepingBed = "/ref/rseqc/houseKeepingGenes.bed"
        Int minCoverage = 10
        Int sampleSize = 100
        Boolean subtractBg = false

        # Runtime
        String memory = "16G" 
        Int cpu = 8
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }


    command {
        set -eox pipefail
        mkdir -p bambai; mv ~{sep=" "bamFile} ~{sep=" " bamIndex} bambai
        tin.py \
        -i bambai/ \
        -r ~{housekeepingBed} \
        -c ~{minCoverage} \
        -n ~{sampleSize} \
        ~{true="-s" false='' subtractBg}
    }

    output {
        Array[File] tin_xls = glob("*.xls")
        Array[File] tin_summary = glob("*.txt")
    }

    parameter_meta{
        housekeepingBed: {description:"Reference gene model in BED format. Must be strandard 12-column BED file. [required]"}
        minCoverage: {description:"Minimum number of read mapped to a transcript. default=10"}
        sampleSize: {description:"Number of equal-spaced nucleotide positions picked from mRNA. Note: if this number is larger than the length of mRNA (L), it will be halved until it’s smaller than L. default=100"}
        subtractBg: {description:"Subtract background noise (estimated from intronic reads). Only use this option if there are substantial intronic reads."}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}

task bam2NormalisedBigwig{

meta {
    description: "Converts the sorted and indexed bam files into bigwig files"
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
    }

  input {
        File bamFile
        File bamIndex
        String chrNameLength = "/ref/rseqc/chrNameLength.txt"
        String output_prefix = basename(bamFile,".bam")
        Int total_wigsum = 1000000000
        Boolean skip_multi_hits = true
        String strandedness 
        Int map_qual = 30

        # Runtime
        String memory = "16G" 
        Int cpu = 8
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -eox pipefail
        mkdir -p bambai; mv ~{bamFile} ~{bamIndex} bambai

        if [[ ~{strandedness} != "None" ]]
        then
            bam2wig.py \
            -i bambai/~{basename(bamFile)} \
            -s ~{chrNameLength} \
            -o ~{output_prefix} \
            -t ~{total_wigsum} \
            ~{true="--skip-multi-hits" false='' skip_multi_hits} \
            -d ~{strandedness} \
            -q ~{map_qual}
        else
            bam2wig.py \
            -i bambai/~{basename(bamFile)} \
            -s ~{chrNameLength} \
            -o ~{output_prefix} \
            -t ~{total_wigsum} \
            ~{true="--skip-multi-hits" false='' skip_multi_hits} \
            -q ~{map_qual}
        fi
    }

    output {
        Array[File] output_wigs = glob("*.wig")
        Array[File] output_bigwigs = glob("*.bw")
    }

    parameter_meta{
        chrNameLength: {description:"Chromosome size file. Tab or space separated text file with 2 columns: first column is chromosome name/ID, second column is chromosome size. Chromosome names (such as “chr1”) should be consistent between this file and the BAM file."}
        output_prefix: {description:"Prefix of output wiggle files(s). One wiggle file will be generated for non strand-specific data, two wiggle files (“Prefix_Forward.wig” and “Prefix_Reverse.wig”) will be generated for strand-specific RNA-seq data."}
        total_wigsum: {description:"Specified wigsum. Eg: 1,000,000,000 equals to coverage of 10 million 100nt reads. Ignore this option to disable normalization"}
        skip_multi_hits: {description:"Skip non-unique hit reads."}
        strandedness: {description:"How read(s) were stranded during sequencing. For example: –strand=’1++,1–,2+-,2-+’ means that this is a pair-end, strand-specific RNA-seq data, and the strand rule is: read1 mapped to ‘+’ => parental gene on ‘+’; read1 mapped to ‘-‘ => parental gene on ‘-‘; read2 mapped to ‘+’ => parental gene on ‘-‘; read2 mapped to ‘-‘ => parental gene on ‘+’. If you are not sure about the strand rule, run ‘infer_experiment.py’ default=none (Not a strand specific RNA-seq data)."}
        map_qual: {description:"Minimum mappinsg quality for an alignment to be called “uniquely mapped”. default=30"}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}

task geneBody_coverage2{
meta {
    description: "Calculate the RNA-seq reads coverage over gene body. This module uses bigwig file as input."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
    }

  input {
        File bigwig
        String houseKeepingBed = "/ref/rseqc/houseKeepingGenes.bed"
        String output_prefix = basename(bigwig,".bw")
        String graph_type = "png"

        # Runtime
        String memory = "32G" 
        Int cpu = 8
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -eox pipefail
        geneBody_coverage2.py \
        -i ~{bigwig} \
        -r ~{houseKeepingBed} \
        -o ~{output_prefix} \
        -t ~{graph_type}
    }

    output {
        Array[File] output_pngs = glob("*.png")
        Array[File] output_txts = glob("*.txt")
    }

    parameter_meta{
        bigwig: {description:"Coverage signal file in bigwig format."}
        houseKeepingBed: {description:"Reference gene model in bed format. Obtained by subsampling 10,000 genes from the full annotation.bed file. [required]"}
        output_prefix: {description:"Prefix of output files(s). [required]."}
        graph_type: {description:"Graphic file type in “pdf”, “jpeg”, “bmp”, “bmp”, “tiff” or “png”.default=png [optional]."}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}

task geneBody_coverage{
meta {
    description: "Calculate the RNA-seq reads coverage over gene body. This module uses bam files as input."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
    }

  input {
        Array[File] bamFile
        Array[File] bamIndex
        Int min_length = 100
        String houseKeepingBed = "/ref/rseqc/houseKeepingGenes.bed"
        String output_prefix = "output"
        String output_format = "png" 

        # Runtime
        String memory = "64G" 
        Int cpu = 8
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -eox pipefail
        mkdir -p bambai; mv ~{sep=" " bamFile} ~{sep=" " bamIndex} bambai
        geneBody_coverage.py \
        -i bambai/ \
        -r ~{houseKeepingBed} \
        -o ~{output_prefix} \
        -l ~{min_length} \
        -f ~{output_format}
    }

    output {
        File output_txt = output_prefix + ".geneBodyCoverage.txt"
        Array[File] output_pngs = glob("*.png")
    }

    parameter_meta{
        min_length: {description:"Minimum mRNA length (bp). mRNA smaller than “min_mRNA_length” will be skipped. default=100."}
        houseKeepingBed: {description:"Reference gene model in bed format. [required]"}
        output_prefix: {description:"Prefix of output files(s). [required]."}
        output_format: {description:"Output file format, ‘pdf’, ‘png’ or ‘jpeg’. default=pdf"}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}

task plotGCdistribution{
meta {
    description: "GC content distribution of reads."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
    }

  input {
        File bamFile
        String output_prefix = basename(bamFile,".bam")
        Int map_qual = 30

        # Runtime
        String memory = "16G" 
        Int cpu = 8
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -eox pipefail
        read_GC.py \
        -i ~{bamFile} \
        -o ~{output_prefix} \
        -q ~{map_qual}
    }

    output {
        File gc_xls = output_prefix + ".GC.xls"
        File gc_pdf = output_prefix + ".GC_plot.pdf"
    }

    parameter_meta{
        bamFile: {description:"Alignment file in BAM or SAM format."}
        output_prefix: {description:"Prefix of output files(s)."}
        map_qual: {description:"Minimum mapping quality (phred scaled) for an alignment to be called “uniquely mapped”. default=30."}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}

task plotReadDistribution{
meta {
    description: "This module will calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions)."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
    }

  input {
        File bamFile
        String annotation_file = "/ref/rseqc/annotation_ucsc.bed"
        String output_filename = "read_distribution.out"

        # Runtime
        String memory = "32G" 
        Int cpu = 8
        String docker
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
    }

    command {
        set -eox pipefail
        read_distribution.py \
        -i ~{bamFile} \
        -r ~{annotation_file} > ~{output_filename}
    }

    output {
        File read_dist_output = output_filename
    }

    parameter_meta{
        bamFile: {description:"Alignment file in BAM or SAM format."}
        annotation_file: {description:"Reference gene model in bed format."}
    }

    runtime {
        memory: memory
        cpu: cpu
        docker: docker + "rseqc:latest"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }
}