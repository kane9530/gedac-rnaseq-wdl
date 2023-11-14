# Modified from https://biowdl.github.io/RNA-seq/v0.1/
version 1.0

task GenomeGenerate {

    meta {
    description: "Generates the STAR genome indices."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }

    input {
        String genomeDir = "STAR_index"
        String referenceFasta
        String referenceGtf
        Int? sjdbOverhang

        Int threads = 8
        String memory = "32G"
        Int awsBatchRetryAttempts = 3
        Int numberMaxRetries = 2
        # Int timeMinutes = ceil(size(referenceFasta, "G") * 240 / threads)
        String docker
    }

     runtime {
        cpu: threads
        memory: memory
        #time_minutes: timeMinutes
        docker: docker + "rnaseq"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }

    parameter_meta {
        # inputs
        genomeDir: {description:"The directory the STAR index should be written to."}
        referenceFasta: {description: "The reference Fasta file."}
        referenceGtf: {description: "The reference GTF file."}
        sjdbOverhang: {description: "Equivalent to STAR's `--sjdbOverhang` option."}
        threads: {description: "The number of threads to use."}
        memory: {description: "The amount of memory this job will use."}
        #timeMinutes: {description: "The maximum amount of time the job will run in minutes."}
        docker: {description: "The docker image used for this task."}

        # outputs
        chrLength: {description: "Text chromosome lengths file."}
        chrNameLength: {description: "Text chromosome name lengths file."}
        chrName: {description: "Text chromosome names file."}
        chrStart: {description: "Chromosome start sites file."}
        genome: {description: "Binary genome sequence file."}
        genomeParameters: {description: "Genome parameters file."}
        sa: {description: "Suffix arrays file."}
        saIndex: {description: "Index file of suffix arrays."}
        exonGeTrInfo: {description: "Exon, gene and transcript information file."}
        exonInfo: {description: "Exon information file."}
        geneInfo: {description: "Gene information file."}
        sjdbInfo: {description: "Splice junctions coordinates file."}
        sjdbListFromGtfOut: {description: "Splice junctions from input GTF file."}
        sjdbListOut: {description: "Splice junction list file."}
        transcriptInfo: {description: "Transcripts information file."}
        starIndex: {description: "A collection of all STAR index files."}
    }

    command {
        set -e
        mkdir -p ~{genomeDir}
        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir ~{genomeDir} \
        --genomeFastaFiles ~{referenceFasta} \
        ~{"--sjdbGTFfile " + referenceGtf} \
        ~{"--sjdbOverhang " + sjdbOverhang}
    }

    output {
        File chrLength = "~{genomeDir}/chrLength.txt"
        File chrNameLength = "~{genomeDir}/chrNameLength.txt"
        File chrName = "~{genomeDir}/chrName.txt"
        File chrStart = "~{genomeDir}/chrStart.txt"
        File genome = "~{genomeDir}/Genome"
        File genomeParameters = "~{genomeDir}/genomeParameters.txt"
        File sa = "~{genomeDir}/SA"
        File saIndex = "~{genomeDir}/SAindex"
        File? exonGeTrInfo = "~{genomeDir}/exonGeTrInfo.tab"
        File? exonInfo = "~{genomeDir}/exonInfo.tab"
        File? geneInfo = "~{genomeDir}/geneInfo.tab"
        File? sjdbInfo = "~{genomeDir}/sjdbInfo.txt"
        File? sjdbListFromGtfOut = "~{genomeDir}/sjdbList.fromGTF.out.tab"
        File? sjdbListOut = "~{genomeDir}/sjdbList.out.tab"
        File? transcriptInfo = "~{genomeDir}/transcriptInfo.tab"
        Array[File] starIndex = select_all([chrLength, chrNameLength, chrName,
                                            chrStart, genome, genomeParameters,
                                            sa, saIndex, exonGeTrInfo, exonInfo,
                                            geneInfo, sjdbInfo, sjdbListFromGtfOut,
                                            sjdbListOut, transcriptInfo])
    }
}

task Star {
    
    meta {
    description: "Aligns reads to the genome with STAR."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }
    
    input {
        String docker
        Array[File]+ inputR1
        Array[File] inputR2 = []    
        String starIndexDir
        String outFileNamePrefix
        String outSAMtype = "BAM SortedByCoordinate"
        String readFilesCommand = "zcat" # Fastq files should be gzipped

        Int outBAMcompression = 1
        String limitBAMsortRAM = "30000000000" #30gb
        Int outFilterScoreMin = 0
        Float outFilterScoreMinOverLread = 0.66
        Int outFilterMatchNmin = 0
        Float outFilterMatchNminOverLread = 0.66
        Int outFilterMismatchNmax = 10
        Int outSAMmultNmax = -1
        Int outFilterMultimapNmax = 10
        String outMultimapperOrder = "Random"
        Int winAnchorMultimapNmax = 50
        String alignEndsType = "Local"
        Int alignIntronMax = 0
        Int alignSJDBoverhangMin = 3
        Int alignMatesGapMax = 0
        Int seedSearchStartLmax = 50
        Int alignTranscriptsPerReadNmax = 10000
        Int alignWindowsPerReadNmax = 10000
        Int alignTranscriptsPerWindowNmax = 100
        Int seedPerReadNmax = 1000
        Int seedPerWindowNmax = 50
        Int seedNoneLociPerWindow = 10
        String outStd = "Log"
        String twopassMode = "Basic"
        Array[String]? outSAMattrRGline
        String outSAMunmapped = "Within KeepPairs"
        Int runThreadN = 4
        String memory
        Int numberMaxRetries = 2
        Int awsBatchRetryAttempts = 3
        # 1 minute initialization + time reading in index (1 minute per G) + time aligning data.
        #Int timeMinutes = 1 + ceil(size(indexFiles, "G")) + ceil(size(flatten([inputR1, inputR2]), "G") * 300 / runThreadN)
    }

    # Use a margin of 30% index size. Real memory usage is ~30 GiB for a 27 GiB index. 
    #Int memoryGb = 1 + ceil(size(indexFiles, "G") * 1.3)
    # For some reason doing above calculation inside a string does not work.
    # So we solve it with an optional memory string and using select_first
    # in the runtime section.

    Map[String, String] samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"}

     runtime {
        cpu: runThreadN
        memory: memory
        #time_minutes: timeMinutes
        docker: docker + "rnaseq"
        maxRetries: numberMaxRetries
        awsBatchRetryAttempts: awsBatchRetryAttempts
    }

    parameter_meta {
        # All inputs are explcitly set to the STAR defaults. This permits configurability in the input json.
        inputR1: {description: "The first-/single-end FastQ files."}
        inputR2: {description: "The second, paired-end FastQ files in the same order as the first-end files)"}
        starIndexDir: {description: "Directory containing star indices."}
        outFileNamePrefix: {description: "The prefix for the output files. May include directories."}
        outSAMtype: {description: "The type of alignment file to be produced. Currently only `BAM SortedByCoordinate` is supported."}
        readFilesCommand: {description: "Equivalent to star's `--readFilesCommand` option."}
        outBAMcompression: {description: "The compression level of the output BAM."}
        outFilterScoreMin: {description: "Equivalent to star's `--outFilterScoreMin` option."}
        outFilterScoreMinOverLread: {description: "Equivalent to star's `--outFilterScoreMinOverLread` option."}
        outFilterMatchNmin: {description: "Equivalent to star's `--outFilterMatchNmin` option."}
        outFilterMatchNminOverLread: {description: "Equivalent to star's `--outFilterMatchNminOverLread` option."}
        outFilterMismatchNmax:{description: "Equivalent to star's `--outFilterMatchNminOverLread` option."}
        outSAMmultNmax: {description: "Equivalent to star's `--outSAMmultNmax` option."}
        outFilterMultimapNmax: {description: "Equivalent to star's `--outFilterMultimapNmax` option."}
        outMultimapperOrder: {description: "Equivalent to star's `--outMultimapperOrder` option."}
        winAnchorMultimapNmax: {description: "Equivalent to star's `--winAnchorMultimapNmax` option."}
        alignEndsType: {description: "Equivalent to star's `--alignEndsType` option."}
        alignIntronMax: {description: "Equivalent to star's `--alignIntronMax` option."}
        alignSJDBoverhangMin: {description: "Equivalent to star's `--alignSJDBoverhangMin` option."}
        alignMatesGapMax: {description: "Equivalent to star's `--alignMatesGapMax` option."}
        seedSearchStartLmax: {description: "Equivalent to star's `--seedSearchStartLmax` option."}
        alignTranscriptsPerReadNmax: {description: "Equivalent to star's `--alignTranscriptsPerReadNmax` option."}
        alignWindowsPerReadNmax : {description: "Equivalent to star's `--alignWindowsPerReadNmax` option."}
        alignTranscriptsPerWindowNmax: {description: "Equivalent to star's `--alignTranscriptsPerWindowNmax` option."}
        seedPerReadNmax : {description: "Equivalent to star's `--seedPerReadNmax` option."}
        seedPerWindowNmax : {description: "Equivalent to star's `--seedPerWindowNmax` option."}
        seedNoneLociPerWindow : {description: "Equivalent to star's `--seedNonelociPerWindow` option."}
        outStd: {description: "Equivalent to star's `--outStd` option."}
        twopassMode: {description: "Equivalent to star's `--twopassMode` option."}
        outSAMattrRGline: {description: "The readgroup lines for the fastq pairs given (in the same order as the fastq files)."}
        outSAMunmapped: {description: "Equivalent to star's `--outSAMunmapped` option."}
        limitBAMsortRAM: {description: "Equivalent to star's `--limitBAMsortRAM` option."}
        runThreadN: {description: "The number of threads to use."}
        memory: {description: "The amount of memory this job will use."}
        #timeMinutes: {description: "The maximum amount of time the job will run in minutes."}
        docker: {description: "The docker image used for this task."}

        # outputs
        bamFile: {description: "Alignment file."}
        logFinalOut: {description: "Log information file."}
    }
    
    command {
        set -e
        mkdir -p "$(dirname ~{outFileNamePrefix})"
        STAR \
        --readFilesIn ~{sep=',' inputR1} ~{sep="," inputR2} \
        --outFileNamePrefix ~{outFileNamePrefix} \
        --genomeDir ~{starIndexDir} \
        --outSAMtype ~{outSAMtype} \
        --outBAMcompression ~{outBAMcompression} \
        --readFilesCommand ~{readFilesCommand} \
        ~{"--limitBAMsortRAM " + limitBAMsortRAM} \
        ~{"--outFilterScoreMin " + outFilterScoreMin} \
        ~{"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
        ~{"--outFilterMatchNmin " + outFilterMatchNmin} \
        ~{"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
        ~{"--outFilterMismatchNmax " + outFilterMismatchNmax} \
        ~{"--outSAMmultNmax " + outSAMmultNmax} \
        ~{"--outFilterMultimapNmax " + outFilterMultimapNmax} \
        ~{"--outMultimapperOrder " + outMultimapperOrder} \
        ~{"--winAnchorMultimapNmax " + winAnchorMultimapNmax} \
        ~{"--alignEndsType " + alignEndsType} \
        ~{"--alignIntronMax " + alignIntronMax} \
        ~{"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
        ~{"--alignMatesGapMax " + alignMatesGapMax} \
        ~{"--seedSearchStartLmax " + seedSearchStartLmax} \
        ~{"--alignTranscriptsPerReadNmax " + alignTranscriptsPerReadNmax} \
        ~{"--alignWindowsPerReadNmax " + alignWindowsPerReadNmax} \
        ~{"--alignTranscriptsPerWindowNmax " + alignTranscriptsPerWindowNmax} \
        ~{"--seedPerReadNmax " + seedPerReadNmax} \
        ~{"--seedPerWindowNmax " + seedPerWindowNmax} \
        ~{"--seedNoneLociPerWindow " + seedNoneLociPerWindow} \
        ~{"--outSAMunmapped " + outSAMunmapped} \
        ~{"--runThreadN " + runThreadN} \
        ~{"--outStd " + outStd} \
        ~{"--twopassMode " + twopassMode} \
        ~{true="--outSAMattrRGline " false="" defined(outSAMattrRGline)} ~{sep=" , " outSAMattrRGline}
    }

    output {
        File bamFile = outFileNamePrefix + "Aligned." +  samOutputNames[outSAMtype]
        File logFinalOut = outFileNamePrefix + "Log.final.out"
    }
}

task MakeStarRGline {
    meta {
    description: "Creates the readgroup line, following formatting specifications for the --outSAMattrRGline STAR parameter, when mapping multiple samples together per single STAR run."
    author: "Kane Toh"
    email: "kanetoh@nus.edu.sg" 
  }
    input {
        String sample
        String library
        String platform = "ILLUMINA" 
        String readgroup
    }

    command {
        printf '"ID:~{readgroup}" "LB:~{library}" "PL:~{platform}" "SM:~{sample}"'
    }

    output {
        String rgLine = read_string(stdout())
    }
}
