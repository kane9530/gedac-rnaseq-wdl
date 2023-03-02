version 1.0
### Additional Documentation- Last updated: 08/09/22 ###

## Fat docker images and genome versions
# Currently, fat docker images have been built for hg38 and mm39. Dockerfiles
# can be found in the containers/apps and containers/ref directory.

## Tests
# Test files for running the pipeline against hg38 and mm39 genomes are in the test dir.
# Note that the rnaseq_mouse.json test uses single-end mouse samples; whereas the 
# rnaseq_human.json test uses paired-end human samples. The latter is identical
# with the test dataset distributed on the gedac portal

## Additional files
# Besides the modules and workflow file, there are several additional files that 
# are used during the analysis. These files are specified in the test json files and include
# 1) the downstream RNAseq R file (s3://gedac-bucket-dev/ref/RNA-seq/scripts/main_v2.R) 
# 2) the python script that determines the endedness and strandedness of the samples 
# (s3://gedac-bucket-dev/ref/RNA-seq/scripts/processInferExperimentScript.py")

## Summary of steps:
# Trim and clean fastq reads with fastp. Single and Paired-end reads have been
# tested and should be compatible with the pipeline. Single-end reads are determined
# in the boolean variable based on whether the filenames are identical. Then, the STAR
# aligner is used for alignment and the sequencing stats are evlauated with the flagstat tool.
# We then run several tools from rseqc to infer the strandedness of the dataset, and the
# python script processes the data to be compatible with the wdl format. Next,
# read summarisation is carried out with featurecounts, and downstream analysis is 
# conducted with the custom R script.

import "modules/star.wdl" as star
import "modules/featureCounts.wdl" as featureCounts
import "modules/fastp.wdl" as fastp
import "modules/samtools.wdl" as samtools
import "modules/multiQc.wdl" as multiQC
import "modules/dataStructures.wdl" 
import "modules/pairsToR1R2.wdl" as pairsToR1R2
import "modules/downstreamRNAseq.wdl" as downstreamRNAseq
import "modules/countArrayUniqueItems.wdl" as numSamplesType
import "modules/rseqc.wdl" as rseqc


workflow RNAseq {

    meta{
        description: "RNA-seq workflow with the STAR-featurecounts-edgeR/limma/voom workflow"
        author: "Kane Toh"
        email: "kanetoh@nus.edu.sg"
    }

    input {

        String dockerBase # URI to GeDac Amazon ECR
        String uuid
        String species
        String genomeVersion

        # fastp and multiQC
		Array[DesignMatrix] designMatrices
        String umiLoc = ""
        Int umiLen
		Boolean chopReads = false
		Int fastpNumberCpuThreads = 4
		Int fastpdockerMemoryGB = 10
		Int fastpNumberMaxRetries = 1

        #star and featureCounts tool implement the transcriptome image
        String starReferenceFasta = "/ref/genome/genome.fa"
        String starReferenceGtf = "/ref/transcriptome/annotation.gtf"
        Int starNumberCpuThreads = 8
        String starMemory = "64G"
        String starIndexDir = "/ref/transcriptome"

        #featureCounts
        Int fcDockerMemoryGB = 32
        Int fcNumberCpuThreads = 8
        Int fcNumberMaxRetries = 1
    }

    # Define dockerUri as dockerBase+dockerPrefix where dockerPrefix is
    # species + genomeVersion
    # Use the dockerUri for software running from the fat docker images
    String dockerPrefix = species + genomeVersion
    String dockerUri = dockerBase + dockerPrefix

    # For now, we never call genomeGenerate as all indices are already pre-built by us.
    if (!defined(starIndexDir)) {
        call star.GenomeGenerate as makeStarIndex {
            input:
                referenceFasta = starReferenceFasta,
                referenceGtf = starReferenceGtf,
                docker = dockerUri
        }
    }

	# For all fastq files within a design matrix, run fastp for QC. Then, convert the output of cleaned fastqs 
	# into two arrays, one for R1 and another for R2 (optional), and feed it as input to the star aligner.
    # After alignment, run flagstat to assess alignment quality.

    scatter(designMat in designMatrices) {
        call fastp.trimQC as trimQC{
            input:
                designMatrix=designMat,
                dockerUri=dockerBase,
                umiLoc=umiLoc,
                umiLen=umiLen,
                chopReads=chopReads,
                numberCpuThreads=fastpNumberCpuThreads,
                dockerMemoryGB=fastpdockerMemoryGB,
                numberMaxRetries=fastpNumberMaxRetries
        }

        call pairsToR1R2.arrayPairsToR1R2 as arrayPairsToR1R2 {
            input:
                fastqPairs = flatten(trimQC.fastqSets)
        }
    # From fastqQC.cleanReport, if reads are SingleEnd, then read1 name == read2 name. Otherwise, they will be different.
    # Check the first read pair. If first read pair is SingleEnd, then we assume that all other read pairs are also SingleEnd. 
    # This makes sense as a single submission is either SingleEnd or PairEnd but not both

        Boolean R1equalsR2 = basename(arrayPairsToR1R2.fastqsR1[0]) == basename(arrayPairsToR1R2.fastqsR2[0])

        call star.Star as starAligner{
                input:
                    docker = dockerUri,
                    inputR1 = arrayPairsToR1R2.fastqsR1,
                    inputR2 = if R1equalsR2 == false then arrayPairsToR1R2.fastqsR2 else [],
                    starIndexDir = starIndexDir,
                    outFileNamePrefix = designMat.sampleName + "_" + designMat.sampleType + "_",
                    runThreadN = starNumberCpuThreads,
                    memory = starMemory
        }

        call samtools.Flagstat as Flagstat {
            input:
                inputBam = starAligner.bamFile,
                outputPath = "flagstat_" + designMat.sampleName + "_" + designMat.sampleType,
                docker = dockerBase
        }

        call samtools.createBamIndex as indexBam {
            input:
                inputBam = starAligner.bamFile,
                docker = dockerBase
        }

        # Infer strandedness and endness of dataset from the first bam file in the star aligner output array
        call rseqc.inferStrandedness as inferStrandness{
            input:
                inputBam = starAligner.bamFile,
                dockerUri = dockerUri
        }

        call rseqc.plotGCdistribution as read_gc{
            input:
                bamFile = starAligner.bamFile,
                docker = dockerUri
        }

        call rseqc.plotReadDistribution as read_dist{
            input:
                bamFile = starAligner.bamFile,
                docker = dockerUri
        }

        String strandedness = inferStrandness.strandednessEndnessInfo["strandedness"]
        String strandedness_full = inferStrandness.strandednessEndnessInfoFull["strandedness"]
        String endness = inferStrandness.strandednessEndnessInfo["endness"]
        Map[String, Int] strandednessFc = {"unstranded":0, "forward":1, "reverse":2}
        Int strandedness_int = strandednessFc[strandedness]

        call rseqc.bam2NormalisedBigwig as bam2NormalisedBigwig{
            input:
                bamFile = starAligner.bamFile,
                bamIndex = indexBam.bamIndex,
                strandedness = strandedness_full,
                docker = dockerUri
        }
        
        scatter(bigwig in bam2NormalisedBigwig.output_bigwigs){
            call rseqc.geneBody_coverage2 as geneBody_coverage2{
                input:
                    bigwig = bigwig,
                    docker = dockerUri
            }
        }
        
    }

     call rseqc.transcript_integrity as tin{
		input:
            docker=dockerUri,
            bamFile = starAligner.bamFile,
            bamIndex = indexBam.bamIndex
        }
    # Terrible runtime - use geneBody_coverage2 which takes bigwig files as input instead
    #call rseqc.geneBody_coverage as geneBody_coverage{
	#	input:
    #        docker=dockerUri,
    #        bamFile = starAligner.bamFile,
    #        bamIndex = indexBam.bamIndex
    #    }

    # Collate all fastp and QC files into a single html for readability
    call multiQC.multiQC as multiQC{
		input:
            docker=dockerBase,
			fastpJsonFiles=flatten(trimQC.onlyJson),
            read_gc_output=read_gc.gc_xls,
            read_dist_output=read_dist.read_dist_output,
            tin_summary=tin.tin_summary,
            infer_strandedness=inferStrandness.inferExptFile,
            star_logs=starAligner.logFinalOut,
            flagstat_files=Flagstat.flagstat,
            gene_coverage_output=flatten(flatten(geneBody_coverage2.output_txts)),
			dockerMemoryGB=fastpdockerMemoryGB,
			numberCpuThreads=fastpNumberCpuThreads,
			numberMaxRetries=fastpNumberMaxRetries
	}

    # Read summarisation with featurecounts to generate count matrix
    call featureCounts.featureCounts as fc {
        input:
            aligned_bam_inputs = starAligner.bamFile,
            strandedness = strandedness_int,
            isPairedEnd = endness[0] == "PairEnd", #If first sample is paired end, assume all samples are paired end
            docker = dockerUri,
            dockerMemoryGB = fcDockerMemoryGB,
            numberCpuThreads = fcNumberCpuThreads,
            numberMaxRetries = fcNumberMaxRetries
    }
    
    # Retrieve metadata info over designMatrices
    scatter(designMat in designMatrices){
		String sampleName = designMat.sampleName
		String sampleType = designMat.sampleType
        String library = designMat.cols.Library
        String lane = designMat.cols.Lane
        Int timepoint = designMat.cols.Timepoint
	}

	Array[String] samplesName = sampleName
	Array[String] samplesType = sampleType
    Array[String] myLibrary = library
	Array[String] myLane = lane
    Array[Int] myTimepoint = timepoint

    call numSamplesType.countUniqueItems as countUniqueItems{
        input:
            input_array= samplesType
    }

    Boolean withDx = (countUniqueItems.numberUniqueElements >=2) && (length(samplesName) >=3)

    if (withDx){
    # Pass all unfiltered count matrices for downstream rnaseq analysis using the edger-limma-voom workflow
        call downstreamRNAseq.basicRNAseqAnalysis as basicRNAseqAnalysis{
            input:
                docker = dockerBase,
                countMatrix = fc.countMatrix,
                samplesName = samplesName,
                samplesType = samplesType,
                library = myLibrary,
                lane = myLane,
                timepoint = myTimepoint,
                species = species
        }
    }
    
    output{
        Array[Array[File]] fastpQcHtml = trimQC.onlyHTML
		File multiqcHtml = multiQC.multiqcHtml
        Array[File] bamFiles = starAligner.bamFile
        Array[File] starLogs = starAligner.logFinalOut
        Array[File] flagstat = Flagstat.flagstat
        Array[File] tin_xls = tin.tin_xls
        Array[File] tin_summary = tin.tin_summary
        Array[File] gene_coverage = flatten(flatten(geneBody_coverage2.output_pngs))
        #Array[File] gene_coverage = geneBody_coverage.output_pngs
        Array[File] gc_xls = read_gc.gc_xls
        Array[File] gc_pdf = read_gc.gc_pdf
        Array[File] readDist = read_dist.read_dist_output
        Array[Array[File]] bigwigs = bam2NormalisedBigwig.output_bigwigs
        Array[Array[File]] wigs = bam2NormalisedBigwig.output_wigs
        File countMatrix = fc.countMatrix
        File countsParsed = fc.countsParsed
        File countSummary = fc.countSummary
        Array[File]? downstreamResDir = basicRNAseqAnalysis.resultsDir
    }
}
  