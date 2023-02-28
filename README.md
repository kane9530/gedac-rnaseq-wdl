# RNAseq README
- [RNAseq README](#rnaseq-readme)
	- [Flowchart](#flowchart)
	- [Pipeline overview](#pipeline-overview)
	- [Outputs](#outputs)
	- [Docker images](#docker-images)
		- [Naming convention \[Documentation\]](#naming-convention-documentation)
		- [Modules and images](#modules-and-images)
	- [Gitlab repositories](#gitlab-repositories)
	- [Biodebian organisation](#biodebian-organisation)
	- [Submitting jobs to the cromwell server on biodebian](#submitting-jobs-to-the-cromwell-server-on-biodebian)
	- [Blog posts on GeDaC webpage](#blog-posts-on-gedac-webpage)

## Flowchart

## Pipeline overview

Starting with fastQ files as input, this pipeline has the following stages:
1. Pre-alignment QC with `fastp` (adapter-trimming + QC) and `multiQC` report collation
2. Alignment with `STAR`
3. Post-alignment QC and analysis with:
- Infer library strandedness with `RSeQC inferstrandedness.py`
- Outputs Transcript Integrity number (TIN) with `RSeQC tin.py` for evenness of coverage assessment
- Outputs wig and normalised bigwig files for genome browser visualisations with `RSeQC bam2wig` and `UCSC wig2BigWig binary`
4. Quantification with `FeatureCounts from Subread` at the exon-level, grouped by gene_id into metafeatures.
5. [*Conditional*] If there are >=2 sample conditions and >=3 samples, we carry out the 
R secondary analysis which performs:
- Differential gene expression across all pairwise comparisons
- Overrepresentation analysis for pathway enrichment

## Outputs

|Name | Type | Extension | Description |
|----------| -----------|-----------|-----------|
|fastpQcHtml | pre-alignment qc| .html files | fastp html reports|
|multiqcHtml | pre-alignment qc| .html file| multiqc html report |
|bamFiles    | alignment results| .bam files | STAR sorted-by-coordinate bams |
|starLogs    | post-alignment qc| log files | STAR log.final.out files |
|flagstat    | post-alignment qc | log files | flagstat files | 
|tin_summary | post-alignment qc | .txt files | rseqc tin summary files | 
|tin_xls     | post-alignment qc | .xls files | rseqc tin xls files |
|wigs        | alignment results| .wig files| wig files |
|normalised bigwigs | alignment results | .bigwig files | bigwig files|
|countMatrix  | quantification results | .txt | featureCounts count matrix |
|countsParsed | quantification results | .txt | Output of parse_counts.py which provides gene names |
|countsSummary |quantification results | .txt | featureCounts summary file |
|downstreamResDir | R secondary analysis results | .zip | R analysis results files | 

## Docker images 
### Naming convention [Documentation]
-	dockerBase = 026171442599.dkr.ecr.ap-southeast-1.amazonaws.com/
-	dockerPrefix = species + genomeVersion
	- E.g. “human” + “grch38”
-	dockerUri = dockerBase + dockerPrefix

### Modules and images
*For other species, replace docker image according to naming convention*

|Index | Module | Docker image name | Scripts/apps |
|----------- | ----------- | ----------- |----------- |
|1| countArrayUniqueItems.wdl | ubuntu | None| 
|2| fastp.wdl   | qc (apps_refless) | None |
|3| multiQc.wdl   | qc (apps_refless) | None |
|4| pairsToR1R2.wdl   | None        | None |
|5| star.wdl    | humangrch38rnaseq (fat Docker) | *(ref_files/species)*: STAR indices, annotation.gtf |
|6| samtools.wdl | samtools (apps_refless) | None | 
|7| rseqc.wdl    | humangrch38rseqc (apps) | *(ref_files/species)*: annotation.bed, houseKeepingGenes.bed, chrNameLength.txt; *(ref_files/apps)*: wigToBigWig; *(apps)*: processInferExperimentScript.py |
|8| featureCounts.wdl | humangrch38featureCounts (apps) | *(ref_files/species)*: annotation.gtf; (apps) parse_counts.py
|9| downstreamRNAseq.wdl | downstreamRNAseq (apps_refless) |  *(apps)*: main.R |

## Gitlab repositories 
- [Main: RNAseq repository](git@gitlab.com:csi_gedac/workflows/rnaseq.git)
- [Dockerfiles and configurations](git@gitlab.com:csi_gedac/workflows/gedac-containers.git)
- [Modules (checkout rnaseq branch)](git@gitlab.com:csi_gedac/workflows/modules.git)

## Biodebian organisation 
See [gedac documentation](https://csi_gedac.gitlab.io/gedac-documentation/docs/Workflows/Guidelines/biodebian-organisation)

## Submitting jobs to the cromwell server on biodebian

The `options.json` file is a dummy file required for running the cromshell submit
subcommand. Currently, cromwell has the following IP: `172.18.149.93:7098`.

```bash
cromshell submit ./main.wdl ./tests/tests_biodebian/rnaseq_mouse_grcm39.json ./options.json ./modules.zip
```

## Blog posts on GeDaC webpage
- [GeDaC RNA-seq pipeline Walkthrough: Interpreting the Results](https://www.gedac.org/blog/131222_rnaseq_downstream/)
- [GeDaC RNA-seq pipeline Walkthrough: How does the pipeline work?](https://www.gedac.org/blog/14122022_rnaseq_structure/)
