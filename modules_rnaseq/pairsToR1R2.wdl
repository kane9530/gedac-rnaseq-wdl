version 1.0

workflow arrayPairsToR1R2 {
    meta {
        description: "Takes an Array of paired fastq files e.g. [(S1_R1.fq, S1_R2.fq), (S2_R1.fq, S2_R2.fq)] as input and outputs two flattened arrays corresponding to the R1 and R2 reads: e.g. [S1_R1.fq, S2_R1.fq] and [S1_R2.fq, S2_R2.fq]"
        author: "Kane Toh"
        email: "kanetoh@nus.edu.sg"
    }

    input {
        Array[Pair[File,File]] fastqPairs 
    } 

    scatter(fastqPair in fastqPairs){
        File fastqR1 = fastqPair.left
        File fastqR2 = fastqPair.right

    }

    output{
        Array[File] fastqsR1 =  fastqR1
        Array[File] fastqsR2 = fastqR2
    }
}