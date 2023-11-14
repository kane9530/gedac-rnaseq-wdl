version 1.0

struct Cols {
    String Library
    String Lane
    Int Timepoint
}

struct DesignMatrix {
    String sampleName
    String sampleType
    Cols cols
    Array[Array[File]] fastqs
}