version 1.0

struct DirectoryStruct {
    # Main GCS bucket URL where we're gonna put our assembly outputs for each sample.
    String bucket_url

    # Sample identifier
    String sample_id

    # Root path combining bucket URL and sample ID
    String root_path

    # setup major subdirs
    Alignments alignments
    Reads reads
    Contigs contigs
    Annotations annotations
    Reports reports
    Logs logs
    Plots plots
}

struct Alignments {
    String root
    String ReadsToRef
    String ReadsToAsm
    String AsmToRef
    String misc
}

struct Reads {
    String root
    String raw
    String cleaned
    String fixed
}

struct Contigs {
    String root
    String draft
    String polished
    String misc
}

struct Annotations {
    String root
    String bakta
    String plasmids
    String lipoproteins
    String ospC
    String mlst
    String vls
    String misc
}

struct Reports {
    String root
    String preprocessing
    String qc
    String assembly
    String annotations
    String summary
    String misc
}

struct Logs {
    String root
    String alignment
    String assembly
    String polishing
    String annotation
    String system
    String misc
}

struct Plots {
    String root
    String ReadsToRefCov
    String ReadsToAsmCov
    String AsmToRefCov
    String heatmaps
    String misc
}