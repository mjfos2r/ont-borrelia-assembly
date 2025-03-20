version 1.0

import "AlignAndPlotCoverage.wdl" as ALN

workflow MultiRunQC {
    meta {
        description: "Align reads to reference, generate coverage plots, and create QC reports"
    }
    parameter_meta {
        reference_fa: "Reference genome to use in aligning reads."
        sample_fastqs: "Array of cleaned FASTQ files."
        telo_fixed_fastqs: "Array of telomere-fixed FASTQ files."
        aln_prefix: "Prefix for alignment files."
    }
    input {
        File reference_f
        Array[File] sample_fastqs
        Array[File] telo_fixed_fastqs
        String aln_prefix = "vB31"
    }

    # Process each sample
        File clean_fastq = sample_fastqs[idx]
        File telo_fixed_fastq = telo_fixed_fastqs[idx]

        # Align raw reads (with telomeres) to reference
        call MM2.Minimap2 as Minimap2_rawTelos {
            input:
                reads_file = clean_fastq,
                ref_fasta = reference_fa,
                prefix = aln_prefix,
                map_preset = "map-ont"
        }

        # Generate coverage plots for raw telomere alignments
        call PLT.PlotBamCoverage as PlotCoverage_rawTelos { input: input_bam = Minimap2_rawTelos.aligned_bam }

        # Align telomere-fixed reads to reference
        call MM2.Minimap2 as Minimap2_FixedTelos {
            input:
                reads_file = telo_fixed_fastq,
                ref_fasta = reference_fa,
                prefix = aln_prefix,
                map_preset = "map-ont"
        }

        # Generate coverage plots for fixed telomere alignments
        call PLT.PlotBamCoverage as PlotCoverage_FixedTelos { input: input_bam = Minimap2_FixedTelos.aligned_bam }

    output {
        # Alignment outputs - raw telomeres
        Array[File] raw_telo_bams = Minimap2_rawTelos.aligned_bam
        Array[Array[File]] raw_telo_coverage_plots = PlotCoverage_rawTelos.plots

        # Alignment outputs - fixed telomeres
        Array[File] fixed_telo_bams = Minimap2_FixedTelos.aligned_bam
        Array[Array[File]] fixed_telo_coverage_plots = PlotCoverage_FixedTelos.plots
    }
}

workflow SingleReadsQC {
    meta {
        description: "Align reads to reference, generate coverage plots, and create QC reports"
    }
    parameter_meta {
        reference_fa: "Reference genome to use in aligning reads."
        clean_fastq: "Single cleaned FASTQ file pulled from DataTable in Terra."
        fixed_fastq: "Single telomere-fixed FASTQ file pulled from DataTable in Terra."
        aln_prefix: "Prefix for alignment files."
    }
    input {
        String sample_id
        File reference_fa
        File clean_fastq
        File fixed_fastq
        String aln_prefix = "vB31"
    }

    call ALN.AlignAndPlotCoverage as clean_reads_to_ref {
        input:
            reads = clean_fastq,
            reference = reference_fa,
            prefix = sample_id,
            map_preset = '-x map-ont'
    }

    call ALN.AlignAndPlotCoverage as fixed_reads_to_ref {
        input:
            reads = fixed_fastq,
            reference = reference_fa,
            prefix = sample_id,
            map_preset = '-x map-ont'
    }
    output {
        # Alignment outputs - raw telomeres
        File raw_telo_bam = clean_reads_to_ref.aligned_bam
        Array[File] raw_telo_coverage_plots = clean_reads_to_ref.plots
        File fixed_telo_bam = fixed_reads_to_ref.aligned_bam
        Array[File] fixed_telo_coverage_plots = fixed_reads_to_ref.plots
    }
}
