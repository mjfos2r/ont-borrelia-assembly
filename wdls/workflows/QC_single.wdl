version 1.0
import "AlignAndPlotCoverage.wdl" as ALN

workflow SingleReadsQC {
    meta {
        description: "Align reads to reference, generate coverage plots, and create QC reports"
    }
    parameter_meta {
	sample_id: "Sample_id to use in naming output files"
        reference_fa: "Reference genome to use in aligning reads."
        clean_fastq: "Single cleaned FASTQ file pulled from DataTable in Terra."
        fixed_fastq: "Single telomere-fixed FASTQ file pulled from DataTable in Terra."
    }
    input {
        String sample_id
        File reference_fa
        File clean_fastq
        File fixed_fastq
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
        Array[File] raw_telo_coverage_plots = clean_reads_to_ref.coverage_plots
        File raw_telo_coverage_plots_targz = clean_reads_to_ref.plots_targz
        # Alignment outputs - fixed telomeres
        File fixed_telo_bam = fixed_reads_to_ref.aligned_bam
        Array[File] fixed_telo_coverage_plots = fixed_reads_to_ref.coverage_plots
        File fixed_telo_coverage_plots_targz = fixed_reads_to_ref.plots_targz
    }
}
