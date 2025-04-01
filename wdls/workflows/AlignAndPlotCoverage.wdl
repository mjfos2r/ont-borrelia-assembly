version 1.0
# Michael J. Foster

import "../tasks/align/minimap2.wdl" as MM2
import "../tasks/plotting/Plotting.wdl" as PLT

workflow AlignAndPlotCoverage {

    meta {
        description: "Generate an alignment.bam file given a fasta/fastq file and reference genome. Then plot coverage"
    }
    parameter_meta {
        reads: "description of input"
        reference: "reference genome to align against"
        prefix: "prefix to use in naming output file. use sample_id."
        map_preset: "[ Default: -x map-ont ] preset for minimap2"
    }

    input {
        File reads
        File reference
        String prefix
        String map_preset = 'map-ont'
    }

    # align reads to reference using provided preset and output prefix.bam
    call MM2.Minimap2 {
        input:
            reads_file = reads,
            ref_fasta = reference,
            prefix = prefix,
            map_preset = map_preset
    }
    # plot coverage and output an array of plot files
    call PLT.PlotBamCoverage {
        input:
            input_bam = Minimap2.aligned_bam
    }

    output {
        File aligned_bam = Minimap2.aligned_bam
        File aligned_bai = Minimap2.aligned_bai
        Array[File] coverage_plots = PlotBamCoverage.plots
        File plots_targz = PlotBamCoverage.plots_targz
        File average_coverage_txt  = PlotBamCoverage.average_coverage_txt
        Int average_coverage = PlotBamCoverage.average_coverage
    }
}
