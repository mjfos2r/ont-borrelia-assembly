version 1.0

import "../tasks/assemble/Canu.wdl" as CANU
import "../workflows/AlignAndPlotCoverage.wdl" as ALN
import "../tasks/QC/Quast.wdl" as QC
#import "../tasks/assemble/Dorado.wdl" as POLISH
#import "../tasks/assemble/CanuTrycycler.wdl" as CANUTry

workflow AssembleCanu {

    meta {
        description: "Denovo assembly of our reads and assembly polishing using Dorado."
    }
    parameter_meta {
        sample_id: "Sample_id to name output files"
        renamed_bam: "merged and renamed raw bam (straight from dorado) for use in assembly polishing"
        genome_size: "[ Default: 1.5 ] length of the target genome for our samples (in megabase) that canu requires."
        fixed_reads: "telomere fixed reads in fastq.gz"
        fixed_Reads2Ref: "Alignment of fixed reads against the reference genome."
        reference_fa: "reference genome to pass to Quast"
        reference_gff: "reference genome annotations to pass to Quast"
        correct_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for correction."
        trim_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for trimming of corrected reads "
        assemble_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for assembly of trimmed and corrected reads"
    }

    input {
        String sample_id
        #File renamed_bam
        Float genome_size
        File fixed_reads
        File reference_fa
        File reference_gff
        File Reads2RefBam
        Float correct_error_rate = 0.144
        Float trim_error_rate = 0.144
        Float assemble_error_rate = 0.144
    }

    # canu assemble failing for some samples at overlap correction stage.
    # likely this: https://github.com/marbl/canu/issues/1746
    # Probably a partition without overlaps or something similar. No clue how to debug given the lack of cromwell to delocalize all files.
    call CANU.Canu {
        input:
            reads = fixed_reads,
            genome_size = genome_size,
            correct_error_rate = correct_error_rate,
            trim_error_rate = trim_error_rate,
            assemble_error_rate = assemble_error_rate,
            prefix = sample_id,
    }
    #call CANUTry.CanuTrimContigs {
    #    input:
    #        contigs = Canu.contigs,
    #        prefix = sample_id
    #}
    # Polishing is disabled until the blank RG issue is solved.
    #call POLISH.Dorado {
    #    input:
    #        reads = renamed_bam,
    #        draft_asm = Canu.contigs,
    #        sample_id = sample_id
    #}
    call ALN.AlignAndPlotCoverage as Reads2Asm {
        input:
            reads = fixed_reads,
            reference = Canu.contigs,
            prefix = sample_id,
            map_preset = 'map-ont'
    }
    call ALN.AlignAndPlotCoverage as Asm2Ref {
        input:
            reads = Canu.contigs,
            reference = reference_fa,
            prefix = sample_id,
            map_preset = 'asm5'
    }
    call QC.Quast {
        input:
            assembly_fa = Canu.contigs,
            reference_fa = reference_fa,
            reference_gff = reference_gff,
            fixed_reads = fixed_reads,
            Reads2AsmBam = Reads2Asm.aligned_bam,
            Reads2RefBam = Reads2RefBam,
    }

    output {
        ## canu output
        File untrimmed_contigs = Canu.contigs
        #File trimmed_contigs = CanuTrimContigs.trimmed
        ## dorado polishing output
        #File ReadsToRawAsm = Dorado.bam
        #File ReadsToRawAsmIndex = Dorado.bai
        #File PolishedContigs = Dorado.polished
        ## minimap2 output
        File Reads2Asm_bam = Reads2Asm.aligned_bam
        File Reads2Asm_bai = Reads2Asm.aligned_bai
        Array[File] Reads2Asm_coverage_plots = Reads2Asm.coverage_plots
        File Reads2Asm_plots_targz = Reads2Asm.plots_targz
        File Reads2Asm_average_depth_txt = Reads2Asm.average_depth_txt
        String Reads2Asm_average_depth = Reads2Asm.average_depth
        File Reads2Asm_average_coverage_txt = Reads2Asm.average_coverage_txt
        String Reads2Asm_average_coverage = Reads2Asm.average_coverage
        ## minimap2 output
        File Asm2Ref_bam = Asm2Ref.aligned_bam
        File Asm2Ref_bai = Asm2Ref.aligned_bai
        Array[File] Asm2Ref_plots = Asm2Ref.coverage_plots
        File Asm2Ref_plots_targz = Asm2Ref.plots_targz
        File Asm2Ref_average_depth_txt = Asm2Ref.average_depth_txt
        String Asm2Ref_average_depth = Asm2Ref.average_depth
        File Asm2Ref_average_coverage_txt = Asm2Ref.average_coverage_txt
        String Asm2Ref_average_coverage = Asm2Ref.average_coverage
        ## quast output
        File QuastReports = Quast.reports
        Int num_contigs = Quast.num_contigs
        Int asm_N50 = Quast.asm_N50
    }
}
