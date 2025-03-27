version 1.0

import "../workflows/QC_single.wdl" as QC
import "../tasks/utilities/GeneralUtils.wdl" as UTILS
import "../tasks/utilities/BamUtils.wdl" as BAM
import "../tasks/utilities/Chopper.wdl" as CHP
import "../tasks/telomeres/ResolveTelomeres.wdl" as TELO

workflow ONT_ReadPrepAndQC {
    meta {
        description: "Begin processing the samples given the DataTable in Terra using outputs from PreprocesingAndRunQC"
    }

    parameter_meta {
        sample_id: "sample_alias for the sample to be processed"
        raw_bams: "array of unmerged bam files to merge and process"
        reference_fa: "reference genome to align reads against"
        contam_fa: "file containing barcode, DCS_fa, and adapter sequences to filter and trim from reads"
    }

    input {
        String sample_id
        File merged_bam
        File reference_fa
        File contam_fa
    }
    # firstly, rename our merged_bam to the sample_id
    call UTILS.RenameFile { input: file = merged_bam, new_name = sample_id}
    # make it a fastq
    call BAM.Bam2Fastq { input: input_bam = RenameFile.renamed_file }
    # Get to the choppah (terminate the bad reads)
    call CHP.Chopper { input: input_reads = Bam2Fastq.fastq, contam_fa = contam_fa }
    # run my bespoke in-silico reverse-telomerase task
    call TELO.ResolveTelomeres { input: reads = Chopper.clean_fq, sample_id = sample_id }

    # now generate our QC for these reads.
    call QC.SingleReadsQC {
        input:
            sample_id = sample_id,
            reference_fa = reference_fa,
            clean_fastq = Chopper.clean_fq,
            fixed_fastq = ResolveTelomeres.fixed_reads,
    }

    output {
        File renamed_bam = RenameFile.renamed_file
        File cleaned_reads = Chopper.clean_fq
        File telomere_bed = ResolveTelomeres.telomere_bed
        File telomere_read_ids = ResolveTelomeres.telomere_read_ids
        File raw_telomere_reads = ResolveTelomeres.telomere_fastq
        File clipped_telomere_reads = ResolveTelomeres.clipped_telomeres
        File fixed_reads = ResolveTelomeres.fixed_reads
        File cleaned_Reads2Ref = SingleReadsQC.raw_telo_bam
        Array[File] raw_telo_coverage_plots = SingleReadsQC.raw_telo_coverage_plots
        File fixed_Reads2Ref = SingleReadsQC.fixed_telo_bam
        Array[File] fixed_telo_coverage_plots = SingleReadsQC.fixed_telo_coverage_plots
    }
}
