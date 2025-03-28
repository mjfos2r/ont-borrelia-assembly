version 1.0

import "../tasks/assemble/Canu.wdl" as CANU
import "../tasks/assemble/CanuTrycycler.wdl" as CANUTry
import "../tasks/assemble/Dorado.wdl" as POLISH
import "../tasks/align/minimap2.wdl" as MM2

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
        correct_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for correction."
        trim_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for trimming of corrected reads "
        assemble_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for assembly of trimmed and corrected reads"
    }

    input {
        String sample_id
        File renamed_bam
        Float genome_size
        File fixed_reads
        File fixed_Reads2Ref
        File reference_fa
        Float correct_error_rate = 0.144
        Float trim_error_rate = 0.144
        Float assemble_error_rate = 0.144
    }

    call CANU.Canu {
        input:
            reads = fixed_reads,
            genome_size = genome_size,
            correct_error_rate = correct_error_rate,
            trim_error_rate = trim_error_rate,
            assemble_error_rate = assemble_error_rate,
            prefix = sample_id,
    }
    call CANUTry.CanuTrimContigs {
        input:
            contigs = Canu.contigs,
            prefix = sample_id
    }
    call POLISH.Dorado {
        input:
            reads = fixed_Reads2Ref,
            draft_asm = CanuTrimContigs.trimmed,
            sample_id = sample_id
    }
    call MM2.Minimap2 as Reads2Asm {
        input:
            reads_file = fixed_reads,
            ref_fasta = reference_fa,
            prefix = sample_id
    }

    output {
        # canu output
        File untrimmed_contigs = Canu.contigs
        File trimmed_contigs = CanuTrimContigs.trimmed
        # dorado polishing output
        File ReadsToRawAsm = Dorado.bam
        File ReadsToRawAsmIndex = Dorado.bai
        File PolishedContigs = Dorado.polished
        # minimap2 output
        File ReadsToPolishedAsm = Reads2Asm.aligned_bam
        File ReadsToPolishedAsmIndex = Reads2Asm.aligned_bai
    }
}
