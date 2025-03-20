version 1.0

import "../tasks/assemble/Canu.wdl" as CANU
import "../tasks/assemble/CanuTrycycler.wdl" as CANUTry
import "../tasks/assemble/Dorado.wdl" as POLISH
import "../tasks/annotation/PlasmidCaller.wdl" as PC
import "../tasks/annotation/Bakta.wdl" as BKT
import "../tasks/align/minimap2.wdl" as MM2
# import "../tasks/annotation/LipoPredict.wdl" as LP
# import "../tasks/annotation/MLST.wdl" as MLST
#import "../tasks/annotation/ospC.wdl" as OSPC

workflow AssembleCanu {

    meta {
        description: "Denovo assembly of our reads and assembly polishing using Dorado."
    }
    parameter_meta {
        sample_id: "Sample_id to name output files"
        raw_reads_bam: "merged bam file containing reads to use in assembly polishing"
        genome_size: "[ Default: 1500000 ] length of the target genome for our samples."
        fixed_reads: "telomere fixed reads in fastq."
        reference_fa: "reference genome to pass to Quast"
        correct_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for correction."
        trim_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for trimming of corrected reads "
        assemble_error_rate: "[ Default: 0.144 ] Error rate to pass to canu for assembly of trimmed and corrected reads"
    }

    input {
        String sample_id
        Int genome_size
        File raw_reads_bam
        File fixed_reads
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
            reads = raw_reads_bam,
            draft_asm = CanuTrimContigs.trimmed,
            sample_id = sample_id
    }
    call MM2.Minimap2 {
        input:
            reads_file = fixed_fastq,
            ref_fasta = reference_fa,
            prefix = sample_id
    }
    call PC.CallPlasmids {
        input:
            contigs = Dorado.polished
    }
    call BKT.bakta {
        input:
            prefix = sample_id,
            contigs = Dorado.polished,
            rename_table = CallPlasmids.rename_table,
            bakta_db = bakta_db
    }
    call QUAST.Quast {
        input:
            assembly_fa = Bakta.fna,
            assembly_gff = Bakta.gff,
            reference_fa = reference_fa,
            reference_gff = reference_gff,
            fixed_reads = fixed_reads,
            Reads2AsmBam = Minimap2.bam,
            Reads2RefBam = Reads2RefBam
    }
    output {
        # canu output
        File untrimmed_contigs = Canu.contigs
        File trimmed_contigs = CanuTrimContigs.trimmed_contigs
        # dorado polishing output
        File ReadsToRawAsm = Dorado.bam
        File ReadsToRawAsmIndex = Dorado.bai
        File PolishedContigs = Dorado.polished
        # minimap2 output
        File ReadsToPolishedAsm = Minimap2.bam
        File ReadsToPolishedAsmIndex = Minimap2.bai
        #plasmid caller output
        File plasmids = PlasmidCaller.calls
        File renaming_table = PlasmidCaller.rename_table
        File pf32_blast_hits = PlasmidCaller.pf32_hits
        File wp_blast_hits = PlasmidCaller.wp_hits
        #bakta output
        File bakta_annotations_tsv = Bakta.annotations_tsv
        File bakta_gff3 = Bakta.gff3
        File bakta_gbff = Bakta.gbff
        File bakta_embl = Bakta.embl
        File bakta_fna = Bakta.fna
        File bakta_ffn = Bakta.ffn
        File bakta_faa = Bakta.faa
        File bakta_inference_metrics_tsv = Bakta.inference_metrics_tsv
        File bakta_hypotheticals_tsv = Bakta.hypotheticals_tsv
        File bakta_hypotheticals_faa = Bakta.hypotheticals_faa
        File bakta_summary_txt = Bakta.summary_txt
        File bakta_png = Bakta.png
        File bakta_svg = Bakta.svg
        File bakta_summary_json = Bakta.summary_json

        # quast output
        File quast_report = Quast.report
    }
}
