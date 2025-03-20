version 1.0

import "../tasks/assemble/Canu.wdl" as CANU
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
        fixed_reads: "telomere fixed reads in fastq."
        polished_contigs: "polished contigs output by dorado"
        bakta_db "tar.gz containing the full bakta database"
        fixed_telo_bam: "Bam file for the aligned reads against the reference geome"
        reference_fa: "Fasta for our reference genome to pass to Quast"
        reference_gff "GFF3 for our reference genome to pass to Quast"
    }

    input {
        String sample_id
        File fixed_reads
        File polished_contigs
        File bakta_db
        File fixed_telo_bam
        File reference_fa
        File reference_gff

    }
    call MM2.Minimap2 {
        input:
            reads_file = fixed_fastq,
            ref_fasta = polished_contigs,
            prefix = sample_id
    }
    call PC.CallPlasmids {
        input:
            contigs = Dorado.polished
    }
    call BKT.bakta {
        input:
            prefix = sample_id,
            contigs = polished_contigs
            rename_table = CallPlasmids.rename_table
            bakta_db = bakta_db
    }
    call QUAST.Quast {
        input:
            assembly_fa = Bakta.fna,
            assembly_gff = Bakta.gff,
            reference_fa = reference_fa,
            reference_gff = reference_gff,
            fixed_reads = fixed_reads,
            Reads2AsmBam = Minimap2.bam
            Reads2RefBam = fixed_telo_bam
    }
    output {
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
