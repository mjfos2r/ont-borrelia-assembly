version 1.0
import "../tasks/utilities/BamUtils.wdl" as BAM

workflow BAMtoFastq {
    meta {
        description: "Simple workflow to convert BAMs to fastqs."
    }
    parameter_meta {
        sample_id: "sample_id for the assembly we're classifying"
        input_bam: "bam file to be converted"
    }
    input {
        String sample_id
        File input_bam
    }
    call BAM.Bam2Fastq {
        input:
        input_bam = input_bam,
        sample_id = sample_id
    }
    output {
        File reads_fastq = Bam2Fastq.fastq
    }
}