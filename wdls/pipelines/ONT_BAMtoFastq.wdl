version 1.0
import "../tasks/utilities/BamUtils.wdl" as BAM

workflow BAMtoFastq {
    meta {
        description: "Simple workflow to convert BAMs to fastqs."
        author: "Michael J. Foster"
    }
    parameter_meta {
        sample_id: "sample_id for the bam file we're converting. Not required, defaults to file basename."
        input_bam: "bam file to be converted"
        st_params: "parameters for 'samtools fastq' command. [default: -T '*']"
    }
    input {
        String? sample_id
        String? postfix
        File input_bam
        String? st_params
    }
    call BAM.Bam2Fastq {
        input:
        input_bam = input_bam,
        sample_id = sample_id,
        st_params = st_params,
        postfix = postfix
    }
    output {
        File reads_fastq = Bam2Fastq.fastq
    }
}