version 1.0

import "../../structs/Structs.wdl"

task Quast {
    meta {
        description: "Task to generate assembly QC metrics using Quast and a reference genome"
    }

    parameter_meta {
        assembly_fa: "fasta for our draft assembly"
        reference_fa: "reference genome for misassembly determination"
        reference_gff: "gff file containing annotations for the reference."
        fixed_reads: "Fastq file containing the reads for our draft assembly"
        Reads2AsmBam: "alignment of reads to the draft assembly"
        Reads2RefBam: "alignment of reads to the reference genome"
    }

    input {
        File assembly_fa
        File reference_fa
        File reference_gff
        File fixed_reads
        File Reads2AsmBam
        File Reads2RefBam
        RuntimeAttr? runtime_attr_override
    }
    Int disk_size = 50 + 7 * ceil(size(fixed_reads, "GB"))

    command <<<
        NPROC=$(awk '/^processor/{print}' /proc/cpuinfo | wc -l)
        mkdir -p quast_output
        quast.py \
            -o quast_output \
            -t "$NPROC" \
            -r "~{reference_fa}" \
            -g "~{reference_gff}" \
            --nanopore "~{fixed_reads}" \
            --ref-bam "~{Reads2RefBam}" \
            --bam "~{Reads2AsmBam}" \
            "~{assembly_fa}"
        cat quast_output/report.tsv | grep "^# contigs\t" | sed 's/# //' | awk '{print $2}' > num_contigs.txt
        cat quast_output/report.tsv | grep "N50" | awk '{print $2}' > asm_N50.txt
        tar -czf quast.tar.gz quast_output/
    >>>

    output {
        File reports = "quast.tar.gz"
        Int num_contigs = read_int("num_contigs.txt")
        Int asm_N50 = read_int("asm_N50.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       50,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "mjfos2r/quast:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}