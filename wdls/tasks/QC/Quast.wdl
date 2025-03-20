task Quast {
    meta {
        description: "Task to generate assembly QC metrics using Quast and a reference genome"
    }

    parameter_meta {
        assembly_fa: "fasta for our draft assembly"
        assembly_gff: "gff3 file for our draft assembly"
        reference_fa: "reference genome for misassembly determination"
        reference_gff: "gff file containing annotations for the reference."
        fixed_reads: "Fastq file containing the reads for our draft assembly"
        Reads2AsmBam: "alignment of reads to the draft assembly"
        Reads2RefBam: "alignment of reads to the reference genome"
    }

    input {
        File assembly_fa
        File assembly_gff
        File reference_fa
        File reference_gff
        File fixed_reads
        File Reads2AsmBam
        File Reads2RefBam
    }

    command <<<
        NPROC=$(awk '/^processor/{print}' /proc/cpuinfo | wc -l)
        mkdir -p quast_output
        quast.py \
            -o quast_output \
            -t "$NPROC" \
            -r "~{reference_fa}" \
            -g "~{reference_gff}" \
            --ref-bam "~{Reads2RefBam}"
            --operons "~{assembly_gff}" \
            --upper-bound-assembly "~{fixed_reads}" \
            --nanopore "~{fixed_reads}"
            --bam "~{Reads2AsmBam}"
        tar -czf quast.tar.gz quast_output/
    >>>

    output {
        File report = "quast.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "mjfos2r/Quast:latest" ##TODO: WRITE MY OWN CONTAINER
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