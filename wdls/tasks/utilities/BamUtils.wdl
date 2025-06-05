version 1.0
# Pulled from github.com/broadinstitute/long-read-pipelines/wdl/tasks/Utility/GeneralUtils.wdl
# With some changes

import "../../structs/Structs.wdl"

task MergeBams {
    meta {
        desciption: "when provided with a directory of bam files for a given barcode, merge them all into a single file, sort it, and index it. also return flagstats"
    }

    parameter_meta {
        input_bams: "list of bams to merge"
        name: "name for merged bam. do not specify .bam"
    }

    input {
        Array[File] input_bams
        String name # should be something like Barcode##
        RuntimeAttr? runtime_attr_override
    }

    String output_bam ="~{name}.bam"

    Int disk_size = 50 + 2*ceil(size(input_bams, "GB"))

    command <<<
    set -euxo pipefail # if anything breaks crash out

    # get the number of procs we have available
    NPROCS=$( grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}' )

    # list our bams that we're merging
    echo "[INFO] :: Merging BAMs ::"
    for bam in "~{sep=' ' input_bams}"; do
        echo "  - $bam"
    done

    # merge and sort em
    samtools merge \
        -f \
        -@ "$NPROCS" \
        -o merged.tmp.bam \
        "~{sep=' ' input_bams}"

    samtools sort \
        -@ "$NPROCS" \
        -m $(( 8 / "$NPROCS"))G \
        -o "~{output_bam}" \
        merged.tmp.bam

    # get the index
    samtools index "~{output_bam}"

    # and get some stats
    samtools flagstat "~{output_bam}"

    # cleanup
    rm merged.tmp.bam
    >>>

    output {
        File merged_bam = "~{output_bam}"
        File merged_bam_index = "~{output_bam}.bai"
        File stats = "~{output_bam}.stats"
    }
    # no preempt.
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "mjfos2r/samtools:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Bam2Fastq {
        meta {
        desciption: "convert bam to fastq.gz file and preserve all tags written by Dorado basecaller."
    }

    parameter_meta {
        input_bam: "Bam file to convert to fastq"
    }

    input {
        File input_bam
        String? sample_id
        Int num_cpus = 8
        Int mem_gb = 32
        RuntimeAttr? runtime_attr_override
    }
    String bam_name = basename(input_bam)
    String filename = select_first([sample_id, bam_name])
    Int disk_size = 365 + 3*ceil(size(input_bam, "GB"))

    command <<<
    set -euo pipefail # if anything breaks crash out

    # get the number of procs we have available
    NPROCS=$( cat /proc/cpuinfo | grep '^processor' | tail -n1 | awk '{print $NF+1}' )

    samtools sort -n ~{input_bam} > sorted.bam
    # preserve all tags that dorado puts in the BAM.
    # and add _R1 to keep things consistent with the thiagen tool.
    samtools fastq -@ "$NPROCS" -T '*' -0 "~{filename}_R1.fastq.gz" sorted.bam
    >>>

    output {
        File fastq = "~{filename}.fastq.gz"
    }
    # no preempt.
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_gb,
        disk_gb:            disk_size,
        boot_disk_gb:       50,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "mjfos2r/samtools:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
