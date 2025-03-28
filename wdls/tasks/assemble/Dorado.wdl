version 1.0

import "../../structs/Structs.wdl"

workflow Dorado {

    meta {
        description: "Align reads and polish a draft assembly using ONT Dorado."
    }
    parameter_meta {
        reads: "merged reads directly from the sequencer"
        draft_asm: "draft assembly to polish"
        prefix: "prefix to output files"
    }

    input {
        File reads
        File draft_asm
        String sample_id
    }

    call Align {
        input:
            reads = reads,
            draft_asm = draft_asm,
            sample_id = sample_id
    }

    call Polish {
        input:
            alignment = Align.bam,
            index = Align.bai,
            draft_asm = draft_asm,
            sample_id = sample_id
    }

    output {
        File bam = Align.bam
        File bai = Align.bai
        File polished = Polish.polished
    }
}

task Align {

    meta {
        description: "Align reads to a draft assembly using Dorado Align. This is required to use Dorado for polishing."
    }
    parameter_meta {
        reads:         "basecalled reads to be used with polishing"
        draft_asm:     "draft assembly to be polished"
        sample_id:     "sample_id to use for output file naming"
    }

    input {
        File reads
        File draft_asm
        String sample_id
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4 * ceil(size([reads, draft_asm], "GB"))

    command <<<
        set -euxo pipefail
        NPROC=$(cat /proc/cpuinfo | awk '/^processor/{print }' | wc -l)

        dorado aligner "~{draft_asm}" "~{reads}" | samtools sort --threads "$NPROC" > "~{sample_id}"_aligned_reads.bam
        samtools index "~{sample_id}"_aligned_reads.bam
    >>>

    output {
        File bam = "~{sample_id}_aligned_reads.bam"
        File bai = "~{sample_id}_aligned_reads.bam.bai"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:              8,
        mem_gb:                 24,
        disk_gb:                disk_size,
        boot_disk_gb:           25,
        preemptible_tries:      0,
        max_retries:            0,
        docker:                 "nanoporetech/dorado:sha3d64678bdcbb70971ee56891c01b9902eab9deea"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks:  "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries, default_attr.max_retries])
        zones:                  ["us-east1-c"]
        cpuPlatform:            "Intel Haswell"
        docker:                 select_first([runtime_attr.docker, default_attr.docker])
    }
}

task Polish {

    meta {
        description: "Polish a draft assembly using Dorado."
    }
    parameter_meta {
        alignment:        "Aligned reads output by the align task. CANNOT BE ALIGNED WITH ANY OTHER TOOL"
        index:            "Index for the aligned reads bam file."
        draft_asm:        "Draft assembly to polish"
        sample_id:        "Sample_id to use in naming output file"
    }

    input {
        File alignment
        File index
        File draft_asm
        String sample_id
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4 * ceil(size([draft_asm], "GB"))

    command <<<
        set -euxo pipefail
        NPROC=$(cat /proc/cpuinfo | awk '/^processor/{print }' | wc -l)
        dorado polish "~{alignment}" "~{draft_asm}" \
            --device cuda:0 \
            --threads "${NPROC}" \
            --bacteria > "~{sample_id}_polished.fasta"

    >>>

    output {
        File polished = "~{sample_id}_polished.fasta"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:              24,
        mem_gb:                 48,
        disk_gb:                disk_size,
        boot_disk_gb:           25,
        preemptible_tries:      0,
        max_retries:            0,
        docker:                 "nanoporetech/dorado:sha3d64678bdcbb70971ee56891c01b9902eab9deea"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks:  "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries, default_attr.max_retries])
        gpuType:                "nvidia-tesla-v100" # let's make sure this all works before we go spinning up this big money GPU...
        gpuCount:               1
        nvidiaDriverVersion:    "450.80.02"
        zones:                  ["us-central1-a"]
        cpuPlatform:            "Intel Haswell"
        docker:                 select_first([runtime_attr.docker, default_attr.docker])
    }
}
