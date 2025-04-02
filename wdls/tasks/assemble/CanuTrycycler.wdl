version 1.0

import "../../structs/Structs.wdl"

# fast mode for trycycler pipeline where it's assembling read subsets.
task CanuFast {
    input {
        File reads
        Float genome_size = 1.5
        String prefix = "canu"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:          "reads that have been canu-corrected-trimmed"
        genome_size:    "[ Default:  1500000  ] Estimate on genome size (parameter to canu's 'genomeSize')."
        prefix:         "[ Default: 'canu' ] Prefix to use for output file"
    }

    Int disk_size = 150 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail
        NUM_CPUS=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        canu -fast \
            -p canu \
            -d canu_temp \
            genomeSize="~{genome_size}m" \
            useGrid=false \
            -nanopore "~{reads}"
    >>>

    output {
        File untrimmed_contigs = "canu_temp/~{prefix}.contigs.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "mjfos2r/canu:latest"
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


task CanuTrimContigs {
    meta {
        description: "task to run trycycler's 'canu_trim.py' script for contig circularization"
    }

    input {
        File contigs
        String? subsample_id
        String? prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        contigs:  "Assembled contigs to be trimmed and circularized"
        subsample_id:   "Indicate which read subset these contigs originated from in the output file."
        prefix:   "Indicate prefix for output file naming."
    }

    Int disk_size = 50 * ceil(size(contigs, "GB"))
    String output_file = if defined(subsample_id) then "assembly_~{subsample_id}.fasta"
                            else "~{prefix}.contigs.trimmed.fasta"
    command <<<
        set -euxo pipefail
        mkdir -p trimmed_asm
        /opt/trycycler/scripts/canu_trim.py "~{contigs}" > trimmed_asm/"~{output_file}"
    >>>

    output {
        File trimmed = output_file # feels broken?
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "mjfos2r/trycycler:latest"
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
