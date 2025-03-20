version 1.0
import "../../structs/Structs.wdl"

workflow Canu {

    meta {
        description: "A workflow that runs the Canu 3-step assembly (correct, trim, assemble). Tested on a small genome (malaria ~23mb), larger genomes may require some changes including tweaks to the default resource allocation. Currently assumes nanopore reads"
        # Adapted from broadinstitute/long-read-pipeline and expanded for use in MJF's Borrelia pipeline(s)
    }
    parameter_meta {
        reads:                "reads to be canu-corrected"
        genome_size:          "estimate on genome size (parameter to canu's 'genomeSize')"
        correct_error_rate:   "parameter to canu's 'correctedErrorRate'"
        trim_error_rate:      "parameter to canu's 'correctedErrorRate'"
        assemble_error_rate:  "parameter to canu's 'correctedErrorRate'"
        prefix:               "[ Default: 'canu' ] prefix to output files"
    }

    input {
        File reads

        Int genome_size
        Float correct_error_rate
        Float trim_error_rate
        Float assemble_error_rate

        String prefix
    }

    call Correct {
        input:
            reads = reads,
            genome_size = genome_size,
            error_rate = correct_error_rate,
            prefix = prefix
    }

    call Trim {
        input:
            genome_size = genome_size,
            corrected_reads = Correct.corrected_reads,
            error_rate = trim_error_rate,
            prefix = prefix
    }

    call Assemble {
        input:
            genome_size = genome_size,
            trimmed_reads = Trim.trimmed_reads,
            error_rate = assemble_error_rate,
            prefix = prefix
    }

    output {
        File contigs = Assemble.canu_contigs_fasta
    }
}

# performs canu correct on raw reads, currently assumes ONT reads
task Correct {
    input {
        File reads
        Int genome_size
        Float error_rate
        String prefix = "canu"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:        "reads to be canu-corrected"
        genome_size:  "estimate on genome size (parameter to canu's 'genomeSize')"
        error_rate:   "parameter to canu's 'correctedErrorRate'"
        prefix:       "[ Default: 'canu' ] prefix to output files"

    }

    Int disk_size = 150 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        canu -correct \
            -p "~{prefix}" -d canu_correct_output \
            genomeSize="~{genome_size}"m \
            corMaxEvidenceErate=0.15 \
            correctedErrorRate="~{error_rate}" \
            -nanopore \
            "~{reads}"
    >>>

    output {
        File corrected_reads = "canu_correct_output/~{prefix}.correctedReads.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
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

# performs canu trim on corrected reads
task Trim {
    input {
        File corrected_reads
        Int genome_size
        Float error_rate
        String prefix = "canu"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        corrected_reads:   "reads that have been canu-corrected"
        genome_size:       "estimate on genome size (parameter to canu's 'genomeSize')"
        error_rate:        "parameter to canu's 'correctedErrorRate'"
        prefix:            "[ Default: 'canu' ] prefix to output files"

    }

    Int disk_size = 50 * ceil(size(corrected_reads, "GB"))

    command <<<
        set -euxo pipefail

        canu -trim \
            -p "~{prefix}" -d canu_trim_output \
            genomeSize="~{genome_size}"m \
            correctedErrorRate="~{error_rate}" \
            -nanopore-corrected \
            "~{corrected_reads}"
    >>>

    output {
        File trimmed_reads = "canu_trim_output/~{prefix}.trimmedReads.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
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

# performs assembly on corrected, then trimmmed reads
task Assemble {
    input {
        Int genome_size
        File trimmed_reads
        Float error_rate
        String prefix = "canu"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        trimmed_reads:  "reads that have been canu-corrected-trimmed"
        genome_size:    "estimate on genome size (parameter to canu's 'genomeSize')"
        error_rate:     "parameter to canu's 'correctedErrorRate'"
        prefix:         "[ Default: 'canu' ] prefix to output files"

    }

    Int disk_size = 50 * ceil(size(trimmed_reads, "GB"))

    command <<<
        set -euxo pipefail

        canu -assemble \
            -p "~{prefix}" -d canu_assemble_output \
            genomeSize="~{genome_size}"m \
            correctedErrorRate="~{error_rate}" \
            -nanopore-corrected \
            "~{trimmed_reads}"
    >>>

    output {
        File canu_contigs_fasta = "canu_assemble_output/~{prefix}.contigs.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-canu:0.1.0"
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

workflow CanuTrycycler {

    meta {
        description: "A workflow that runs the Canu 3-step assembly (correct, trim, assemble). Tested on a small genome (malaria ~23mb), larger genomes may require some changes including tweaks to the default resource allocation. Currently assumes nanopore reads"
        # Adapted from broadinstitute/long-read-pipeline and expanded for use in MJF's Borrelia pipeline(s)
    }
    parameter_meta {
        reads:                "reads to be canu-corrected"
        genome_size:          "estimate on genome size (parameter to canu's 'genomeSize')"
        prefix:               "[ Default: 'canu' ] prefix to output files"
    }

    input {
        File reads
        Int genome_size = 1500000
        Int subsample_id
        String prefix = "canu"
    }

    call CanuFast { input: reads=reads, genome_size=genome_size, prefix=prefix }
    call CanuTrimContigs { input: contigs=CanuFast.untrimmed_contigs, subsample_id=subsample_id }

    output {
        File trimmed_contigs = CanuTrimContigs.trimmed_contigs
    }
}

# fast mode for trycycler pipeline where it's assembling read subsets.
task CanuFast {
    input {
        File reads
        Int genome_size = 1500000
        String prefix = "canu"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:          "reads that have been canu-corrected-trimmed"
        genome_size:    "[ Default:  1500000  ] Estimate on genome size (parameter to canu's 'genomeSize')."
        prefix:         "[ Default: 'canu' ] Prefix to use for output file"
    }

    Int disk_size = 150 * ceil(size(trimmed_reads, "GB"))

    command <<<
        set -euxo pipefail
        NUM_CPUS=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        canu -fast \
            -p canu \
            -d canu_temp \
            genomeSize="$genome_size" \
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
        docker:             "us.gcr.io/terra-942df462/canu:latest" ## TODO: Canu Container!
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
        subset_id:   "Indicate which read subset these contigs originated from in the output file."
    }

    Int disk_size = 50 * ceil(size(trimmed_reads, "GB"))
    String output_file = if defined(subsample_id) then "assembly_~{subset_id}.fasta"
                            else "~{prefix}.contigs.trimmed.fasta"
    command <<<
        set -euxo pipefail
        mkdir -p trimmed_asm
        /opt/trycycler/scripts/canu_trim.py "~{contigs}" > trimmed_asm/"~{output_file}"
    >>>

    output {
        File trimmed = glob("trimmed_asm/*.fasta")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/terra-942df462/trycycler:latest" ##TODO: WRITE MY OWN CONTAINER
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