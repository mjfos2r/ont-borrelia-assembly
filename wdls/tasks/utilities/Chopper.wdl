version 1.0

import "../../structs/Structs.wdl"

task Chopper {
    input {
        File input_reads
        File contam_fa
        Int min_quality = 10
        Int min_length = 500
        Boolean compress_output = true
        String? extra_args

        Int num_cpus = 4
        RuntimeAttr? runtime_attr_override
    }

    String basename = sub(basename(input_reads), "\\..*$", "")
    Boolean is_input_gzipped = sub(input_reads, ".*\\.", "") == "gz"
    String output_filename = basename + ".clean.fq" + (if compress_output then ".gz" else "")

    Int disk_size = 3*ceil(size(input_reads, "GB")) + 10

    command <<<
        set -euxo pipefail

        NPROCS=$( grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}' )

        # Determine input type and run appropriate Chopper command
        if [ "~{is_input_gzipped}" == "true" ]; then
            # Input is gzipped
            gunzip -c ~{input_reads} | \
            chopper --contam ~{contam_fa} -q ~{min_quality} -l ~{min_length} ~{extra_args} --threads "$NPROCS" | \
            ~{if compress_output then "gzip" else "cat"} > ~{output_filename}
        else
            # Input is not gzipped
            chopper --contam ~{contam_fa} -q ~{min_quality} -l ~{min_length} -i ~{input_reads} ~{extra_args} --threads "$NPROCS" | \
            ~{if compress_output then "gzip" else "cat"} > ~{output_filename}
        fi

        # Calculate read stats before and after filtering
        if [ "~{is_input_gzipped}" == "true" ]; then
            echo "Input reads: $(gunzip -c ~{input_reads} | grep -c "^@")" > stats.txt
        else
            echo "Input reads: $(grep -c "^@" ~{input_reads})" > stats.txt
        fi

        if [ "~{compress_output}" == "true" ]; then
            echo "Output reads: $(gunzip -c ~{output_filename} | grep -c "^@")" >> stats.txt
        else
            echo "Output reads: $(grep -c "^@" ~{output_filename})" >> stats.txt
        fi

        # Convert to map format
        sed 's/: /\t/' stats.txt > stats_map.txt
    >>>

    output {
        File clean_fq = "~{output_filename}"
        File stats = "chopper_stats.tsv"
        Map[String, Int] stats_map = read_map("chopper_stats.tsv")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/terra-942df462/chopper:latest" #FIX
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