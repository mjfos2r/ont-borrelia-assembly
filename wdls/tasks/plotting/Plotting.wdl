version 1.0
import "../../structs/Structs.wdl"

task PlotBamCoverage {
    input {
        File input_bam
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        input_bam: "input alignment to plot coverage of"
    }

    Int disk_size = 50 + ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail
        bam2plot from_bam \
            --bam "~{input_bam}" \
            --outpath "BamCovPlots" \
            --rolling_window 50 \
            --threshold 10 -s -c \
            -n 22
        tar -zcvf BamCovPlots.tar.gz BamCovPlots/

        # Get our total genome size and then get the average coverage across all positions.
        GSIZE=$(samtools view -H "~{input_bam}" | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
        samtools depth "~{input_bam}" | awk -v gsize="$GSIZE" '{sum+=$3} END { print sum/gsize "X"}' > average_depth.txt
        samtools coverage "~{input_bam}" | awk 'NR>1 {sum+=$6*$3; total+=$3} END {print sum/total "%"}' > average_coverage.txt
    >>>

    output {
        Array[File] plots = glob("BamCovPlots/*")
        File plots_targz = "BamCovPlots.tar.gz"
        File average_depth_txt  = "average_depth.txt"
        String average_depth = read_string("average_depth.txt")
        File average_coverage_txt  = "average_coverage.txt"
        String average_coverage = read_string("average_coverage.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "mjfos2r/plotting:latest"
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