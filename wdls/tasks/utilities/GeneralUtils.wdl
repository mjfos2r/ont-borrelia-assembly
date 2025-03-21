version 1.0
# Some tasks pulled from github.com/broadinstitute/long-read-pipelines/wdl/tasks/Utility/GeneralUtils.wdl
# Others written by Michael J. Foster (github.com/mjfos2r)

import "../../structs/Structs.wdl"

workflow ValidateAndDecompressSeqRun {
    meta {
        description: "This is a workflow to validate and decompress a sequencing run."
    }

    parameter_meta {
        run_tarball: "description of input"
        md5sum: "file containing the md5sum for the tarball"
    }

    input {
        File run_tarball
        File md5sum
    }

    call ValidateMd5sum { input: file = run_tarball, checksum = md5sum }

    if (ValidateMd5sum.is_valid) { call DecompressRunTarball { input: tarball = run_tarball } }

    output {
        # Instead of maps, output parallel arrays
        Array[String]? barcodes = DecompressRunTarball.barcodes
        Array[File]? bam_lists = DecompressRunTarball.bam_lists
        Array[Int]? bam_counts = DecompressRunTarball.bam_counts
    }
}

task CompressTarPigz {
    meta {
        description: "compress files into a tarball using parallel gzip and generate an md5 checksum"
    }

    parameter_meta {
        files: "List of files to zip up."
        name: "Name of the tar.gz file."
    }

    input {
        Array[File] files
        String name

        Int num_cpus = 16
        Int mem_gb = 32
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(files, "GB")) + 15

    command <<<
        NUM_CPUS=$(cat /proc/cpuinfo | awk '/^processor/{print }' | wc -l)
        set -euxo pipefail # crash out
        mkdir -p tarpit/
        for ff in ~{sep=' ' files}; do cp "${ff}" tarpit/ ; done
        tar -cf - -C tarpit/ . | pigz -p "${NUM_CPUS}" > ~{name}.tar.gz
        md5sum ~{name}.tar.gz >~{name}.tar.gz.md5
    >>>

    output {
        File tarball = "~{name}.tar.gz"
        File checksum = "~{name}.tar.gz.md5"
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_gb,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "mjfos2r/basic:latest"
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

task GetTodayDate {
    meta {
        description: "Generates a YYYY-MM-DD date of today (when this task is called). UTC."
        volatile: true
    }

    input {
        Int num_cpus = 1
        Int mem_gb = 4
        Int disk_size = 10
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        date '+%Y-%m-%d'
    >>>

    output {
        String yyyy_mm_dd = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_gb,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "mjfos2r/basic:latest"
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

task ValidateMd5sum {
    meta {
        description: "simple task to validate checksum of a file"
    }

    parameter_meta {
        file: "input_file to validate"
        checksum: "file containing checksum"
    }

    input {
        File file
        File checksum

        RuntimeAttr? runtime_attr_override
    }

    # Estimate disk size - the compressed archive plus space for extraction and output
    Int disk_size = ceil(size(file, "GB")) + 15

    command <<<
    set -euxo pipefail # crash out

    EXPECTED_MD5=$(awk '{print $1}' ~{checksum})
    ACTUAL_MD5=$(md5sum ~{file} | awk '{print $1}')
    if [ "$EXPECTED_MD5" != "$ACTUAL_MD5" ]; then
        echo "ERROR: CHECKSUM VALIDATION FAILED FOR ~{file}"
        echo "ERROR: Expected: $EXPECTED_MD5"
        echo "ERROR:   Actual: $ACTUAL_MD5"
        echo "false" > valid.txt
    else
        echo "Checksum validation successful for ~{file}"
        echo "true" > valid.txt
    fi
    >>>

    output {
        Boolean is_valid = read_boolean("valid.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "mjfos2r/basic:latest"
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

task DecompressRunTarball {
    meta {
        description: "Decompress a validated run tarball using pigz"
    }

    parameter_meta {
        tarball: "validated tarball to decompress"
    }

    input {
        File tarball

        # Runtime parameters
        Int num_cpus = 4
        Int mem_gb = 8
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(size(tarball, "GB") * 3)

    command <<<
        set -euxo pipefail

        mkdir -p extracted

        tar --use-compress-program=pigz -xf ~{tarball} -C extracted
        # Get a list of our directories, pull the barcode ID, and make a list of files for each
        find extracted -mindepth 1 -maxdepth 1 -type d | sort > directory_list.txt
        cut -d'/' -f2 directory_list.txt > barcodes.txt
        mkdir -p file_lists

        # Create/clear the counts file before the loop
        true > bam_counts.txt

        while read -r dir_path; do
            BARCODE=$(basename "$dir_path")
            find "$dir_path" -name "*.bam" | sort > "file_lists/${BARCODE}_bams.txt"
            # Count the files and append to the counts file
            wc -l < "file_lists/${BARCODE}_bams.txt" >> bam_counts.txt
        done < directory_list.txt

        # Count the number of directories for verification
        wc -l directory_list.txt | awk '{print $1}' > directory_count.txt

        # Use python to dump a JSON so we can immediately parse it into a damn map.
        python3 - <<EOF
        import os
        import json

        bc_to_bams = {}

        with open('barcodes.txt', 'r') as f:
            barcodes = [ line.strip() for line in f ]

        for barcode in barcodes:
            bam_file_list = f"file_lists/{barcode}_bams.txt"
            with open(bam_file_list, 'r') as f:
                bam_files = [line.strip() for line in f]
            barcode_to_bams[barcode] = bam_files

        # write to file and dump to stdout so we can create a map
        # and use that for downstream scatter logic
        with open('barcode_to_bams.json', 'w') as f:
            json.dump(barcode_to_bams, f, indent=2)

        print(json.dumps(barcode_to_bams, indent=2))
        EOF
    >>>

    output {
        Int directory_count = read_int("directory_count.txt")
        Array[Int] bam_counts = read_lines("bam_counts.txt")
        Array[String] barcode_dirs = read_lines("directory_list.txt")
        Array[String] barcodes = read_lines("barcodes.txt")
        Array[File] bam_lists = glob("file_lists/*_bams.txt")
        File barcode_to_bams_json = "barcode_to_bams.json"
        Map[String, Array[File]] barcode_to_bams_map = read_json(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_gb,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "mjfos2r/basic-python:3.11-slim"
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

task CreateMap {
    # gfys WDL1.0
    input {
        Array[String] keys
        Array[File] values
    }

    command <<<
        python3 <<CODE
        import json
        keys = ["~{sep='","' keys}"]
        values = ["~{sep='","' values}"]
        result = {}
        for i in range(len(keys)):
            result[keys[i]] = values[i]
        with open("map.json", "w") as f:
            json.dump(result, f)
        CODE
    >>>

    output {
        Map[String, File] mmap = read_json("map.json")
    }

    runtime {
        docker: "python:3.11-slim"
    }
}