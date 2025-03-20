version 1.0
import "../../structs/Structs.wdl"

task CreateBarcodeMap {
    # this is super scuffed. might be entirely pointless given the parsed samplesheet output below.
    input {
        File sample_sheet
    }

    command <<<
        # Extract barcode and sample columns from TSV and convert to key=value format
        awk -F'\t' 'NR>1 {print $1 "\t" $2}' ~{sample_sheet} > barcode_map.txt
    >>>

    output {
        Map[String, String] barcode_to_sample_map = read_map("barcode_map.txt")
    }

    runtime {
        docker: "MY_BASE_DOCKER_IMAGE"
        disks: "local-disk 10 HDD"
    }
}

task RenameFile {
    # This is broken and needs rewriting.
    input {
        File original_file
        Map[String, String] barcode_map
    }

    # Extract the barcode from the filename
    String filename = basename(original_file)

    command <<<
        set -e

        # Get original filename
        FILENAME="~{filename}"

        # Loop through all possible barcodes in the map
        # This approach handles the mapping in bash since WDL can't directly index maps with variables
        FOUND=false
        NEW_NAME=""

        # Create a temporary file with the mapping via a cursed herefile.
        cat <<EOF > temp_map.txt
        ~{write_map(barcode_map)}
        EOF

        # For each barcode in our mapping
        while IFS=$'\t' read -r barcode sample; do
            # If the filename contains this barcode
            if [[ "$FILENAME" == *"$barcode"* ]]; then
                # Replace barcode with sample name
                NEW_NAME="${FILENAME/$barcode/$sample}"
                FOUND=true
                break
            fi
        done < temp_map.txt

        # If no match found, keep original name
        if [ "$FOUND" = false ]; then
            NEW_NAME="$FILENAME"
            echo "Warning: No matching barcode found for $FILENAME"
        fi

        # Copy the file with the new name
        cp "~{original_file}" "$NEW_NAME"
    >>>

    output {
        File renamed_file = glob("*")[0]  # This assumes only one output file
    }

    runtime {
        docker: "MY_BASE_DOCKER_IMAGE"
        disks: "local-disk 10 HDD"
    }
}


task ParseSamplesheet {
    meta {
        description: "Parse our run's samplesheet so that we can rename our files to the sample_alias"
    }

    parameter_meta {
        samplesheet: "csv samplesheet"
    }

    input {
        File samplesheet
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        python3 << EOF
        from collections import defaultdict

        with open("~{samplesheet}", 'r') as infile:
            lines = infile.readlines()

        keys = lines[0].strip().split(',')
        vals = defaultdict(list)

        for row in lines[1:]:
            line = row.strip().split(',')
            for idx, key in enumerate(keys):
                if idx < len(line):  # no bad rows...
                    vals[key].append(line[idx])

        # Also write individual files for debugging and use later.
        for key in vals.keys():
            cleaned_key = key.lower().replace(' ', '_')
            with open(f"{cleaned_key}.txt", 'w') as out:
                out.write('\n'.join(vals[key]))
        EOF
    >>>

    output {
        # Optional individual outputs for backward compatibility
        Array[String] flow_cell_id = read_lines("flow_cell_id.txt")
        Array[String] position_id = read_lines("position_id.txt")
        Array[String] experiment_id = read_lines("experiment_id.txt")
        Array[String] flow_cell_product_code = read_lines("flow_cell_product_code.txt")
        Array[String] kit = read_lines("kit.txt")
        Array[String] barcodes = read_lines("barcode.txt")
        Array[String] sample_id = read_lines("alias.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            4,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "python:3.9-slim"
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