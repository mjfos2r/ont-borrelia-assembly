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


task ParseSamplesheetToDataTable {
    meta {
        description: "Parse the samplesheet and extract per-column arrays, also create a full TSV with an added raw_bams column for easy use as a Cromwell DataTable."
    }

    parameter_meta {
        samplesheet: "CSV-formatted samplesheet with sample metadata. Formatted per ONT's guidelines except renamed column: alias => sample_id and output a DataTable for Cromwell input."
        raw_bam_paths: "Array of BAM file paths per sample, each being an array of BAMs."
    }

    input {
        File samplesheet
        Array[Array[File]] raw_bam_paths

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        python3 << EOF
        import csv
        import json
        from collections import defaultdict

        # Get our bam paths and format this as a json dump.
        bam_paths = ~{sep=',' raw_bam_paths}
        bam_path_strings = [json.dumps([str(p) for p in paths]) for paths in bam_paths]

        # Read in our samplesheet CSV
        with open("~{samplesheet}", 'r') as infile:
            reader = list(csv.reader(infile))
            header = reader[0]
            rows = reader[1:]

        # Transpose into column-wise defaultdict(list)
        vals = defaultdict(list)
        for row in rows:
            for idx, key in enumerate(header):
                if idx < len(row):
                    vals[key].append(row[idx])

        # dump each col to a file. (which can be input as Array[String])
        for key, values in vals.items():
            outname = key.lower().replace(" ", "_") + ".txt"
            with open(outname, 'w') as out:
                out.write('\n'.join(values))

        # Write bam_paths.txt (for downstream input as Array[File])
        with open("bam_paths.txt", 'w') as out:
            out.write('\n'.join(bam_path_strings))

        # Get the experiment id and set the output filename for clarity once downloaded.
        DataTable_out = f"{vals['experiment_id'][1]}__DataTable.tsv"

        # Write combined TSV with bam_paths as new column that can be uploaded and directly used as a DataTable in Cromwell.
        new_header = header + ["raw_bams"]
        with open(f"{DataTable_out}", 'w', newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(new_header)
            for i, row in enumerate(rows):
                bp = bam_path_strings[i] if i < len(bam_path_strings) else "[]" # this is bad but oh well.
                writer.writerow(row + [bp])
        EOF
    >>>

    output {
        Array[String] flow_cell_id = read_lines("flow_cell_id.txt")
        Array[String] position_id = read_lines("position_id.txt")
        Array[String] experiment_id = read_lines("experiment_id.txt")
        Array[String] flow_cell_product_code = read_lines("flow_cell_product_code.txt")
        Array[String] kit = read_lines("kit.txt")
        Array[String] barcode = read_lines("barcode.txt")
        Array[String] sample_id = read_lines("alias.txt")
        Array[String] bam_paths = read_lines("bam_paths.txt")
        File samplesheet_with_bams = glob("*__DataTable.tsv")[0]
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
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

