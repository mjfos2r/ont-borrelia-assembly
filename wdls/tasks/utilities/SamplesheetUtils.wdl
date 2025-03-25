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
        # I know this is busted but we've gotta see how this changes things.
    }

    parameter_meta {
        samplesheet: "CSV-formatted samplesheet with sample metadata. Formatted per ONT's guidelines except renamed column: alias => sample_id and output a DataTable for Cromwell input."
        merged_bams: "Array[File]: List of bam file paths (delocalized) for each barcode."
    }

    input {
        File samplesheet
        Array[File] merged_bams

        Int num_cpu = 4
        Int mem_gb = 16
        RuntimeAttr? runtime_attr_override
    }
    Int disk_size = 50 + ceil(size(merged_bams, "GB"))
    command <<<
        set -euxo pipefail

        cat ~{sep='\n' merged_bams} > merged_bams.txt

        python3 - <<EOF
        import os
        import csv
        import json

        from collections import defaultdict

        # if this actually outputs delocalized paths I'll be so surprised.
        # let's try it anyway!
        barcode_to_bam = {}
        with open("merged_bams.txt", 'r') as f:
            for line in f:
                path = line.strip()
                print(f"Input Path From Cromwell: {path}")
                #gcs_path = path.replace("cromwell_root/", "gs://")
                filename = os.path.basename(path)
                barcode = filename.split(".")[0] # barcode01.merged.bam
                barcode_to_bam[barcode] = path

        experiment_id = ""
        rows = []
        # Read in our samplesheet CSV
        with open("~{samplesheet}", 'r', newline='') as infile:
            reader = list(csv.DictReader(infile, delimiter=','))
            for row in reader:
                experiment_id = row.get("experiment_id", "")
                barcode = row["barcode"]
                row["merged_bam"] = barcode_to_bam.get(barcode, "")
                rows.append(row)

        DataTable_out_tsv = f"{experiment_id}__DataTable.tsv"
        with open(DataTable_out_tsv, 'w') as outf:
            fieldnames = list(row[0].keys())
            writer = csv.DictWriter(outf, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        DataTable_out_json = f"{experiment_id}__DataTable.json"
        with open(DataTable_out_json, "w") as outf:
            json.dump(rows, outf, indent=2)

        EOF
    >>>

    output {
        File samplesheet_with_bams = glob("*__DataTable.tsv")[0]
        File samplesheet_with_bams_json = glob("*__DataTable.json")[0]
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpu,
        mem_gb:             mem_gb,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
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

