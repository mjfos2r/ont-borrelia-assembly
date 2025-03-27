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

task ParseSamplesheetToDataTable {
    meta {
        description: "Parse the samplesheet and extract per-column arrays, also create a full TSV with an added raw_bams column for easy use as a Cromwell DataTable."
        # I know this is busted but we've gotta see how this changes things.
    }

    parameter_meta {
        samplesheet: "CSV-formatted samplesheet with sample metadata. Formatted per ONT's guidelines except renamed column: alias => sample_id and output a DataTable for Cromwell input."
        file_paths: "file containing GCS paths to the decompressed and merged bam files for addition to our DataTable"
    }

    input {
        File samplesheet
        File file_paths

        Int num_cpu = 4
        Int mem_gb = 16
        RuntimeAttr? runtime_attr_override
    }
    Int disk_size = 50 + ceil(size(samplesheet, "GB"))
    command <<<
        set -euxo pipefail

        cat "~{file_paths}" > gcs_paths.txt

        python3 <<EOF
        import os
        import csv
        import json
        barcode_to_bam = {}
        with open("gcs_paths.txt", 'r') as f:
            for line in f:
                path = line.strip()#.replace('"', '') # sanitize any potential quote issues.
                print(f"Input Path From Cromwell: {path}")
                filename = os.path.basename(path)
                barcode = filename.split(".")[0] # barcode01.merged.bam
                barcode_to_bam[barcode] = path
        experiment_id = ""
        rows = []
        # Read in our samplesheet CSV (utf-8-sig since excel has to be quirky)
        with open("~{samplesheet}", 'r', newline='', encoding='utf-8-sig') as infile:
            reader = csv.DictReader(infile, delimiter=',')
            for row in reader:
                # do this to check for bad samplesheet column naming.
                if "Bb_sample_id" not in row.keys():
                    row['Bb_sample_id'] = row.pop('sample_id')
                experiment_id = row.get("experiment_id", "")
                barcode = row["barcode"]
                merged_bam = barcode_to_bam.get(barcode, "")
                row["merged_bam"] = merged_bam
                rows.append(row)
                print(f"experiment_id: {experiment_id}")
                print(f"barcode: {barcode}")
                print(f"merged_bam: {merged_bam}\n")
                print(f"SampleSheet Columns: {row.keys()}")
        DataTable_out_tsv = "DataTable.tsv"
        print(DataTable_out_tsv)
        with open(DataTable_out_tsv, 'w') as outf:
            fieldnames = [ 'Bb_sample_id', 'barcode', 'experiment_id', 'flow_cell_id', 'position_id', 'flow_cell_product_code', 'kit', 'merged_bam']
            #fieldnames = rows[0].keys()
            print(f"Fieldnames: {fieldnames}")
            writer = csv.DictWriter(outf, delimiter='\t', fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        DataTable_out_json = "DataTable.json"
        print(DataTable_out_json)
        with open(DataTable_out_json, "w") as outf:
            json.dump(rows, outf, indent=2)
        EOF
    >>>

    output {
        File samplesheet_with_bams = "DataTable.tsv"
        File samplesheet_with_bams_json = "DataTable.json"
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

