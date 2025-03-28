version 1.0
import "../../structs/Structs.wdl"

task ResolveTelomeres {
    meta {
        description: "Bespoke task for linear plasmid telomere processing: find, extract, clip, filter and merge"
        author: "Michael J. Foster"
    }

    parameter_meta {
        reads: "Input reads file in FASTQ format, gzipped"
        sample_id: "sample_id"
        motif: "Telomere motif to search for [ Default: 'AWWTWATTWTTTATTAGTATA' ] "
    }

    input {
        File reads
        String sample_id
        String motif = "AWWTWATTWTTTATTAGTATA"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 40 + 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        NPROCS=$( grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}' )

        # Step 1: Find telomere reads
        echo "Finding telomere reads..."
        zcat "~{reads}" | seqkit locate -j "$NPROCS" -i -I -d -p "~{motif}" --bed > "~{sample_id}.telomere_reads.bed"

        # Step 2: Extract read_ids for telomeric reads
        echo "Extracting read IDs..."
        awk 'NR>1 {print $1}' "~{sample_id}.telomere_reads.bed" > "~{sample_id}.telomere_read_ids.txt"

        # Step 3: Extract reads by read_id and also extract non-telo reads by inverted seqkit grep
        echo "Extracting telomere reads..."
        zcat "~{reads}" | seqkit grep -f "~{sample_id}.telomere_read_ids.txt" > "~{sample_id}.raw_telomeres.fastq"
        zcat "~{reads}" | seqkit grep -v -f "~{sample_id}.telomere_read_ids.txt" > "~{sample_id}.no_telomeres.fastq"

        # Step 4: Clip em
        echo "Clipping telomere reads..."
        python3 /opt/clip_telomere_reads.py "~{sample_id}.telomeres.bed" "~{sample_id}.raw_telomeres.fastq" "~{sample_id}.clipped_telomeres.fastq"

        # Step 5: Merge clipped telomere reads and remove the unclipped parents
        echo "Filtering out telomere reads..."
        python3 /opt/merge_reads.py "~{sample_id}.no_telomeres.fastq" "~{sample_id}.clipped_telomeres.fastq" "~{sample_id}.fixed.fastq"

        # step 6: gzip reads.
        cat "~{sample_id}.fixed.fastq" | gzip -9 > "~{sample_id}.fixed.fastq.gz"
    >>>

    output {
        File telomere_bed = "~{sample_id}.telomere_reads.bed"
        File telomere_read_ids = "~{sample_id}.telomere_read_ids.txt"
        File telomere_fastq = "~{sample_id}.raw_telomeres.fastq"
        File clipped_telomeres = "~{sample_id}.clipped_telomeres.fastq"
        File fixed_reads = "~{sample_id}.fixed.fastq.gz"
    }
    # do not preempt
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "mjfos2r/telomere-tools:latest"
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