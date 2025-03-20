# fast mode for trycycler pipeline where it's assembling read subsets.
task Bakta {

    parameter_meta {
        prefix:         "Prefix to pass to bakta"
        contigs:        "contigs to annotate"
        rename_table:   "contig renaming table using our plasmid calls"
        bakta_db:       "tar.gz containing the full bakta database for annotation"
        genus:          "[ Default: Borrelia ] genus to pass to bakta"
        species:        "[ Default: burgdorferi ] species to pass to bakta"
    }

    input {
        String prefix
        File contigs
        File rename_table
        File bakta_db
        String genus = "Borrelia"
        String species = "burgdorferi"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 15 + 2 * ceil(size(trimmed_reads, "GB"))

    command <<<
        set -euxo pipefail
        NPROC=$(awk '/^processor/{print}' /proc/cpuinfo | wc -l)
        tar --use-compress-program=pigz -xf "~{tarball}" -C bakta_db
        BAKTA_DB="bakta_db"

        bakta \
            --threads "$NPROC"
            --prefix "~{prefix}" \
            --db "$BAKTA_DB" \
            --replicons "~{rename_table}" \
            --genus "~{genus}" \
            --species "~{species}"
            --gram - \
            "~{contigs}"
    >>>

    output {
        File annotations_tsv = "~{prefix}.tsv"
        File gff3 = "~{prefix}.gff3"
        File gbff = "~{prefix}.gbff"
        File embl = "~{prefix}.embl"
        File fna = "~{prefix}.fna"
        File ffn = "~{prefix}.ffn"
        File faa = "~{prefix}.faa"
        File inference_metrics_tsv = "~{prefix}.inference.tsv"
        File hypotheticals_tsv = "~{prefix}.hypotheticals.tsv"
        File hypotheticals_faa = "~{prefix}.hypotheticals.faa"
        File summary_txt = "~{prefix}.txt"
        File png = "~{prefix}.png"
        File svg = "~{prefix}.svg"
        File summary_json = "~{prefix}.json"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          32,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "mjfos2r/bakta:latest"
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