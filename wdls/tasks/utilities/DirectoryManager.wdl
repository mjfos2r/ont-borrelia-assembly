version 1.0

import "../../structs/Structs.wdl"
import "../../structs/DirectoryStructure.wdl"

# workflow to expose the directory structure
workflow SampleDirs {
    input {
        String gcs_bucket
        String sample_id
    }

    call InitializeDirectoryStructure {
        input:
            gcs_bucket = gcs_bucket,
            sample_id = sample_id
    }

    output {
        DirectoryStruct dirs = InitializeDirectoryStructure.dirs
    }
}

# init the structure given a bucket url to output to and the sample_id
task InitializeDirectoryStructure {
    input {
        String gcs_bucket
        String sample_id
    }

    command {
        echo "Initializing directory structure for ${sample_id}"
    }

    output {
        DirectoryStruct dirs = object {
            bucket_url: gcs_bucket,
            sample_id: sample_id,
            root_path: sub(gcs_bucket + "/" + sample_id, "/+", "/"),

            alignments: object {
                root: sub(gcs_bucket + "/" + sample_id + "/alignments", "/+", "/"),
                ReadsToRef: sub(gcs_bucket + "/" + sample_id + "/alignments/ReadsToRef", "/+", "/"),
                ReadsToAsm: sub(gcs_bucket + "/" + sample_id + "/alignments/ReadsToAsm", "/+", "/"),
                AsmToRef: sub(gcs_bucket + "/" + sample_id + "/alignments/AsmToRef", "/+", "/"),
                misc: sub(gcs_bucket + "/" + sample_id + "/alignments/misc", "/+", "/")
            },

            reads: object {
                root: sub(gcs_bucket + "/" + sample_id + "/reads", "/+", "/"),
                raw: sub(gcs_bucket + "/" + sample_id + "/reads/raw", "/+", "/"),
                cleaned: sub(gcs_bucket + "/" + sample_id + "/reads/cleaned", "/+", "/"),
                fixed: sub(gcs_bucket + "/" + sample_id + "/reads/fixed", "/+", "/")
            },

            contigs: object {
                root: sub(gcs_bucket + "/" + sample_id + "/contigs", "/+", "/"),
                draft: sub(gcs_bucket + "/" + sample_id + "/contigs/draft", "/+", "/"),
                polished: sub(gcs_bucket + "/" + sample_id + "/contigs/polished", "/+", "/"),
                misc: sub(gcs_bucket + "/" + sample_id + "/contigs/misc", "/+", "/")
            },

            annotations: object {
                root: sub(gcs_bucket + "/" + sample_id + "/annotations", "/+", "/"),
                bakta: sub(gcs_bucket + "/" + sample_id + "/annotations/bakta", "/+", "/"),
                plasmids: sub(gcs_bucket + "/" + sample_id + "/annotations/plasmids", "/+", "/"),
                lipoproteins: sub(gcs_bucket + "/" + sample_id + "/annotations/lipoproteins", "/+", "/"),
                ospC: sub(gcs_bucket + "/" + sample_id + "/annotations/ospC", "/+", "/"),
                mlst: sub(gcs_bucket + "/" + sample_id + "/annotations/mlst", "/+", "/"),
                vls: sub(gcs_bucket + "/" + sample_id + "/annotations/vls", "/+", "/"),
                misc: sub(gcs_bucket + "/" + sample_id + "/annotations/misc", "/+", "/")
            },

            reports: object {
                root: sub(gcs_bucket + "/" + sample_id + "/reports", "/+", "/"),
                preprocessing: sub(gcs_bucket + "/" + sample_id + "/reports/preprocessing", "/+", "/"),
                qc: sub(gcs_bucket + "/" + sample_id + "/reports/qc", "/+", "/"),
                assembly: sub(gcs_bucket + "/" + sample_id + "/reports/assembly", "/+", "/"),
                annotations: sub(gcs_bucket + "/" + sample_id + "/reports/annotations", "/+", "/"),
                summary: sub(gcs_bucket + "/" + sample_id + "/reports/summary", "/+", "/"),
                misc: sub(gcs_bucket + "/" + sample_id + "/reports/misc", "/+", "/")
            },

            logs: object {
                root: sub(gcs_bucket + "/" + sample_id + "/logs", "/+", "/"),
                alignment: sub(gcs_bucket + "/" + sample_id + "/logs/alignment", "/+", "/"),
                assembly: sub(gcs_bucket + "/" + sample_id + "/logs/assembly", "/+", "/"),
                polishing: sub(gcs_bucket + "/" + sample_id + "/logs/polishing", "/+", "/"),
                annotation: sub(gcs_bucket + "/" + sample_id + "/logs/annotation", "/+", "/"),
                system: sub(gcs_bucket + "/" + sample_id + "/logs/system", "/+", "/"),
                misc: sub(gcs_bucket + "/" + sample_id + "/logs/misc", "/+", "/")
            },

            plots: object {
                root: sub(gcs_bucket + "/" + sample_id + "/plots", "/+", "/"),
                ReadsToRefCov: sub(gcs_bucket + "/" + sample_id + "/plots/ReadsToRefCov", "/+", "/"),
                ReadsToAsmCov: sub(gcs_bucket + "/" + sample_id + "/plots/ReadsToAsmCov", "/+", "/"),
                AsmToRefCov: sub(gcs_bucket + "/" + sample_id + "/plots/AsmToRefCov", "/+", "/"),
                heatmaps: sub(gcs_bucket + "/" + sample_id + "/plots/heatmaps", "/+", "/"),
                misc: sub(gcs_bucket + "/" + sample_id + "/plots/misc", "/+", "/")
            }
        }
    }

    runtime {
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 1 HDD"
        docker: "ubuntu:latest"
    }
}