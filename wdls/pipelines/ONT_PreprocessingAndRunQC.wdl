version 1.0

import "../tasks/QC/NanoPlot.wdl" as NP
import "../tasks/utilities/GeneralUtilities.wdl" as GenUtils
import "../workflows/Preprocessing.wdl" as Preproc

workflow ONT_PreprocessingAndRunQC {
    meta {
        description: "Take in the tarball of the bam_pass file for our run, decompress it, merge all of the bams for each barcode, rename the barcode, then trim, filter, and generate QC report for all of our samples. Also generate a NanoPlot from the ONT summary.txt file"
    }

    parameter_meta {
        RunTarball: "Gzipped tarball of bam_pass generated by the sequencer"
        RunChecksum: "md5sum of the tarball containing all of the directories"
        summary_files: "sequencing summary file(s) generated by dorado"
        summary_checksums: "checksum file(s) for our sequencing summary file(s) generated by dorado"
        samplesheet: "samplesheet containing sample information, most importantly, the barcode_id and sample_alias"
    }

    input {
        File RunTarball
        File RunChecksum
        Array[File] summary_files
        Array[File] summary_checksums
        File samplesheet
    }

    # Run Preprocessing workflow
    call Preproc.Preprocessing {
        input:
            RunTarball = RunTarball,
            RunChecksum = RunChecksum,
            samplesheet = samplesheet,
    }

    scatter (idx in range(length(summary_files))) {
        call GenUtils.ValidateMd5sum as summary_validation {
            input:
                file = summary_files[idx],
                checksum = summary_checksums[idx]
        }
    }
    Array[Boolean] summary_valid = summary_validation.is_valid
    Array[Boolean] false_values = select_all(boolean_array, false)
    if (length(false_values) = 0 ) {
        call NP.NanoPlotFromSummary { input: summary_file = summary_files }
    }
    Array[Pair[String, Boolean]] summary_integrity = zip(summary_files, summary_valid)

    output {
        # Run-level metadata from Preprocessing
        ## Metadata parsed from samplesheet
        Array[String] position_id = Preprocessing.position_id
        Array[String] experiment_id = Preprocessing.experiment_id
        Array[String] flow_cell_product_code = Preprocessing.flow_cell_product_code
        Array[String] kit = Preprocessing.kit
        Array[String] barcode = Preprocessing.barcode
        Array[String] sample_id = Preprocessing.sample_id
        # Validation output
        Boolean md5_validation_passed = Preprocessing.is_valid
        Array[Pair[String, Boolean]] summary_file_integrity = summary_integrity
        # decompressed outputs from DecompressRunTarball if md5 is valid.
        Int? directory_count = Preprocessing.directory_count
        Array[Int]? bam_counts = Preprocessing.bam_counts
        Array[String]? barcode_dirs = Preprocessing.barcode_dirs
        Array[File]? bam_lists = Preprocessing.bam_lists
        Map[String, Array[File]]? barcode_to_bams_map = Preprocessing.barcode_to_bams_map

        # NanoPlot outputs
        File? nanoplot_map = NanoPlotFromSummary.map
        File? nanoplot_tarball = NanoPlotFromSummary.tarball
        Map[String, Float]? nanoplot_stats_map = NanoPlotFromSummary.stats_map
    }
}
