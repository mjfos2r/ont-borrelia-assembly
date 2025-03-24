version 1.0

import "../tasks/QC/NanoPlot.wdl" as NP
import "../tasks/utilities/GeneralUtils.wdl" as GenUtils
import "../tasks/utilities/SamplesheetUtils.wdl" as SSUtils

workflow ONT_DemuxAndQC {
    meta {
        description: "Take in the tarball of the bam_pass file for our run, decompress it, merge all of the bams for each barcode, rename the barcode, then trim, filter, and generate QC report for all of our samples. Also generate a NanoPlot from the ONT summary.txt file"
        author: "Michael J. Foster"
    }

    parameter_meta {
        RunTarball: "Gzipped tarball of bam_pass generated by the sequencer"
        RunChecksum: "md5sum of the tarball containing all of the directories"
        summary_files: "sequencing summary file(s) generated by dorado"
        summary_checksums: "checksum file(s) for our sequencing summary file(s) generated by dorado"
        samplesheet: "samplesheet containing sample information, most importantly, the barcode_id and sample_alias"
    }

    input {
        File samplesheet
        File RunTarball
        File RunChecksum
        Array[File] summary_files
        Array[File] summary_checksums
    }

    ##### workflow ONTSummaryValidateAndQC {
    #####     input {
    #####         Array[File] summary_files = summary_files
    #####         Array[File] checksums = summary_checksums
    #####     }
    # Get the filenames for our summary files. create an array to use in generating the final validation array..
    scatter (file in summary_files){
        String summary_filename = basename(file)
    }
    Array[String] summary_filenames = summary_filename

    # now we validate our summary file(s)
    scatter (idx in range(length(summary_files))) {
        call GenUtils.ValidateMd5sum as summary_validation {
            input:
                file = summary_files[idx],
                checksum = summary_checksums[idx]
        }
        Boolean summary_is_valid = read_boolean(summary_validation.is_valid)
    }
    # gather our validation statuses. it'll be useful when it fails due to corruption that can't be dealt with in a conditional check.
    Array[Boolean] summary_validity = summary_is_valid
    Array[Pair[String, Boolean]] summary_integrity = zip(summary_filenames, summary_validity)

    # begone validation checks. If it's invalid just check after the execution fails...
    call NP.NanoPlotFromSummary { input: summary_files = summary_files }

    # Validate tarball's integrity but power through regardless since cromwell hates me and won't do a simple boolean conditional...
    call GenUtils.ValidateMd5sum as run_bams_validation { input: file = RunTarball, checksum = RunChecksum }
    Boolean run_is_valid = read_boolean(run_bams_validation.is_valid)
    call GenUtils.DecompressRunTarball { input: tarball = RunTarball }
    call SSUtils.ParseSamplesheetToDataTable {
        input:
            samplesheet = samplesheet,
            raw_bam_paths = DecompressRunTarball.raw_bam_paths
        }
    }

    output {
        # Metadata parsed from samplesheet
        Array[String] flow_cell_id = ParseSamplesheetToDataTable.flow_cell_id
        Array[String] position_id = ParseSamplesheetToDataTable.position_id
        Array[String] experiment_id = ParseSamplesheetToDataTable.experiment_id
        Array[String] flow_cell_product_code = ParseSamplesheetToDataTable.flow_cell_product_code
        Array[String] kit = ParseSamplesheetToDataTable.kit
        Array[String] barcode = ParseSamplesheetToDataTable.barcode
        Array[String] sample_id = ParseSamplesheetToDataTable.sample_id
        Array[String] bam_paths = ParseSamplesheetToDataTable.bam_paths
        File? samplesheet_with_bams = ParseSamplesheetToDataTable.samplesheet_with_bams

        # Validation output
        File run_validation_file = run_bams_validation.is_valid
        Boolean run_validation_status = run_is_valid

        # decompressed outputs from DecompressRunTarball if md5 is valid.
        Int directory_count = DecompressRunTarball.directory_count
        Array[Int] bam_counts = DecompressRunTarball.bam_counts
        Array[String] barcode_dirs = DecompressRunTarball.barcode_dirs
        Array[File] bam_lists = DecompressRunTarball.bam_lists
        Array[Array[File]] raw_bam_paths = DecompressRunTarball.raw_bam_paths
        File? bam_paths_json = DecompressRunTarball.bam_paths_json

        Array[Pair[String, Boolean]] summary_file_integrity = summary_integrity
        # NanoPlot outputs
        File nanoplot_map = NanoPlotFromSummary.map
        File  nanoplot_tarball = NanoPlotFromSummary.tarball
        Map[String, Float]  nanoplot_stats_map = NanoPlotFromSummary.stats_map
    }
}
