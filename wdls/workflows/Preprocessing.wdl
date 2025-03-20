version 1.0

import "../tasks/utilities/GeneralUtils.wdl" as GenUtils
import "../tasks/utilities/SamplesheetUtils.wdl" as SSUtils

workflow Preprocessing {
    meta {
        description: "Take in the tarball of the bam_pass file for our run, decompress it, merge all of the bams for each barcode, rename the barcode, then trim and clean the reads"
    }
    parameter_meta {
        RunTarball: "Gzipped tarball of bam_pass generated by the sequencer"
        RunChecksum: "md5sum of the tarball containing all of the directories"
        samplesheet: "Samplesheet containing sample information, most importantly, the barcode_id and sample_alias"
    }
    input {
        File RunTarball
        File RunChecksum
        File samplesheet
    }

    # Parse samplesheet to get barcode to sample name mapping
    call SSUtils.ParseSamplesheet { input: samplesheet = samplesheet }

    # Validate tarball, then decompress and process each barcode directory below
    call GenUtils.ValidateMd5sum { input: file = RunTarball, checksum = RunChecksum }

    # Only proceed with processing if tarball is valid
    if (ValidateMd5sum.is_valid) {
        # Decompress the run tarball to get access to barcode dirs and BAM files
        call GenUtils.DecompressRunTarball { input: tarball = RunTarball }
    }

    output {
        # Metadata parsed from samplesheet
        Array[String] position_id = ParseSamplesheet.position_id
        Array[String] experiment_id = ParseSamplesheet.experiment_id
        Array[String] flow_cell_product_code = ParseSamplesheet.flow_cell_product_code
        Array[String] kit = ParseSamplesheet.kit
        Array[String] barcode = ParseSamplesheet.barcodes
        Array[String] sample_id = ParseSamplesheet.sample_id

        # Validation output
        Boolean is_valid = ValidateMd5sum.is_valid
        # decompressed outputs from DecompressRunTarball if md5 is valid.
        Int? directory_count = DecompressRunTarball.directory_count
        Array[Int]? bam_counts = DecompressRunTarball.bam_counts
        Array[String]? barcode_dirs = DecompressRunTarball.barcode_dirs
        Array[File]? bam_lists = DecompressRunTarball.bam_lists
        Map[String, Array[File]]? barcode_to_bams_map = DecompressRunTarball.barcode_to_bams_map
    }
}
