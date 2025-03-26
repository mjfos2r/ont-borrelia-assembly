version 1.0
workflow CallPlasmids { meta {
			author: "Michael J. Foster"
			description: "Bespoke task to classify and annotate Borrelia burgdorferi plasmids in draft genome assemblies. Written by Michael J. Foster in the Lemieux Lab @ MGH/BroadInstitute"
		}
		parameter_meta {
			sample_id: "sample_id for the contigs we're classifying."
			contigs: "Fasta file of contigs/sequences to classify against our WP_db and PF32_db."
		}
		input {
			String sample_id
			File contigs
		}
		command <<<
		call_plasmids.py --contigs "~{contigs}"
		>>>
		output {

		}

}
