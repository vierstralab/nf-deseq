#!/usr/bin/env nextflow


process pydeseq_main {
	conda '/home/acote/miniconda3/envs/pydeseq2'

	input:
		tuple val(output_prefix), path(output_dir)

	output:
		tuple val(output_prefix), path(dds_file)

	tag "${output_prefix}"

	script:
	"""
	python run_deseq.py ${params.output_prefix} -c ${params.count_matrix} -m ${params.meta} -d ${params.design_factors} -n ${params.n_cpus} -o ${output_dir}
	
	"""
}

process get_all_contrasts {
	conda '/home/acote/miniconda3/envs/pydeseq2'

	input:
		tuple path(meta)

	tag "${params.output_prefix}:${params.design_factors}"

	output:
		stdout

	script:
	"""
	python generate_contrasts.py ${meta} ${params.design_factors}
	
	"""
}

process pydeseq_contrast {
	conda '/home/acote/miniconda3/envs/pydeseq2'

	input:
		tuple val(contrast)

	tag "${params.output_prefix}:${contrast}"

	script:
	"""
	python run_contrasts.py ${params.deseq_dataset} ${params.design_factors} ${contrast} ${params.output_dir} ${params.output_prefix}
	
	"""
}

workflow all_contrasts {
	Channel.fromPath(params.meta)
	| get_all_contrasts
	| flatten()
	| pydeseq_contrast()
}



