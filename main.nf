#!/usr/bin/env nextflow


process pydeseq_main {
	conda '/home/acote/miniconda3/envs/nf-pydeseq'

	input:
		tuple val(output_prefix), path(output_dir)

	output:
		tuple val(output_prefix), path(dds_file)

	tag "${output_prefix}"

	script:
	"""
	python $launchDir/run_deseq.py ${params.output_prefix} -c ${params.count_matrix} -m ${params.meta} -d ${params.design_factors} -n ${params.n_cpus} -o ${output_dir}
	
	"""
}

process get_all_contrasts {
	conda '/home/acote/miniconda3/envs/nf-pydeseq'

	input:
		path(meta)

	tag "${params.output_prefix}:${params.design_factors}"

	output:
		stdout

	script:
	"""
	python $launchDir/generate_contrasts.py ${meta} ${params.design_factors}
	
	"""
}

process pydeseq_contrast {
	label 'bigmem'
	conda '/home/acote/miniconda3/envs/nf-pydeseq'

	input:
		each(contrast)

	tag "${params.output_prefix}"

	script:
	"""
	python $launchDir/run_contrasts.py ${params.deseq_dataset} ${params.design_factors} "${contrast}" ${params.output_dir} ${params.output_prefix}
	
	"""
}


process get_all_contrasts_r {
	conda '/home/acote/miniconda3/envs/rdeseq'

	input:
		path(meta)

	tag "${meta}:${params.meta_column}"

	output:
		stdout

	script:
	"""
	Rscript $launchDir/generate_contrasts.R ${meta} ${params.meta_column}
	
	"""
}

process run_rdeseq {
	conda '/home/acote/miniconda3/envs/rdeseq'

	input:
		each(contrast)

	tag "${contrast}"
	
	script:
	"""
	echo "${contrast}"
	Rscript $launchDir/run_rdeseq.R ${params.index_dir} ${params.meta_file} ${params.meta_column} "${contrast}" ${params.n_cpus} ${params.output_dir}
	"""
}


workflow all_contrasts {
        Channel.fromPath(params.meta)
        | get_all_contrasts
	| flatMap(n -> n.split('\n'))
	| pydeseq_contrast
}

workflow rdeseq_allcontrasts {
        Channel.fromPath(params.meta_file)
        | get_all_contrasts_r
	| flatMap(n -> n.split('\n'))
	| run_rdeseq
}

workflow rdeseq_specifiedcontrast {
        Channel.fromPath(params.contrast_file)
	| flatMap(n -> n.split('\n'))
	| run_rdeseq
}

workflow rdeseq_generatecontrasts {
        Channel.fromPath(params.meta_file)
        | get_all_contrasts_r
	| view()
}

// workflow all_contrasts {
// 	Channel.fromPath(params.meta)
// 	| get_all_contrasts
// 	| flatten
// 	| pydeseq_contrast
// }



