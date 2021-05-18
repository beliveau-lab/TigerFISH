configfile: "config.yaml"

SAMPLES, = glob_wildcards("data/" + "{sample}.fa")

rule all:
    input: 
        expand("kmer_specificity_out/{sample}_kmer_filter_probes.txt",sample=SAMPLES)

rule repeat_ID:
    input:
        fa="data/{sample}.fa",
        jf=config["jf_file"],
        chr_path=config["chr_path"] + "{sample}.fa"

    params: 
        sliding_window = config["sliding_window"],
        threshold = config["threshold"],
        composition = config["composition"],
        file_start = config["file_start"],
        mfree="40G",
        h_rt="24:0:0"

    output:
        'repeat_id_out/{sample}_index.txt',
        'repeat_id_out/{sample}_regions.bed',
        'repeat_id_out/{sample}_regions.fa',
        'repeat_id_out/{sample}_jf_temp.txt'
    shell: 
        "python scripts/repeat_ID.py -j {input.jf} -f {input.fa} -chr {wildcards.sample} -s {params.sliding_window} -t {params.threshold} -c {params.composition} -st {params.file_start} -schr {input.chr_path}"
        
rule design_probes:
    input: 
        region_fa="repeat_id_out/{sample}_regions.fa"

    output:
       "designed_probes_out/{sample}_blockParse_probe_df.bed"

    params:
        mfree="25G",
        h_rt="24:0:0"

    shell: 
        "python scripts/design_probes.py -f {input.region_fa} -chr {wildcards.sample}"

rule specificity: 
    input:
       	jf='repeat_id_out/{sample}_jf_temp.txt',
       	idx="repeat_id_out/{sample}_index.txt",
       	probes="designed_probes_out/{sample}_blockParse_probe_df.bed",
        region_fa="repeat_id_out/{sample}_regions.fa"

    params:
        mfree="40G",
        h_rt="24:0:0",
        merlength=config['merlength'],
        enrich_score=config['enrich_score'],
        copy_num=config['copy_num']

    output: 
        "initial_specificity_out/{sample}_probes_final_score.txt"
    
    shell: 
        "python scripts/specificity_refactor.py -p {input.probes} -j {input.jf} -k {input.idx} -c {wildcards.sample} -f {input.region_fa} -m {params.merlength} -e {params.enrich_score} -cn {params.copy_num}"

rule kmer_specificity:
    input:
        filtered_probes="initial_specificity_out/{sample}_probes_final_score.txt"

    params:
        k_val=config["k_val"],
        mfree="25G",
        h_rt="24:0:0"
         
    output: 
        "kmer_specificity_out/{sample}_kmer_filter_probes.txt"
    
    shell: 
        "python scripts/kmer_decomp_filter.py -f {input.filtered_probes} -k {params.k_val} -c {wildcards.sample}"
