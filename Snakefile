configfile: "config.yml"

rule all:
    input:
        expand("results/lda_specificity_out/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}_l{local_val}_g{global_val}/{sample}_global_kept.txt", sample=config["samples"], window=config["sliding_window"], threshold=config["threshold"], composition=config["composition"],enrich_score=config["enrich_score"], copy_num = config["copy_num"], local_val = config["local_val"], global_val = config["global_val"])

rule generate_jf_idx:
    input:
        fasta_file = "../../../../../../../reference/Assemblies/hg38_bp_noheader/{sample}.fa",
        jf = config["jf_file"],
        chr_path = config["chr_path"] + "{sample}.fa"

    conda:
        'envs/tigerfish_pipeline.yml'

    params:
        chrom_name = "{sample}",
        file_start = config["file_start"],
        mfree="100G",
        h_rt="120:0:0"

    benchmark:
        "results/benchmarks/generate_jf_idx/{sample}_log.log"

    output:
        "results/jf_index_files/{sample}_jf_temp.txt",
        "results/jf_index_files/{sample}_index.txt"

    shell:
        "python ../../bin/generate_jf_idxs.py -f {input.fasta_file} -j {input.jf} -chr {params.chrom_name} -schr {input.chr_path} -st {params.file_start}"

rule repeat_ID:
    input:
        jf_count = "results/jf_index_files/{sample}_jf_temp.txt",
        chrom_index = "results/jf_index_files/{sample}_index.txt",
        chr_path = config["chr_path"] + "{sample}.fa"

    conda:
        'envs/tigerfish_pipeline.yml'

    params:
        window = "{window}",
        threshold = "{threshold}",
        composition = "{composition}",
        file_start = config["file_start"],
        chrom_name = "{sample}",
        mfree="100G",
        h_rt="120:0:0"

    benchmark:
        "results/benchmarks/repeat_ID/w{window}_t{threshold}_c{composition}/{sample}_log.log"

    output:
        out_bed = "results/repeat_id_out/w{window}_t{threshold}_c{composition}/{sample}_regions.bed",
        out_fasta = "results/repeat_id_out/w{window}_t{threshold}_c{composition}/{sample}_regions.fa"

    shell:
        "python ../../bin/refactor_repeatID.py -j {input.jf_count} -i {input.chrom_index} -t {params.threshold} -c {params.composition} -chr {params.chrom_name}  -st {params.file_start} -schr {input.chr_path} -o_b {output.out_bed} -o_f {output.out_fasta}"

rule design_probes:
    input:
        region_fa = "results/repeat_id_out/w{window}_t{threshold}_c{composition}/{sample}_regions.fa"

    output:
        "results/designed_probes_out/w{window}_t{threshold}_c{composition}/{sample}_blockParse_probe_df.bed"

    conda:
        'envs/tigerfish_pipeline.yml'

    params:
        mfree="100G",
        h_rt="120:0:0",
        window = "{window}",
        threshold = "{threshold}",
        composition = "{composition}",
        chrom_name = "{sample}"

    benchmark:
        "results/benchmarks/design_probes/w{window}_t{threshold}_c{composition}/{sample}_log.log"

    shell:
        "python ../../bin/design_probes.py -f {input.region_fa} -chr {params.chrom_name} -win {params.window} -thresh {params.threshold} -comp {params.composition}"

rule specificity:
    input:
        jf = "results/jf_index_files/{sample}_jf_temp.txt",
        probes = "results/designed_probes_out/w{window}_t{threshold}_c{composition}/{sample}_blockParse_probe_df.bed",
        region_fa = "results/repeat_id_out/w{window}_t{threshold}_c{composition}/{sample}_regions.fa"

    conda:
        'envs/tigerfish_pipeline.yml'

    params:
        mfree="100G",
        h_rt="120:0:0",
        merlength = config["merlength"],
        enrich_score = "{enrich_score}",
        copy_num = "{copy_num}",
        window = "{window}",
        threshold = "{threshold}",
        composition = "{composition}",
        chrom_name = "{sample}"

    benchmark:
        "results/benchmarks/specificity/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}/{sample}_log.log"

    output:
        "results/initial_specificity_out/pre_filter/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}/{sample}_probes_final_score.txt",
        "results/initial_specificity_out/post_filter/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}/{sample}_probes_final_score.txt"

    shell:
        "python ../../bin/kmer_filter.py -p {input.probes} -j {input.jf} -ch {params.chrom_name} -f {input.region_fa}  -w {params.window} -t {params.threshold} -c {params.composition} -m {params.merlength} -e {params.enrich_score} -cn {params.copy_num}"

rule lda_specificity:
    input:
        filtered_probes = "results/initial_specificity_out/post_filter/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}/{sample}_probes_final_score.txt",
        pre_filtered_probes = "results/initial_specificity_out/pre_filter/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}/{sample}_probes_final_score.txt"

    conda:
       'envs/tigerfish_pipeline.yml'

    params:
        local_val = "{local_val}",
        global_val = "{global_val}",
        mfree="80G",
        h_rt="120:0:0",
        enrich_score = "{enrich_score}",
        copy_num = "{copy_num}",
        window = "{window}",
        threshold = "{threshold}",
        composition = "{composition}",
        chrom_name = "{sample}"

    benchmark:
        "results/benchmarks/lda_specificity/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}_l{local_val}_g{global_val}/{sample}_log.log"

    output:
        "results/lda_specificity_out/w{window}_t{threshold}_c{composition}_e{enrich_score}_cn{copy_num}_l{local_val}_g{global_val}/{sample}_global_kept.txt"

    shell:
        "python ../../bin/bt2_allvall.py -f {input.filtered_probes} -o {params.chrom_name}_global_kept.txt  -w {params.window} -t {params.threshold} -c {params.composition} -e {params.enrich_score} -cn {params.copy_num} -t_l {params.local_val} -t_g {params.global_val}"
