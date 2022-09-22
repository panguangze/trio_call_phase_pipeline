def parse_sample_names(w,input):
    return list(joint_calling_group_lists.loc[w.joint_calling_group])
rule deepvariant_gvcf:
    input:
        bams=lambda w: expand(
            "results/mapped/{sample}.bam",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
        idxs=lambda w: expand(
            "results/mapped/{sample}.bam.csi",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
        ref=config['ref']['fasta'],
    output:
        # gvcfs=expand("results/individual_calls/{joint_calling_group}_gvcf"),
        # vcfs=directory("results/individual_calls/{joint_calling_group}_vcf"),
        out_dir=directory("results/deepvariant_gvcf/{joint_calling_group}_dgcall"),
        # scratch=directory("results/individual_calls/{joint_calling_group}_interm"),
        # vcfs=directory("results/all_group_samples_joint_calls/{joint_calling_group}_vcf"),
        # gvcfs=directory("results/all_group_samples_joint_calls/{joint_calling_group}_gvcf"),
    params:
        sample_names=parse_sample_names,
        # config=config["deepvariant_gvcf"]["config"],
        ref=config['ref']['fasta'],
    threads: config["deepvariant_gvcf"]["threads"]
    log:
        "results/logs/deepvariant_gvcf/{joint_calling_group}/stdout.log",
    singularity:
        "docker://google/deepvariant:deeptrio-latest"
    shell:
        "run_deeptrio "
        "--model_type WGS "
        "--ref {params.ref} "
        "--reads_child {input.bams[0]} "
        "--reads_parent1 {input.bams[1]} "
        "--reads_parent2 {input.bams[2]} "
        "--output_vcf_child {output.out_dir}/{params.sample_names[0]}.vcf.gz "
        "--output_vcf_parent1 {output.out_dir}/{params.sample_names[1]}.vcf.gz "
        "--output_vcf_parent2 {output.out_dir}/{params.sample_names[2]}.vcf.gz "
        "--sample_name_child {params.sample_names[0]} "
        "--sample_name_parent1 {params.sample_names[1]} "
        "--sample_name_parent2 {params.sample_names[2]} "
        "--num_shards {threads}  "
        "--intermediate_results_dir {output}/intermediate_results_dir "
        "--output_gvcf_child {output.out_dir}/{params.sample_names[0]}.g.vcf.gz "
        "--output_gvcf_parent1 {output.out_dir}/{params.sample_names[1]}.g.vcf.gz "
        "--output_gvcf_parent2 {output.out_dir}/{params.sample_names[2]}.g.vcf.gz "
        


rule glnexus:
    input:
        gvcfs=rules.deepvariant_gvcf.output.out_dir,
    output:
        vcf="results/glnexus_gvcf/{joint_calling_group}.vcf.gz",
        scratch=temp(
            directory("results/glnexus_gvcf/{joint_calling_group}.DB")
        ),
    params:
        config=config["glnexus"]["config"],
        sample_names=parse_sample_names,
    threads: config["glnexus"]["threads"]
    log:
        "results/logs/glnexus/{joint_calling_group}/stdout.log",
    container:
        "docker://quay.io/mlin/glnexus:v1.3.1"
    shell:
        "glnexus_cli "
        "--config DeepVariantWGS "
        "--dir {output.scratch} "
        "--threads {threads} "
        "{input.gvcfs}/{params.sample_names[0]}.g.vcf.gz "
        "{input.gvcfs}/{params.sample_names[1]}.g.vcf.gz "
        "{input.gvcfs}/{params.sample_names[2]}.g.vcf.gz "
        "2> {log} "
        "| bcftools view - "
        "| bgzip -c "
        "> {output.vcf} "


rule bcftools_index:
    input:
        "{vcffile}.vcf.gz",
    output:
        "{vcffile}.vcf.gz.csi",
    params:
        extra=config["bcftools_index"]["extra"] + " --threads {}".format(
            config["bcftools_index"]["threads"]
        ),
    log:
        "results/logs/bcftools_index/{vcffile}.log",
    threads: config["bcftools_index"]["threads"]
    wrapper:
        "master/bio/bcftools/index"


rule bcftools_nomiss:
    input:
        vcf="results/glnexus_gvcf/{joint_calling_group}.vcf.gz",
    output:
        out="results/glnexus_gvcf/{joint_calling_group}.nomiss.vcf.gz",
    log:
        "results/logs/bcftools_nomiss_{joint_calling_group}.log",
    params:
        extra="-g ^miss -Oz",
    wrapper:
        "master/bio/bcftools/view"

# rule bcftools_filter:
#     input:
#         vcf="results/glnexus_gvcf/{joint_calling_group}.nomiss.vcf.gz",
#     output:
#         "results/glnexus_gvcf/{joint_calling_group}.filter.vcf.gz",
#     log:
#         "results/logs/bcftools_filter_{joint_calling_group}.log",
#     params:
#         filter=config["bcftools_filter"]["filter"],
#         extra="",
#     wrapper:
#         "master/bio/bcftools/filter"
