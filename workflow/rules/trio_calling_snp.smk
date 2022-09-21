rule deepvariant_gvcf:
    input:
        bam=unpack(get_bam),
        ref=config['ref']['fasta'],
    output:
        vcf="results/individual_calls/{sample}.vcf.gz",
        gvcf="results/individual_calls/{sample}.g.vcf.gz",
        tmp_dir=temp(directory("results/individual_calls/{sample}_tmp")),
    params:
        model=config["deepvariant_gvcf"]["model"],
        extra=config["deepvariant_gvcf"]["extra"],
    threads: config["deepvariant_gvcf"]["threads"]
    log:
        "results/logs/deepvariant_gvcf/{sample}/stdout.log",
    singularity:
        "docker://google/deepvariant:latest"
    shell:
        "run_deepvariant "
        "--model_type {params.model} "
        "--ref {input.ref} "
        "--reads {input.bam} "
        "--output_vcf {output.vcf} "
        "--num_shards {threads}  "
        "--intermediate_results_dir {output.tmp_dir} "
        "--output_gvcf {output.gvcf} 2> {log}"

rule glnexus:
    input:
        gvcfs=lambda w: expand(
            "results/individual_calls/{sample}.g.vcf.gz",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
    output:
        vcf=temp("results/all_group_samples_joint_calls/{joint_calling_group}.vcf.gz"),
        scratch=temp(
            directory("results/all_group_samples_joint_calls/{joint_calling_group}.DB")
        ),
    params:
        config=config["glnexus"]["config"],
    threads: config["glnexus"]["threads"]
    log:
        "results/logs/glnexus/{joint_calling_group}/stdout.log",
    singularity:
        "docker://quay.io/mlin/glnexus:v1.3.1"
    shell:
        "glnexus_cli "
        "--config DeepVariantWGS "
        "--dir {output.scratch} "
        "--threads {threads} "
        "{input} "
        "2> {log} "
        "| bcftools view - "
        "| bgzip -c "
        "> {output.vcf} "


rule create_reheader_sample_file:
    input:
        joint_calling_groups=config["joint_calling_groups"],
    output:
        samples=temp("results/joint_calls/{joint_calling_group}_sample_names.tsv"),
    log:
        "results/logs/reheader_sample_file/{joint_calling_group}.log",
    run:
        (
            joint_calling_groups.assign(
                group_sample=lambda x: x.group.str.cat(x.sample_id, sep=":")
            )
            .loc[
                lambda x: x.group == wildcards.joint_calling_group,
                ["sample_id", "group_sample"],
            ]
            .to_csv(output.samples, sep="\t", index=False, header=None)
        )


rule update_sample_names:
    input:
        vcf=rules.glnexus.output.vcf,
        samples=rules.create_reheader_sample_file.output.samples,
    output:
        vcf="results/joint_calls/{joint_calling_group}.vcf.gz",
    log:
        "results/logs/update_sample_names/{joint_calling_group}.log",
    params:
        extra="",
        view_extra="-O z",
    wrapper:
        "master/bio/bcftools/reheader"

rule bcftools_index:
    input:
        "results/joint_calls/{joint_calling_group}.vcf.gz",
    output:
        "results/joint_calls/{joint_calling_group}.vcf.gz.csi",
    params:
        extra=config["bcftools_index"]["extra"] + " --threads {}".format(
            config["bcftools_index"]["threads"]
        ),
    log:
        "results/logs/bcftools_index/{joint_calling_group}.log",
    threads: config["bcftools_index"]["threads"]
    wrapper:
        "master/bio/bcftools/index"


# rule bcftools_merge:
#     input:
#         calls=[
#             *expand(
#                 "results/joint_calls/{joint_calling_group}.vcf.gz",
#                 joint_calling_group=joint_calling_group_lists.index,
#             ),
#         ],
#         idxs=[
#             *expand(
#                 "results/joint_calls/{joint_calling_group}.vcf.gz.csi",
#                 joint_calling_group=joint_calling_group_lists.index,
#             ),
#         ],
#     output:
#         calls=temp("results/merged_calls/all.unfiltered.vcf.gz"),
#     log:
#         "results/logs/bcftools_merge/bcftools_merge.log",
#     params:
#         config["bcftools_merge"]["params"] + " -Oz",  # optional parameters for bcftools concat (except -o)
#     wrapper:
#         "master/bio/bcftools/merge"



rule bcftools_nomiss:
    input:
        vcf="results/joint_calls/{joint_calling_group}.vcf.gz",
        idx=rules.bcftools_index.output
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
