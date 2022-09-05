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
        samples=lambda w: expand(
            "{sample}",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
        ref=config['ref']['fasta'],
    output:
        gvcfs=expand("results/individual_calls/{sample}.g.vcf.gz"),
        vcfs=expand("results/individual_calls/{sample}.vcf.gz"),
        scratch=directory("results/individual_calls/{joint_calling_group}_interm"),
        # vcfs=directory("results/all_group_samples_joint_calls/{joint_calling_group}_vcf"),
        # gvcfs=directory("results/all_group_samples_joint_calls/{joint_calling_group}_gvcf"),
    params:
        config=config["glnexus"]["config"],
    threads: config["glnexus"]["threads"]
    log:
        "results/logs/glnexus/{joint_calling_group}/stdout.log",
    container:
        "docker://hub.docker.com/google/deepvariant:deeptrio-latest"
    shell:
        "echo 'text'"

        


rule glnexus:
    input:
        gvcfs=lambda w: expand(
            "results/individual_calls/{sample}.g.vcf.gz",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
    output:
        vcf="results/individual_calls/{joint_calling_groupple}.vcf.gz",
        scratch=temp(
            directory("results/individual_calls/{joint_calling_group}.DB")
        ),
    params:
        config=config["glnexus"]["config"],
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
        "{input} "
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
        "0.75.0/bio/bcftools/index"


rule bcftools_merge:
    input:
        calls=[
            *expand(
                "results/individual_calls/{sample}.vcf.gz",
            ),
        ],
        idxs=[
            *expand(
                "results/individual_calls/{sample}.vcf.gz.csi",
            ),
        ],
    output:
        calls=temp("results/merged_calls/all.unfiltered.vcf.gz"),
    log:
        "results/logs/bcftools_merge/bcftools_merge.log",
    params:
        config["bcftools_merge"]["params"] + " -Oz",  # optional parameters for bcftools concat (except -o)
    wrapper:
        "0.75.0/bio/bcftools/merge"


rule bcftools_filter:
    input:
        rules.bcftools_merge.output.calls,
    output:
        "results/merged_calls/all.vcf.gz",
    log:
        "results/logs/bcftools_filter_all.log",
    params:
        filter=config["bcftools_filter"]["filter"],
        extra="",
    wrapper:
        "0.75.0/bio/bcftools/filter"
