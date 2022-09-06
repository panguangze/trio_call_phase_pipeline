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
        # gvcfs=expand("results/individual_calls/{joint_calling_group}_gvcf"),
        # vcfs=directory("results/individual_calls/{joint_calling_group}_vcf"),
        gvcfs=expand(
            "results/individual_calls/{sample}.g.vcf.gz",
            sample=(
                samples.loc[
                    ~samples.sample_id.isin(joint_calling_groups.sample_id)
                ].sample_id.unique()
            ),
        ),
        # scratch=directory("results/individual_calls/{joint_calling_group}_interm"),
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
        gvcfs=ruels.deepvariant_gvcf.gvcfs,
    output:
        vcf="results/individual_calls/{joint_calling_group}.vcf.gz",
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


rule bcftools_filter:
    input:
        vcf="results/individual_calls/{joint_calling_group}.vcf.gz"
    output:
        "results/merged_calls/all.vcf.gz",
    log:
        "results/logs/bcftools_filter_all.log",
    params:
        filter=config["bcftools_filter"]["filter"],
        extra="",
    wrapper:
        "0.75.0/bio/bcftools/filter"
