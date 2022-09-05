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
        vcf="results/individual_calls/{joint_calling_group}.vcf.gz",
        report=report(
            "results/individual_calls/{joint_calling_group}.visual_report.html",
            caption="../report/vcf.rst",
            category="Calls",
        ),
        gvcf="results/individual_calls/{joint_calling_group}.g.vcf.gz",
        scratch=temp(
            directory("results/all_group_samples_joint_calls/{joint_calling_group}_interm")
        ),
    params:
        config=config["glnexus"]["config"],
    threads: config["glnexus"]["threads"]
    log:
        "results/logs/glnexus/{joint_calling_group}/stdout.log",
    container:
        "docker://hub.docker.com/google/deepvariant:deeptrio-latest"
    shell:
        "/opt/deepvariant/bin/deeptrio/run_deeptrio "
        "--model_type WGS "
        "--ref {input.ref} "
        "--reads_child {input.bams[0]} "
        "--reads_parent1 {input.bams[1]} "
        "--reads_parent2 {input.bams[2]} "
        "--output_vcf_child {output.vcf} "
        "--output_vcf_parent1 {output.vcf} "
        "--output_vcf_parent2 {output.vcf} "
        "--sample_name_child {input.samples[0]} "
        "--sample_name_parent1 {input.samples[0]} "
        "--sample_name_parent2 {input.samples[0]} "
        "--num_shards {threads}  "
        "--intermediate_results_dir {output.scratch} "
        "--output_gvcf_child {output.gvcf} "
        "--output_gvcf_parent1 {output.gvcf} "
        "--output_gvcf_parent2 {output.gvcf}"
        


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
        "0.75.0/bio/bcftools/reheader"


rule bcftools_merge:
    input:
        calls=[
            *expand(
                "results/individual_calls/{sample}.vcf.gz",
                sample=(
                    samples.loc[
                        samples.sample_id.isin(joint_calling_groups.sample_id)
                    ].sample_id.unique()
                ),
            ),
            *expand(
                "results/joint_calls/{joint_calling_group}.vcf.gz",
                joint_calling_group=joint_calling_group_lists.index,
            ),
        ],
        idxs=[
            *expand(
                "results/individual_calls/{sample}.vcf.gz.csi",
                sample=(
                    samples.loc[
                        samples.sample_id.isin(joint_calling_groups.sample_id)
                    ].sample_id.unique()
                ),
            ),
            *expand(
                "results/joint_calls/{joint_calling_group}.vcf.gz.csi",
                joint_calling_group=joint_calling_group_lists.index,
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
