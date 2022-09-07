def parse_sample_names(w,input):
    return list(joint_calling_group_lists.loc[w.joint_calling_group])
rule dysgu_run:
    input:
        bam="results/mapped/{sample}.bam",
        idx="results/mapped/{sample}.bam.csi",
    output:
        # gvcfs=expand("results/individual_calls/{joint_calling_group}_gvcf"),
        # vcfs=directory("results/individual_calls/{joint_calling_group}_vcf"),
        out_vcf="results/dysgu_sv/{sample}.dysgu.vcf",
        tmp_dir=temp(directory("results/dysgu_sv/{sample}_dysgu_tmp")),
        # scratch=directory("results/individual_calls/{joint_calling_group}_interm"),
        # vcfs=directory("results/all_group_samples_joint_calls/{joint_calling_group}_vcf"),
        # gvcfs=directory("results/all_group_samples_joint_calls/{joint_calling_group}_gvcf"),
    params:
        run_params=config["dysgu"]["run_params"],
        ref=config['ref']['fasta'],
        bin_path=config["dysgu"]["bin_path"],
    threads: config["dysgu"]["threads"]
    log:
        "results/logs/{sample}_dysgu.log",
    shell:
        "{params.bin_path} run "
        "-p {threads} "
        "{params.run_params} "
        "{params.ref} "
        "{output.tmp_dir} "
        "{input.bam} "
        "-o {output.out_vcf}"
        

rule dysgu_merge:
    input:
        vcfs=lambda w: expand(
            "results/dysgu_sv/{sample}.dysgu.vcf",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
    output:
        out_vcf="results/dysgu_sv/{joint_calling_group}.dysgu.merge.vcf",
    params:
        bin_path=config["dysgu"]["bin_path"],
        merge_trio_params=config["dysgu"]["merge_trio_params"]
    threads: config["dysgu"]["threads"]
    log:
        "results/logs/dysgu_sv/{joint_calling_group}_dysgu.merge.log",
    shell:
        "{params.bin_path} merge "
        "{params.merge_trio_params} "
        "{input.vcfs} "
        "-o {output.out_vcf}"


rule bcftools_sv_filter:
    input:
        vcf="results/dysgu_sv/{joint_calling_group}.dysgu.merge.vcf",
    output:
        out="results/dysgu_sv/{joint_calling_group}.dysgu.merge.filter.vcf.gz",
    log:
        "results/logs/dysgu_sv/bcftools_sv_filter_{joint_calling_group}.log",
    params:
        filter=config["bcftools_filter"]["filter"],
        extra="",
    wrapper:
        "0.75.0/bio/bcftools/filter"

rule bcftools_sv_index:
    input:
        "results/dysgu_sv/{joint_calling_group}.dysgu.merge.filter.vcf.gz",
    output:
        idx="results/dysgu_sv/{joint_calling_group}.dysgu.merge.filter.vcf.gz.csi",
    params:
        extra=config["bcftools_index"]["extra"] + " --threads {}".format(
            config["bcftools_index"]["threads"]
        ),
    log:
        "results/logs/bcftools_index/{joint_calling_group}_sv.log",
    threads: config["bcftools_index"]["threads"]
    wrapper:
        "0.75.0/bio/bcftools/index"