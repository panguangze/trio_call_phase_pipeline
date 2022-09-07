# concat sv and snp
rule concat_sv_snp:
    input:
        calls=[rules.bcftools_sv_filter.output.out, rules.bcftools_nomiss.output.out],
        idx=rules.bcftools_sv_index.output.idx,
    output:
        concat_sv_snp_vcf="results/concat_sv_snp/{joint_calling_group}.concat.vcf.gz",
    log:
        "results/logs/concat_sv_snp/concat_sv_snp_{joint_calling_group}.log",
    params:
        extra="-a",
    wrapper:
        "0.75.0/bio/bcftools/concat"

rule bcftools_index:
    input:
        "results/concat_sv_snp/{joint_calling_group}.concat.vcf.gz",
    output:
        "results/concat_sv_snp/{joint_calling_group}.concat.vcf.gz.csi",
    params:
        extra=config["bcftools_index"]["extra"] + " --threads {}".format(
            config["bcftools_index"]["threads"]
        ),
    log:
        "results/logs/bcftools_index/{joint_calling_group}.log",
    threads: config["bcftools_index"]["threads"]
    wrapper:
        "0.75.0/bio/bcftools/index"

rule extractHairs:
    input:
        vcf="results/concat_sv_snp/{joint_calling_group}.concat.vcf.gz",
        vcf_csi="results/concat_sv_snp/{joint_calling_group}.concat.vcf.gz.csi",
        bams=lambda w: expand(
            "results/mapped/{sample}.bam",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
        idxs=lambda w: expand(
            "results/mapped/{sample}.bam.csi",
            sample=joint_calling_group_lists.loc[w.joint_calling_group],
        ),
    output:
        lst_done=touch("results/concat_sv_snp/{joint_calling_group}_lst/lst.done"),
    log:
        "results/logs/concat_sv_snp/concat_sv_snp_{joint_calling_group}.log",
    params:
        bin_path=config["extractHairs"]["bin_path"],
        ngs_params=config["extractHairs"]["ngs_params"],
        lst_dir="results/concat_sv_snp/{joint_calling_group}_lst",
    run:
        shell(
            "mkdir -p {params.lst_dir} && "
            "{params.bin_path} "
            "{params.ngs_params} "
            "--bam {input.bams[0]} "
            "--VCF {input.vcf} "
            "--idx 0 "
            "--out {params.lst_dir}/0.lst && sort -k3 -n {params.lst_dir}/0.lst > {params.lst_dir}/0.s.lst"
        )
        shell(
            "{params.bin_path} "
            "{params.ngs_params} "
            "--bam {input.bams[1]} "
            "--VCF {input.vcf} "
            "--idx 1 "
            "--out {params.lst_dir}/1.lst && sort -k3 -n {params.lst_dir}/1.lst > {params.lst_dir}/1.s.lst"
        )
        shell(
            "{params.bin_path} "
            "{params.ngs_params} "
            "--bam {input.bams[2]} "
            "--VCF {input.vcf} "
            "--idx 2 "
            "--out {params.lst_dir}/2.lst && sort -k3 -n {params.lst_dir}/2.lst > {params.lst_dir}/2.s.lst"
        )

rule spechap_trio:
    input:
        vcf="results/concat_sv_snp/{joint_calling_group}.concat.vcf.gz",
        lst_done=rules.extractHairs.output.lst_done,
    output:
        phased_dir=directory("results/spechap/{joint_calling_group}"),
    log:
        "results/logs/spechap/{joint_calling_group}.spec.log",
    params:
        lst_dir=lambda w,input: os.path.dirname(input.lst_done),
        bin_path=config["specHap"]["bin_path"],
        ngs_params=config["specHap"]["ngs_params"],
        tabix_path=config["other_bin"]["tabix"],
        bgzip_path=config["other_bin"]["bgzip"],
    run:
        shell(
            "mkdir -p {output.phased_dir} && "
            "{params.bin_path} "
            "{params.ngs_params} "
            "-v {input.vcf} "
            "-f {params.lst_dir}/0.s.lst "
            "-o  {output.phased_dir}/0.spec.vcf "
            "--idx 0  && {params.bgzip_path} {output.phased_dir}/0.spec.vcf" 
        )
        shell(
            "{params.bin_path} "
            "{params.ngs_params} "
            "-v {input.vcf} "
            "-f {params.lst_dir}/1.s.lst "
            "-o  {output.phased_dir}/1.spec.vcf "
            "--idx 1  && {params.bgzip_path} {output.phased_dir}/1.spec.vcf" 
        )
        shell(
            "{params.bin_path} "
            "{params.ngs_params} "
            "-v {input.vcf} "
            "-f {params.lst_dir}/2.s.lst "
            "-o  {output.phased_dir}/2.spec.vcf "
            "--idx 2  && {params.bgzip_path} {output.phased_dir}/2.spec.vcf" 
        )

rule tabix_spec:
    input:
        phased_dir=rules.spechap_trio.output.phased_dir
    output:
        temp("results/spechap/{joint_calling_group}_temp.txt")
    params:
        tabix_path=config["other_bin"]["tabix"],
    run:
        shell("{params.tabix_path} -f {input.phased_dir}/0.spec.vcf.gz")
        shell("{params.tabix_path} -f {input.phased_dir}/1.spec.vcf.gz")
        shell("{params.tabix_path} -f {input.phased_dir}/2.spec.vcf.gz")
        shell("touch {output}")

    