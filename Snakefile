chrN = [*map(str, range(1,22)), "M", "X", "Y"]
samples, = glob_wildcards("data/vcf/{sample}.vcf.gz")

rule all:
    input: "output/all-variants.exonic_variant_function"

rule clean:
    shell: "rm -rf output"

rule index_vcf_file:
    input: "{file}.vcf.gz"
    output: "{file}.vcf.gz.tbi"
    shell: "bcftools index {input:q} -t -o {output:q}"

rule download_annovar:
    output: "vendor/annovar.latest.tar.gz"
    shell: "open https://www.openbioinformatics.org/annovar/annovar_download_form.php"

rule download_hg19_ref_files:
    output: expand("vendor/hg19/chr{N}.fa.gz", N=chrN)
    shell: "rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ vendor/hg19"

rule uncompress_hg19_ref_files:
    input: "vendor/hg19/{chrX}.fa.gz"
    output: "vendor/hg19/uncompressed/{chrX}.fa"
    shell: "gunzip -c {input:q} > {output:q}"

rule merge_hg19_ref_files:
    input: expand("vendor/hg19/uncompressed/chr{N}.fa", N=chrN)
    output: "vendor/hg19/uncompressed/hg19_ref_genome.fa"
    shell: "cat {input:q} > {output:q}"

rule unpack_annovar:
    input: "vendor/annovar.latest.tar.gz"
    output: directory("vendor/annovar")
    shell: "tar -xvf {input} -C vendor"

rule sort_vcf_file:
    input: "data/vcf/{sample}.vcf.gz"
    output: "output/{sample}.sorted.vcf.gz"
    shell: "bcftools sort -Oz {input:q} -o {output:q}"

rule normalize_vcf_file:
    input: "output/{sample}.sorted.vcf.gz"
    output: "output/{sample}.normalized.vcf.gz"
    shell: "bcftools norm -m-both -Oz {input:q} -o {output:q}"

rule left_align_vcf_file:
    input:
        hg19="vendor/hg19/uncompressed/hg19_ref_genome.fa",
        vcf="output/{sample}.normalized.vcf.gz"
    output: "output/{sample}.leftaligned.vcf.gz"
    shell: "bcftools norm -f {input.hg19:q} {input:q}"

# HACK: removing wrong FORMAT/PL values
rule remove_format_pl_from_vcf_file:
    input: "output/{sample}.normalized.vcf.gz"
    output: "output/{sample}.noformatpl.vcf.gz"
    shell: "bcftools annotate -x FORMAT/PL -Oz {input:q} -o {output:q}"

rule merge_vcf_files:
    input:
        vcf = expand("output/{sample}.noformatpl.vcf.gz", sample=samples),
        tbi = expand("output/{sample}.noformatpl.vcf.gz.tbi", sample=samples)
    output: "output/all-variants.vcf"
    shell: "bcftools merge {input.vcf:q} -o {output:q}"

rule convert_vcf_to_annovar:
    input:
        av_c2a="vendor/annovar/convert2annovar.pl",
        vcf="{file}.vcf"
    output: "{file}.avinput"
    shell: "{input.av_c2a} -format vcf4 -allsample -withfreq -includeinfo {input.vcf:q} > {output:q}"

rule annotate_merged_vcf_file:
    input:
        av_hg19_knowngene="vendor/annovar/humandb/hg19_refGene.txt",
        av_annotate="vendor/annovar/annotate_variation.pl",
        avfile="output/all-variants.avinput"
    output:
        "output/all-variants.variant_function",
        "output/all-variants.exonic_variant_function"
    shell: "{input.av_annotate} -out output/all-variants -build hg19 {input.avfile} $(dirname {input.av_hg19_knowngene})"