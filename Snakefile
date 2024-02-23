chrN = [*map(str, range(1,23)), "M", "X", "Y"]
samples, = glob_wildcards("data/vcf/{sample}.vcf.gz")

rule all:
    input: multiext("output/all-variants", ".variant_function", ".exonic_variant_function")

rule clean:
    shell: "rm -rf output"

rule index_vcf_file:
    input: "{file}.vcf.gz"
    output: "{file}.vcf.gz.tbi"
    shell: "bcftools index {input:q} -t -o {output:q}"

rule download_annovar:
    output: "vendor/annovar.latest.tar.gz"
    shell: "open https://www.openbioinformatics.org/annovar/annovar_download_form.php"

rule download_ucsc_hg38_ref_files:
    output: directory("vendor/hg38")
    shell: "wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/*' -P vendor/hg38"

rule uncompress_genome_ref_files:
    input: "vendor/{genome}/{chrX}.fa.gz"
    output: "vendor/{genome}/uncompressed/{chrX}.fa"
    shell: "gunzip -c {input:q} > {output:q}"

rule merge_hg38_ref_files:
    input: expand("vendor/hg38/uncompressed/chr{N}.fa", N=chrN)
    output: "vendor/hg38/uncompressed/hg38_ref_genome.fa"
    shell: "cat {input:q} > {output:q}"

rule unpack_annovar:
    input: "vendor/annovar.latest.tar.gz"
    output: directory("vendor/annovar")
    shell: "tar -xvf {input} -C vendor"

rule sort_vcf_file:
    input: "data/vcf/{sample}.vcf.gz"
    output: "output/{sample}.sorted.vcf.gz"
    shell: "bcftools sort -Oz {input:q} -o {output:q}"

rule spread_multiline_vcf_file:
    input: "output/{sample}.sorted.vcf.gz"
    output: "output/{sample}.spread.vcf.gz"
    shell: "bcftools norm -m-both -Oz {input:q} -o {output:q}"

rule left_align_vcf_file:
    input:
        hg38="vendor/hg38/uncompressed/hg38_ref_genome.fa",
        vcf="output/{sample}.spread.vcf.gz"
    output: "output/{sample}.normalized.vcf.gz"
    shell: "bcftools norm -f {input.hg38:q} -Oz {input.vcf:q} -o {output:q}"

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

rule download_annovar_dbfile:
    output: "vendor/annovar/humandb/{dbfile}"
    shell: "wget -q -O - http://www.openbioinformatics.org/annovar/download/{wildcards.dbfile}.gz | gunzip > {output:q}"

rule annotate_merged_vcf_file_genes:
    input:
        av_annotate="vendor/annovar/annotate_variation.pl",
        av_hg38_knowngene="vendor/annovar/humandb/hg38_knownGene.txt",
        av_hg38_knowngene_kgxref="vendor/annovar/humandb/hg38_kgXref.txt",
        av_hg38_knowngene_mrna="vendor/annovar/humandb/hg38_knownGeneMrna.fa",
        avfile="output/all-variants.avinput"
    output:
        "output/all-variants.variant_function",
        "output/all-variants.exonic_variant_function"
    shell: "{input.av_annotate} -out output/all-variants -geneanno -build hg38 {input.avfile} vendor/annovar/humandb -dbtype knownGene"