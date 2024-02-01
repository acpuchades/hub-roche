samples, = glob_wildcards("data/vcf/{sample}.vcf.gz")

rule clean:
    shell: "rm -rf output"

rule sort_vcf_file:
    input: "data/vcf/{sample}.vcf.gz"
    output: "output/{sample}.sorted.vcf.gz"
    shell: "bcftools sort -Oz {input:q} -o {output:q}"

rule index_vcf_file:
    input: "output/{sample}.sorted.vcf.gz"
    output: "output/{sample}.sorted.vcf.gz.tbi"
    shell: "bcftools index {input:q} -t -o {output:q}"

rule merge_vcf_files:
    input:
        vcf = expand("output/{sample}.sorted.vcf.gz", sample=samples),
        tbi = expand("output/{sample}.sorted.vcf.gz.tbi", sample=samples)
    output: "output/variants.merged.vcf.gz"
    shell: "bcftools merge -Oz {input.vcf:q} -o {output:q}"

rule convert_vcf_to_plink:
    input: "output/variants.merged.vcf"
    output: multiext("output/variants", ".bed", ".bim", ".ped")
    shell: "plink --make-bed --vcf {input:q} --out output/variants.merged"

rule convert_vcf_to_plink2:
    input: "output/variants.merged.vcf"
    output: multiext("output/variants", ".pgen", ".pvar", ".psam")
    shell: "plink2 --make-pgen --vcf {input:q} --out output/variants"