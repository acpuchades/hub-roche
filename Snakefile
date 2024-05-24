samples, = glob_wildcards("data/vcf/{sample}.vcf.gz")

auto_chrN = list(map(str, range(1, 23)))
nuclear_chrN = auto_chrN + ['X', 'Y']
all_chrN = nuclear_chrN + ['M']

sys.path.insert(1, './src')

rule all:
    input:
        "output/debug/vcf.txt",
        "output/debug/vcf-normalized.txt",
        "output/debug/vcf.unmerged.txt",
        "output/phenotypes/ms-assoc.pheno",
        "output/phenotypes/als-assoc.pheno",
        "output/common-variants/plink.bed",
        "output/common-variants/plink.eigenvec"

rule clean:
    shell: "rm -rf output"

rule generate_control_phenotypes_file:
    input:
        helpers="src/helpers.py",
        script="src/make-pheno-ufela.py",
        samples="data/ufela/samples-20240201.xlsx"
    output: "output/phenotypes/controls.pheno"
    shell: "{input.script} --samples {input.samples:q} --output controls > {output:q}"

rule generate_als_phenotypes_file:
    input:
        helpers="src/helpers.py",
        script="src/make-pheno-ufela.py",
        samples="data/ufela/samples-20240201.xlsx",
        ufela_db="data/ufela/formulario_2023-11-15.sqlite"
    output: "output/phenotypes/als.pheno"
    shell: "{input.script} --samples {input.samples:q} --database {input.ufela_db:q} --output patients > {output:q}"

rule generate_ms_phenotypes_file:
    input:
        helpers="src/helpers.py",
        script="src/make-pheno-ufem.py",
        samples="data/ufem/samples-20240201.xlsx",
        nhc="data/ufem/FClinica.xlsx",
        edmus_db="data/ufem"
    output: "output/phenotypes/ms.pheno"
    shell: "{input.script} --samples {input.samples:q} --nhc {input.nhc:q} --database {input.edmus_db:q} > {output:q}"

rule generate_assoc_phenotype_file:
    input:
        cases="output/phenotypes/{pheno}.pheno",
        controls="output/phenotypes/controls.pheno"
    output: "output/phenotypes/{pheno}-assoc.pheno"
    shell: "\
        printf FID\\\\tIID\\\\t{wildcards.pheno}\\\\n | tr '[:lower:]' '[:upper:]' > {output:q} && \
        tail -n +2 {input.cases:q} | awk 'BEGIN {{OFS=\"\\t\"}} {{print $1, $2, 2}}' >> {output:q} && \
        tail -n +2 {input.controls:q} | awk 'BEGIN {{OFS=\"\\t\"}} {{print $1, $2, 1}}' >> {output:q} \
    "

rule generate_vcf_debug_file:
    input: "output/renamed-variants/samples.original.txt"
    output: "output/debug/vcf.txt"
    shell: "cp {input:q} {output:q}"

rule generated_vcf_normalized_debug_file:
    input: "output/renamed-variants/samples.normalized.txt"
    output: "output/debug/vcf-normalized.txt"
    shell: "cp {input:q} {output:q}"

rule generate_pheno_debug_file:
    input: "output/phenotypes/{pheno}.pheno"
    output: "output/debug/{pheno}.txt"
    shell: "tail -n +2 {input:q} | cut -d $'\t' -f2 | sort -k 1,1 > {output:q}"

rule generate_vcf_unmerged_debug_file:
    input:
        als_pheno="output/debug/als.txt",
        ms_pheno="output/debug/ms.txt",
        controls_pheno="output/debug/controls.txt",
        vcf_normalized="output/debug/vcf-normalized.txt"
    output: "output/debug/vcf.unmerged.txt"
    shell: "cat {input.vcf_normalized:q} | grep -f {input.als_pheno:q} -v | grep -f {input.ms_pheno:q} -v | grep -f {input.controls_pheno:q} -v | sort -n -k1,1 > {output:q}"

rule compress_file_bgzip:
    input: "output/{path}"
    output: "output/{path}.gz"
    shell: "bgzip --keep {input:q}"

rule index_vcf_file_tbi:
    input: "output/{path}.vcf.{ext}"
    output: "output/{path}.vcf.{ext,b?gz}.tbi"
    shell: "bcftools index {input:q} -t -o {output:q}"

rule index_bcf_file_tbi:
    input: "output/{path}.bcf"
    output: "output/{path}.bcf.tbi"
    shell: "bcftools index {input:q} -t -o {output:q}"

rule download_annovar:
    output: "vendor/annovar.latest.tar.gz"
    shell: "open https://www.openbioinformatics.org/annovar/annovar_download_form.php"

rule download_ucsc_hg38_ref_files:
    output: directory("vendor/genomes/{genome,^uncompressed}")
    shell: "wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.genome}/chromosomes/*' -P {output:q}"

rule concat_hg38_ref_files:
    input: "vendor/genomes/{genome}"
    output: "vendor/genomes/uncompressed/{genome}.fa"
    shell: "find {input:q} -name '*.fa.gz' -exec gunzip -c {{}} >> {output:q} \\;"

rule unpack_annovar:
    input: "vendor/annovar.latest.tar.gz"
    output: directory("vendor/annovar")
    shell: "tar -xvf {input} -C vendor"

rule sort_vcf_file:
    input: "data/vcf/{sample}.vcf.gz"
    output: "output/sorted-variants/{sample}.vcf.gz"
    shell: "bcftools sort -Oz {input:q} -o {output:q}"

rule spread_multiline_vcf_file:
    input: "output/sorted-variants/{sample}.vcf.gz"
    output: "output/spread-variants/{sample}.vcf.gz"
    shell: "bcftools norm -m-both -Oz {input:q} -o {output:q}"

rule left_align_vcf_file:
    input:
        hg38="vendor/genomes/uncompressed/hg38.fa",
        vcf="output/spread-variants/{sample}.vcf.gz"
    output: "output/normalized-variants/{sample}.vcf.gz"
    shell: "bcftools norm -f {input.hg38:q} -Oz {input.vcf:q} -o {output:q}"

rule remove_format_pl_from_vcf_file:
    input: "output/normalized-variants/{sample}.vcf.gz"
    output: "output/fixed-variants/{sample}.vcf.gz"
    shell: "bcftools annotate -x FORMAT/PL -Oz {input:q} -o {output:q}"

rule merge_vcf_files:
    input:
        vcf = expand("output/fixed-variants/{sample}.vcf.gz", sample=samples),
        tbi = expand("output/fixed-variants/{sample}.vcf.gz.tbi", sample=samples)
    output: "output/merged-variants/merged.vcf"
    shell: "bcftools merge {input.vcf:q} -Ov -o {output:q}"

rule download_annovar_dbfile:
    output: "vendor/annovar/humandb/{dbfile}"
    shell: "wget -q -O - http://www.openbioinformatics.org/annovar/download/{wildcards.dbfile}.gz | gunzip > {output:q}"

rule annovar_annotate_merged_variants:
    input:
        av_tableav="vendor/annovar/table_annovar.pl",
        av_hg38_knowngene="vendor/annovar/humandb/hg38_knownGene.txt",
        av_hg38_knowngene_kgxref="vendor/annovar/humandb/hg38_kgXref.txt",
        av_hg38_knowngene_mrna="vendor/annovar/humandb/hg38_knownGeneMrna.fa",
        av_hg38_dbnsfp30a="vendor/annovar/humandb/hg38_dbnsfp30a.txt",
        vcf="output/merged-variants/merged.vcf"
    output: multiext("output/annovar-annotated-variants/merged.hg38_multianno", ".vcf", ".txt")
    shell: "{input.av_tableav} -buildver hg38 {input.vcf:q} \
            -out output/annovar-annotated-variants/merged \
            -protocol refGene,cytoBand,dbnsfp30a \
            -operation g,r,f -nastring . -vcfinput -polish \
            vendor/annovar/humandb"

rule download_vep_cache_files:
    output: "vendor/vep/homo_sapiens_vep_111_GRCh38.tar.gz"
    shell: "curl https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz -o {output:q}"

rule extract_vep_cache_files:
    input: "vendor/vep/homo_sapiens_vep_111_GRCh38.tar.gz"
    output: directory("vendor/vep/homo_sapiens")
    shell: "tar xvf {input:q} -C vendor/vep"

rule vep_annotate_merged_variants:
    input:
        vep_cache="vendor/vep/homo_sapiens",
        hg38="vendor/genomes/uncompressed/hg38.fa",
        vcf="output/merged-variants/merged.vcf"
    output: "output/vep-annotated-variants/merged.vcf"
    shell: "vep --cache --dir vendor/vep -i {input.vcf:q} --fasta {input.hg38:q} --sift b --polyphen b --hgvs --vcf -o {output:q}"

rule split_vep_annotated_variant_files:
    input:
        vcf="output/vep-annotated-variants/merged.vcf.gz",
        tbi="output/vep-annotated-variants/merged.vcf.gz.tbi"
    output: "output/vep-annotated-variants/split.{chrN}.vcf.gz"
    shell: "bcftools view -Oz {input.vcf:q} -o {output:q} --regions {wildcards.chrN}"

rule download_gnomad_vcf_files:
    output: "vendor/gnomad/{version}/{type}/{file}"
    shell: "curl 'https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.version}/vcf/{wildcards.type}/{wildcards.file}' -o {output:q}"

rule download_clinvar_files:
    output: "vendor/clinvar/{build}/{file}"
    shell: "curl 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{wildcards.build}/{wildcards.file}' -o {output:q}"

rule download_dbsnp_files:
    output: "vendor/dbsnp/GRCh38/{file}"
    shell: "curl 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/{wildcards.file}' -o {output:q}"

rule interpolate_vcfanno_gnomad_file_nuclear:
    input:
        config="config/gnomad.conf.in",
        gnomad_vcf="vendor/gnomad/4.0/genomes/gnomad.genomes.v4.0.sites.{chrN}.vcf.bgz"
    output: "output/gnomad-annotated-variants/vcfanno.{chrN}.conf"
    shell: "GNOMAD_VCF_FILE={input.gnomad_vcf:q} envsubst < {input.config:q} > {output:q}"

rule gnomad_annotate_nuclear_variants:
    input:
        gnomad_vcf="vendor/gnomad/4.0/genomes/gnomad.genomes.v4.0.sites.{chrN}.vcf.bgz",
        gnomad_tbi="vendor/gnomad/4.0/genomes/gnomad.genomes.v4.0.sites.{chrN}.vcf.bgz.tbi",
        vcfanno_config="output/gnomad-annotated-variants/vcfanno.{chrN}.conf",
        vcf="output/vep-annotated-variants/split.{chrN}.vcf.gz"
    output: "output/gnomad-annotated-variants/annotated.{chrN}.vcf"
    shell: "vcfanno -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule concat_gnomad_annotated_variant_files:
    input:
        vcf=expand("output/gnomad-annotated-variants/annotated.chr{N}.vcf.gz", N=auto_chrN),
        tbi=expand("output/gnomad-annotated-variants/annotated.chr{N}.vcf.gz.tbi", N=auto_chrN),
        chrM_vcf="output/vep-annotated-variants/split.chrM.vcf.gz",
        chrM_tbi="output/vep-annotated-variants/split.chrM.vcf.gz.tbi",
    output: "output/gnomad-annotated-variants/annotated.bcf"
    shell: "bcftools concat {input.vcf:q} {input.chrM_vcf:q} -Ob -o {output:q}"

rule generate_chr_conversion_file:
    output: "output/gnomad-annotated-variants/chr-conv.txt"
    run:
        with open(output[0], "w") as out:
            for n in all_chrN:
                print(f'chr{n}\t{n}', file=out)

rule rename_chrs_in_gnomad_annotate_variants:
    input:
        vcf="output/gnomad-annotated-variants/annotated.bcf",
        chr_conv="output/gnomad-annotated-variants/chr-conv.txt"
    output: "output/gnomad-annotated-variants/chr_renamed.vcf.gz"
    shell: "bcftools annotate --rename-chrs {input.chr_conv:q} {input.vcf:q} -Oz -o {output:q}"

rule interpolate_vcfanno_clinvar_dbsnp_file:
    input:
        config="config/clinvar_dbsnp.conf.in",
        clinvar_vcf="vendor/clinvar/GRCh38/clinvar_20240221.vcf.gz",
        dbsnp_vcf="vendor/dbsnp/GRCh38/All_20180418.vcf.gz"
    output: "output/clinvar-annotated-variants/vcfanno.conf"
    shell: "CLINVAR_VCF_FILE={input.clinvar_vcf:q} DBSNP_VCF_FILE={input.dbsnp_vcf:q} envsubst < {input.config:q} > {output:q}"

rule clinvar_annotate_variants:
    input:
        clinvar_vcf="vendor/clinvar/GRCh38/clinvar_20240221.vcf.gz",
        clinvar_tbi="vendor/clinvar/GRCh38/clinvar_20240221.vcf.gz.tbi",
        dbsnp_vcf="vendor/dbsnp/GRCh38/All_20180418.vcf.gz",
        dbsnp_tbi="vendor/dbsnp/GRCh38/All_20180418.vcf.gz.tbi",
        vcfanno_config="output/clinvar-annotated-variants/vcfanno.conf",
        vcf="output/gnomad-annotated-variants/chr_renamed.vcf.gz"
    output: "output/clinvar-annotated-variants/annotated.vcf"
    shell: "vcfanno -p {workflow.cores} -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule extract_variant_file_sample_names:
    input: "output/clinvar-annotated-variants/annotated.vcf"
    output: "output/renamed-variants/samples.original.txt"
    shell: "bcftools query --list-samples {input:q} > {output:q}"

rule normalize_sample_names:
    input:
        helpers="src/helpers.py",
        sample_names="output/renamed-variants/samples.original.txt"
    output: "output/renamed-variants/samples.normalized.txt"
    run:
        from helpers import normalize_sample_id
        with open(input.sample_names, 'r') as f:
            with open(output[0], 'w') as out:
                for line in f.readlines():
                    name = line.rstrip('\n')
                    name = normalize_sample_id(name)
                    print(name, file=out)

rule rename_variant_file_samples:
    input:
        vcf="output/clinvar-annotated-variants/annotated.vcf",
        normalized_names="output/renamed-variants/samples.normalized.txt"
    output: "output/renamed-variants/renamed.vcf"
    shell: "bcftools reheader -s {input.normalized_names:q} -o {output:q} {input.vcf:q}"

rule sort_normalized_sample_names:
    input: "output/renamed-variants/samples.normalized.txt"
    output: "output/renamed-variants/samples.sorted.txt"
    shell: "sort -k1,1 {input:q} > {output:q}"

rule sort_renamed_variant_file_samples:
    input:
        vcf="output/renamed-variants/renamed.vcf",
        sorted_names="output/renamed-variants/samples.sorted.txt"
    output: "output/renamed-variants/sorted.bcf"
    shell: "bcftools view -S {input.sorted_names} -Ob -o {output:q} {input.vcf:q}"

rule filter_common_variants_by_maf:
    input: "output/renamed-variants/sorted.bcf"
    output: "output/common-variants/maf_filtered.ge.0.05.bcf"
    shell: "bcftools view -i 'gno_af_all > 0.05' -Ob -o {output:q} {input:q}"

rule filter_annotated_variants_with_maf_lesser_than:
    input: "output/renamed-variants/sorted.bcf"
    output: "output/rare-variants/maf_filtered.lt.{threshold}.bcf"
    shell: "bcftools view -i 'gno_af_all < {wildcards.threshold}' -Ob -o {output:q} {input:q}"

rule filter_rare_variants_maf_unknown:
    input: "output/renamed-variants/sorted.bcf"
    output: "output/rare-variants/maf_filtered.unknown.bcf"
    shell: "bcftools view -e 'gno_af_all' -Ob -o {output:q} {input:q}"

rule convert_common_variants_to_plink:
    input: "output/common-variants/maf_filtered.ge.0.05.bcf"
    output: multiext("output/common-variants/plink", ".bed", ".bim", ".fam")
    shell: "plink --bcf {input:q} --allow-extra-chr --out output/common-variants/plink"

rule extract_principal_components_from_variants:
    input: "output/common-variants/plink.bed"
    output: multiext("output/common-variants/plink", ".eigenvec", ".eigenval")
    shell: "plink --bfile output/common-variants/plink --pca --out output/common-variants/plink"