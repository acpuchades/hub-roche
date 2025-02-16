SPLICEAI_GENOME_SCORES_DIR = "vendor/spliceai/genome_scores_v1.3-194103939/genome_scores_v1.3-ds.20a701bc58ab45b59de2576db79ac8d0"

samples = glob_wildcards("data/laboratorio/vcf/{sample}.vcf.gz")[0]

auto_chrN = [str(x) for x in range(1,23)]
nuclear_chrN = auto_chrN + ['X', 'Y']
all_chrN = nuclear_chrN + ['M']

all_pheno = ["ALS", "MS", "ALS_MS"]

cva_suffixes = [".qq.png", ".volcano.png", ".manhattan.png", ".bar.png"]

rva_tests = ["CMC", "VariableThresholdPrice", "Skat", "SkatO"]

rva_filters = ["impact_HIGH", "cadd_20", "sift4g_D", "polyphen2_D"]
rva_suffixes = [".qq.png", ".bar.png"]
rva_mafs = ["any", 0.01, 0.02, 0.05]

rule all:
    input:
        expand("output/phenotype-info/{pheno}.samples", pheno=["CONTROL", "ALS", "MS"]),
        "output/variant-analysis/ALL_PCA.eigenval.scree.png",
        "output/variant-analysis/CONTROL_HWE.hardy.ternary.png",
        expand(
           "output/variant-analysis/cva/results.{pheno}.glm.logistic.hybrid{suffix}",
           pheno=all_pheno, suffix=cva_suffixes
        ),
        expand(
            "output/variant-analysis/rva.{filter}.maf_{maf}/{pheno}.{test}.assoc{suffix}",
            pheno=all_pheno, filter=rva_filters, maf=rva_mafs,
             test=rva_tests, suffix=rva_suffixes
        )

rule clean:
    shell: "rm -rf output"

rule compress_output_file_bgzip:
    input: "output/{path}"
    output: "output/{path}.gz"
    shell: "bgzip --keep {input:q}"

rule convert_bcf_file_to_vcf:
    input: "output/{path}.bcf"
    output: "output/{path}.vcf"
    shell: "bcftools view {input:q} -Ov -o {output:q}"

rule convert_bcf_file_to_compressed_vcf:
    input: "output/{path}.bcf"
    output: "output/{path}.vcf.gz"
    shell: "bcftools view {input:q} -Oz -o {output:q}"

ruleorder:
    convert_bcf_file_to_compressed_vcf > compress_output_file_bgzip

rule index_vcf_file_tbi:
    input: "output/{path}.{ext}"
    output: "output/{path}.{ext}.tbi"
    shell: "bcftools index {input:q} -t -o {output:q}"

rule download_ucsc_hg38_ref_genome:
    output: "vendor/genomes/hg38/hg38.fa.gz"
    shell: "wget -O - 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' | gunzip > {output:q}"

rule left_align_vcf_file:
    input:
        hg38="vendor/genomes/hg38/hg38.fa",
        vcf="output/sorted-variants/{sample}.vcf.gz"
    output: "output/leftaligned-variants/{sample}.vcf.gz"
    shell: "bcftools norm -f {input.hg38:q} -Oz {input.vcf:q} -o {output:q}"

rule spread_multiallelic_vcf_file:
    input: "output/leftaligned-variants/{sample}.vcf.gz"
    output: "output/normalized-variants/{sample}.vcf.gz"
    shell: "bcftools norm -m -both -Oz {input:q} -o {output:q}"

rule filter_variants_qc:
    input: "data/laboratorio/vcf/{sample}.vcf.gz"
    output: "output/filtered-variants/{sample}.bcf"
    shell: "bcftools view -f PASS {input:q} -Ob -o {output:q}"

rule sort_vcf_files:
    input: "output/filtered-variants/{sample}.bcf"
    output: "output/sorted-variants/{sample}.bcf"
    shell: "bcftools sort {input:q} -Ob -o {output:q}"

rule remove_format_pl_from_vcf_file:
    input: "output/sorted-variants/{sample}.bcf"
    output: "output/noformatpl-variants/{sample}.vcf.gz"
    shell: "bcftools annotate -x FORMAT/PL {input:q} -Oz -o {output:q}"

rule merge_vcf_files:
    input:
        vcf=expand("output/noformatpl-variants/{sample}.vcf.gz", sample=samples),
        tbi=expand("output/noformatpl-variants/{sample}.vcf.gz.tbi", sample=samples)
    output: "output/merged-variants/merged.bcf"
    shell: "bcftools merge -m none --missing-to-ref {input.vcf:q} -Ob -o {output:q}"

rule extract_variant_file_sample_names:
    input: "output/merged-variants/merged.bcf"
    output: "output/renamed-variant-samples/samples.original.txt"
    shell: "bcftools query --list-samples {input:q} > {output:q}"

rule normalize_sample_names:
    input: "output/renamed-variant-samples/samples.original.txt"
    output: "output/renamed-variant-samples/samples.normalized.txt"
    shell: "sed -e 's/^EMADN/EM/' \
                -e 's/-ambBED$//' \
                -e 's/ELA1912ADN1/ELA1912-2/' \
                -e 's/ELA1923ADN1/ELA1923-2/' \
                -e 's/ELA1924ADN1/ELA1924-2/' \
                -e 's/ADN[12]$//' \
                -e 's/-//g' \
                {input:q} > {output:q}"

rule normalize_vcf_sample_headers:
    input:
        vcf="output/merged-variants/merged.bcf",
        normalized_names="output/renamed-variant-samples/samples.normalized.txt"
    output: "output/renamed-variant-samples/normalized.bcf"
    shell: "bcftools reheader -s {input.normalized_names} {input.vcf:q} -o {output:q}"

rule sort_sample_names:
    input: "output/renamed-variant-samples/samples.normalized.txt"
    output: "output/renamed-variant-samples/samples.sorted.txt"
    shell: "sort -k1,1 {input:q} > {output:q}"

rule reorder_vcf_sample_names_in_header:
    input:
        vcf="output/renamed-variant-samples/normalized.bcf",
        sorted_names="output/renamed-variant-samples/samples.sorted.txt"
    output: "output/renamed-variant-samples/reordered.vcf"
    shell: "bcftools view -S {input.sorted_names} -Ov -o {output:q} {input.vcf:q}"

rule download_vep_cache_files:
    output: "vendor/vep/homo_sapiens_vep_{ver}_{genome}.tar.gz"
    shell: "curl -k https://ftp.ensembl.org/pub/release-{wildcards.ver}/variation/indexed_vep_cache/homo_sapiens_vep_{wildcards.ver}_{wildcards.genome}.tar.gz -o {output:q}"

rule download_vep_refseq_files:
    output: "vendor/vep/homo_sapiens_vep_{ver}_GRCh38.dna.primary_assembly.fa.gz"
    shell: "curl -k https://ftp.ensembl.org/pub/release-{wildcards.ver}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -o {output:q}"

rule bgzip_compress_vep_refseq_files:
    input: "vendor/vep/homo_sapiens_vep_{ver}_GRCh38.dna.primary_assembly.fa.gz"
    output: "vendor/vep/homo_sapiens_vep_{ver}_GRCh38.dna.primary_assembly.fa.bgz"
    shell: "gunzip -c {input:q} | bgzip > {output:q}"

rule extract_vep_cache_files:
    input: "vendor/vep/homo_sapiens_vep_{version}.tar.gz"
    output: directory("vendor/vep/homo_sapiens/{version}")
    shell: "tar xvf {input:q} -C vendor/vep"

rule download_dbnsfp_database_files:
    output: "vendor/dbnsfp/dbNSFP4.8a.zip"
    shell: "wget https://dbnsfp.s3.amazonaws.com/dbNSFP4.8a.zip -O {output:q}"

rule extract_dbnsfp_database_files:
    input: "vendor/dbnsfp/dbNSFP4.8a.zip"
    output: directory("vendor/dbnsfp/4.8a")
    shell: "unzip -d {output:q} {input:q}"

rule process_dbnsfp_database_files_header:
    input: "vendor/dbnsfp/{version}/dbNSFP{version}_variant.chr1.gz"
    output: "vendor/vep/dbnsfp/dbNSFP{version}_hdr.txt"
    shell: "gzcat {input:q} | head -n1 -q > {output:q}"

rule process_dbnsfp_database_files_GRCh38:
    input:
        dbnsfp_hdr="vendor/vep/dbnsfp/dbNSFP{version}_hdr.txt",
        dbnsfp_chr=expand("vendor/dbnsfp/{version}/dbNSFP{version}_variant.chr{N}.gz", version="4.8a", N=all_chrN)
    output: "vendor/vep/dbnsfp/dbNSFP{version}_GRCh38.gz"
    shell: "zgrep -h -v '^#chr' {input.dbnsfp_chr:q} | sort -k1,1 -k2,2n | cat {input.dbnsfp_hdr:q} - | bgzip -c > {output:q}"

rule index_dbnsfp_database_processed_file_GRCh38:
    input: "vendor/vep/dbnsfp/dbNSFP{version}_GRCh38.gz"
    output: "vendor/vep/dbnsfp/dbNSFP{version}_GRCh38.gz.tbi"
    shell: "tabix -s 1 -b 2 -e 2 {input:q}"

rule vep_annotate_merged_variants:
    input:
        go_terms="vendor/vep/go_terms",
        vep_cache="vendor/vep/homo_sapiens/112_GRCh38",
        vep_refseq="vendor/vep/homo_sapiens_vep_112_GRCh38.dna.primary_assembly.fa.bgz",
        dbnsfp_db="vendor/vep/dbnsfp/dbNSFP4.8a_GRCh38.gz",
        dbnsfp_tbi="vendor/vep/dbnsfp/dbNSFP4.8a_GRCh38.gz.tbi",
        spliceai_snv=f"{SPLICEAI_GENOME_SCORES_DIR}/spliceai_scores.raw.snv.hg38.vcf.gz",
        spliceai_indel=f"{SPLICEAI_GENOME_SCORES_DIR}/spliceai_scores.raw.indel.hg38.vcf.gz",
        vcf="output/renamed-variant-samples/reordered.vcf"
    output: "output/vep-annotated-variants/annotated.vcf"
    shell: "vep --cache --dir vendor/vep -i {input.vcf:q} --fasta {input.vep_refseq:q} --everything \
                --plugin dbNSFP,{input.dbnsfp_db:q},SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_pred,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,REVEL_score,REVEL_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore \
                --plugin SpliceAI,snv={input.spliceai_snv:q},indel={input.spliceai_indel:q} \
                --plugin GO,dir={input.go_terms:q} --fork 4 --force_overwrite --vcf -o {output:q}"

rule split_vep_annotated_variant_files:
    input:
        vcf="output/vep-annotated-variants/annotated.vcf.gz",
        tbi="output/vep-annotated-variants/annotated.vcf.gz.tbi"
    output: "output/vep-annotated-variants/annotated.{chrN}.vcf.gz"
    shell: "bcftools view -Oz {input.vcf:q} -o {output:q} --regions {wildcards.chrN}"

rule download_gnomad_exome_files:
    output: "vendor/gnomad/{version}/exomes/{file}"
    shell: "curl -k 'https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.version}/vcf/exomes/{wildcards.file}' -o {output:q}"

rule download_gnomad_genome_files:
    output: "vendor/gnomad/{version}/genomes/{file}"
    shell: "curl -k 'https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.version}/vcf/genomes/{wildcards.file}' -o {output:q}"

rule download_gnomad_sv_files:
   output: "vendor/gnomad/{version}/genome_sv/{file}"
   shell: "curl -k 'https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.version}/genome_sv/{wildcards.file}' -o {output:q}"

rule download_clinvar_files:
    output: "vendor/clinvar/{build}/{file}"
    shell: "curl -k 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{wildcards.build}/{wildcards.file}' -o {output:q}"

rule download_dbsnp_files:
    output: "vendor/dbsnp/GRCh38/{file}"
    shell: "curl -k 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/{wildcards.file}' -o {output:q}"

rule interpolate_vcfanno_gnomad_file_nuclear:
    input:
        config="config/gnomad-chrN.conf.in",
        gnomad_vcf="vendor/gnomad/4.1/genomes/gnomad.genomes.v4.1.sites.{chrN}.vcf.bgz"
    output: "output/gnomad-annotated-variants/vcfanno_gnomad_nsv.{chrN}.conf"
    shell: "GNOMAD_VCF_FILE={input.gnomad_vcf:q} envsubst < {input.config:q} > {output:q}"

rule interpolate_vcfanno_gnomad_file_mitochondrial:
    input:
        config="config/gnomad-chrM.conf.in",
        gnomad_vcf="vendor/gnomad/3.1/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
    output: "output/gnomad-annotated-variants/vcfanno_gnomad_nsv.chrM.conf"
    shell: "GNOMAD_VCF_FILE={input.gnomad_vcf:q} envsubst < {input.config:q} > {output:q}"

rule interpolate_vcfanno_gnomad_file_sv:
   input:
       config="config/gnomad-sv.conf.in",
       gnomad_vcf="vendor/gnomad/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz"
   output: "output/gnomad-annotated-variants/vcfanno_gnomad_sv.conf"
   shell: "GNOMAD_VCF_FILE={input.gnomad_vcf:q} envsubst < {input.config:q} > {output:q}"

rule gnomad_annotate_nuclear_variants:
    input:
        gnomad_vcf="vendor/gnomad/4.1/genomes/gnomad.genomes.v4.1.sites.{chrN}.vcf.bgz",
        gnomad_tbi="vendor/gnomad/4.1/genomes/gnomad.genomes.v4.1.sites.{chrN}.vcf.bgz.tbi",
        vcfanno_config="output/gnomad-annotated-variants/vcfanno_gnomad_nsv.{chrN}.conf",
        vcf="output/vep-annotated-variants/annotated.{chrN}.vcf.gz"
    output: "output/gnomad-annotated-variants/annotated_nsv.{chrN}.vcf"
    shell: "vcfanno -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule gnomad_annotate_mitochondrial_variants:
    input:
        gnomad_vcf="vendor/gnomad/3.1/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
        gnomad_tbi="vendor/gnomad/3.1/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi",
        vcfanno_config="output/gnomad-annotated-variants/vcfanno_gnomad_nsv.chrM.conf",
        vcf="output/vep-annotated-variants/annotated_go.chrM.vcf.gz"
    output: "output/gnomad-annotated-variants/annotated_nsv.chrM.vcf"
    shell: "vcfanno -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule gnomad_annotate_structural_variants:
    input:
       gnomad_sv_vcf="vendor/gnomad/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz",
       gnomad_sv_tbi="vendor/gnomad/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi",
       vcfanno_config="output/gnomad-annotated-variants/vcfanno_gnomad_sv.conf",
       vcf="output/gnomad-annotated-variants/annotated_nsv.{chrN}.vcf"
    output: "output/gnomad-annotated-variants/annotated.{chrN}.vcf"
    shell: "vcfanno -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule generate_chr_conversion_file:
    output: "output/clinvar-annotated-variants/chr-conv.txt"
    run:
        with open(output[0], "w") as out:
            for n in all_chrN:
                if n == "M":
                    print(f'chr{n}\tMT', file=out)
                else:
                    print(f'chr{n}\t{n}', file=out)

rule rename_chrs_in_gnomad_annotate_variants:
    input:
        vcf="output/gnomad-annotated-variants/annotated.{chrN}.vcf",
        chr_conv="output/clinvar-annotated-variants/chr-conv.txt"
    output: "output/clinvar-annotated-variants/chr_renamed.{chrN}.vcf"
    shell: "bcftools annotate --rename-chrs {input.chr_conv:q} {input.vcf:q} -Ov -o {output:q}"

rule interpolate_vcfanno_clinvar_dbsnp_file:
    input:
        config="config/clinvar_dbsnp.conf.in",
        clinvar_vcf="vendor/clinvar/GRCh38/clinvar_20240221.vcf.gz",
        dbsnp_vcf="vendor/dbsnp/GRCh38/All_20180418.vcf.gz"
    output: "output/clinvar-annotated-variants/vcfanno_clinvar.conf"
    shell: "CLINVAR_VCF_FILE={input.clinvar_vcf:q} DBSNP_VCF_FILE={input.dbsnp_vcf:q} envsubst < {input.config:q} > {output:q}"

rule clinvar_annotate_variant_files:
    input:
        clinvar_vcf="vendor/clinvar/GRCh38/clinvar_20240221.vcf.gz",
        clinvar_tbi="vendor/clinvar/GRCh38/clinvar_20240221.vcf.gz.tbi",
        dbsnp_vcf="vendor/dbsnp/GRCh38/All_20180418.vcf.gz",
        dbsnp_tbi="vendor/dbsnp/GRCh38/All_20180418.vcf.gz.tbi",
        vcfanno_config="output/clinvar-annotated-variants/vcfanno_clinvar.conf",
        vcf="output/clinvar-annotated-variants/chr_renamed.{chrN}.vcf"
    output: "output/clinvar-annotated-variants/annotated.{chrN}.vcf"
    shell: "vcfanno -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule concat_clinvar_annotated_variant_files:
    input: expand("output/clinvar-annotated-variants/annotated.chr{N}.vcf", N=all_chrN)
    output: "output/clinvar-annotated-variants/annotated.bcf"
    shell: "bcftools concat {input:q} -Ob -o {output:q}"

rule split_vep_tags_in_variant_files:
    input: "output/clinvar-annotated-variants/annotated.bcf"
    output: "output/split-tag-variants/split.bcf"
    shell: "bcftools +split-vep -c IMPACT,Polyphen2_HVAR_pred,SIFT_pred,CADD_phred:Float,SIFT4G_pred,GO -s worst {input:q} -Ob -o {output:q}"

rule fill_tags_in_variant_files:
    input: "output/split-tag-variants/split.bcf"""
    output: "output/split-tag-variants/annotated.bcf"
    shell: "bcftools +fill-tags {input:q} -Ob -o {output:q}"

rule export_variant_stats:
    input: "output/clinvar-annotated-variants/annotated.bcf"
    output: "output/analysis-report/variant-stats.tsv"
    shell: "bcftools +split-vep -H -d \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%INFO/gno_id\t%INFO/dbsnp_id\t%INFO/clinvar_id\t%INFO/gno_af_all\t%CSQ\n' -A tab \
        {input:q} | sort | uniq > {output:q}"

rule filter_annotated_variants_maf_lt:
    input: "output/tag-split-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/maf_lt_{threshold}.vcf.gz"
    shell: "bcftools view -i 'MAF < {wildcards.threshold}' -Oz -o {output:q} {input:q}"

rule filter_annotated_variants_maf_gt:
    input: "output/split-tag-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/maf_gt_{threshold}.vcf.gz"
    shell: "bcftools view -i 'MAF > {wildcards.threshold}' -Oz -o {output:q} {input:q}"

rule filter_annotated_variants_by_impact:
    input: "output/split-tag-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/impact_{category}.vcf.gz"
    shell: "bcftools view -i 'INFO/IMPACT ~ \"{wildcards.category}\"' -Oz -o {output:q} {input:q}"

rule filter_annotated_variants_by_sift:
    input: "output/split-tag-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/sift_{category}.vcf.gz"
    shell: "bcftools view -i 'INFO/SIFT_pred ~ \"{wildcards.category}\"' -Oz -o {output:q} {input:q}"

rule filter_annotated_variants_by_polyphen2:
    input: "output/split-tag-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/polyphen2_{category}.vcf.gz"
    shell: "bcftools view -i 'INFO/Polyphen2_HVAR_pred ~ \"{wildcards.category}\"' -Oz -o {output:q} {input:q}"

rule filter_annotated_variants_by_cadd:
    input: "output/split-tag-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/cadd_{threshold}.vcf.gz"
    shell: "bcftools view -i 'INFO/CADD_phred >= {wildcards.threshold}' -Oz -o {output:q} {input:q}"

rule filter_annotated_variants_by_sift4g:
    input: "output/split-tag-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/sift4g_{category}.vcf.gz"
    shell: "bcftools view -i 'INFO/SIFT4G_pred ~ \"{wildcards.category}\"' -Oz -o {output:q} {input:q}"

rule extract_sample_information_for_ms_subjects:
    input:
        sample_names="output/renamed-variant-samples/samples.normalized.txt",
        rscript="src/extract-ms-samples-info.r"
    output: "output/phenotype-info/MS.tsv"
    shell: "Rscript {input.rscript:q} {input.sample_names:q} > {output:q}"

rule extract_sample_information_for_als_subjects:
    input:
        sample_names="output/renamed-variant-samples/samples.normalized.txt",
        rscript="src/extract-als-samples-info.r"
    output: "output/phenotype-info/ALS.tsv"
    shell: "Rscript {input.rscript:q} {input.sample_names:q} > {output:q}"

rule generate_overall_sample_information:
    input:
        als_info="output/phenotype-info/ALS.tsv",
        ms_info="output/phenotype-info/MS.tsv"
    output: "output/phenotype-info/ALL.tsv"
    shell: "echo -e 'SAMPLE\tNHC\tSEX\tPHENO' > {output:q} && \
            cut -f 1-4 <(tail -n +2 {input.als_info:q}) >> {output:q} && \
            cut -f 1-4 <(tail -n +2 {input.ms_info:q})  >> {output:q}"

rule extract_samples_by_pheno:
    input: "output/phenotype-info/ALL.tsv"
    output: "output/phenotype-info/{pheno}.samples"
    shell: "tail -n +2 {input:q} | gawk '$4 == \"{wildcards.pheno}\" {{ print $1 }}' > {output:q}"

rule generate_psam_file:
    input: "output/phenotype-info/ALL.tsv"
    output: "output/variant-analysis/input/samples.psam"
    shell: "echo -e '#IID\tSEX\tALS\tMS\tALS_MS' > {output:q} && \
            gawk ' \
                BEGIN {{ OFS=\"\\t\"; }}\
                {{ \
                    print $1, \
                        $3 == \"M\" ? 1 \
                            : $3 == \"F\" ? 2 \
                            : \"NA\", \
                        $4 == \"CONTROL\" ? 1 \
                            : $4 == \"ALS\" ? 2 \
                            : \"NA\", \
                        $4 == \"CONTROL\" ? 1 \
                            : $4 == \"MS\" ? 2 \
                            : \"NA\", \
                        $4 == \"CONTROL\" ? 1 \
                            : ($4 == \"ALS\" || $4 == \"MS\") ? 2 \
                            : \"NA\" \
                }} \
            ' <(tail -n +2 {input:q}) >> {output:q}"

rule generate_ped_file:
    input: "output/phenotype-info/ALL.tsv"
    output: "output/variant-analysis/input/samples.ped"
    shell: "echo -e 'FID\tIID\tFATID\tMATID\tSEX\tALS\tMS\tALS_MS' > {output:q} && \
            gawk ' \
                BEGIN {{ OFS=\"\\t\"; }}\
                {{ \
                    print $1, $1, 0, 0, \
                        $3 == \"M\" ? 1 \
                            : $3 == \"F\" ? 2 \
                            : \"NA\", \
                        $4 == \"CONTROL\" ? 1 \
                            : $4 == \"ALS\" ? 2 \
                            : \"NA\", \
                        $4 == \"CONTROL\" ? 1 \
                            : $4 == \"MS\" ? 2 \
                            : \"NA\", \
                        $4 == \"CONTROL\" ? 1 \
                            : ($4 == \"ALS\" || $4 == \"MS\") ? 2 \
                            : \"NA\" \
                }} \
            ' <(tail -n +2 {input:q}) >> {output:q}"

rule convert_unfiltered_variants_to_plink:
    input:
        vcf="output/clinvar-annotated-variants/annotated.vcf.gz",
        psam="output/variant-analysis/input/samples.psam"
    output: multiext("output/variant-analysis/input/samples", ".pgen", ".pvar")
    shell: "plink2 --make-pgen --vcf {input.vcf:q} --psam {input.psam:q} --double-id \
                   --mind 0.1 --geno 0.05 --hwe 1e-8 --snps-only \
                   --set-all-var-ids @:#_\\$r_\\$a \
                   --out output/variant-analysis/input/plink2 && \
            mv output/variant-analysis/input/plink2.pgen output/variant-analysis/input/samples.pgen && \
            mv output/variant-analysis/input/plink2.pvar output/variant-analysis/input/samples.pvar && \
            rm output/variant-analysis/input/plink2.psam"

rule extract_variants_pca_for_samples:
    input: "output/variant-analysis/input/samples.pgen"
    output: multiext("output/variant-analysis/ALL_PCA", ".eigenvec", ".eigenval")
    shell: "plink2 --pfile output/variant-analysis/input/samples \
                   --pca --maf 0.05 --indep-pairwise 200 50 0.2 \
                   --out output/variant-analysis/ALL_PCA"

rule make_pca_scree_plot:
    input:
        eigenval="output/{path}.eigenval",
        rscript="src/make-scree-plot.r"
    output: "output/{path}.eigenval.scree.png"
    shell: "Rscript {input.rscript:q} {input.eigenval:q} {output:q}"

rule calculate_hwe_for_controls:
    input: multiext("output/variant-analysis/input/samples", ".pgen", ".pvar", ".psam")
    output: "output/variant-analysis/CONTROL_HWE.hardy"
    shell: "plink2 --pfile output/variant-analysis/input/samples \
                   --hardy --snps-only --autosome --require-pheno ALS MS \
                   --out output/variant-analysis/CONTROL_HWE"

rule make_hwe_ternary_plot:
    input:
        hwe="output/{path}.hardy",
        rscript="src/make-ternary-plot.r"
    output: "output/{path}.hardy.ternary.png"
    shell: "Rscript {input.rscript:q} {input.hwe:q} {output:q}"

rule test_common_variants_association:
    input:
        pbinary_set=multiext("output/variant-analysis/input/samples", ".pgen", ".pvar", ".psam"),
        samples_pca="output/variant-analysis/ALL_PCA.eigenvec"
    output: expand("output/variant-analysis/cva/results.{pheno}.glm.logistic.hybrid", pheno=all_pheno)
    shell: "plink2 --pfile output/variant-analysis/input/samples --glm --maf 0.05 \
                   --covar {input.samples_pca:q} --covar-col-nums 2-5 \
                   --out output/variant-analysis/cva/results"

rule download_ucsc_ref_flat_genes_hg38:
    output: "vendor/rvtests/refFlat_hg38.txt.gz"
    shell: "wget -O - 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz' > {output:q}"

rule test_rare_variants_association_nofilter_maf_any:
    input:
        vcf="output/clinvar-annotated-variants/annotated.vcf.gz",
        tbi="output/clinvar-annotated-variants/annotated.vcf.gz.tbi",
        ped="output/variant-analysis/input/samples.ped",
        samples_pca="output/variant-analysis/ALL_PCA.eigenvec",
        ucsc_gene_file="vendor/rvtests/refFlat_hg38.txt.gz"
    output: multiext("output/variant-analysis/rva.nofilter.maf_any/{pheno}", \
                     ".CMC.assoc", ".Skat.assoc", ".SkatO.assoc", ".VariableThresholdPrice.assoc")
    shell: "rvtest --noweb --inVcf {input.vcf:q} \
                   --pheno {input.ped:q} --pheno-name {wildcards.pheno} \
                   --geneFile {input.ucsc_gene_file:q} \
                   --burden cmc --vt price --kernel skat,skato \
                   --out output/variant-analysis/rva.nofilter.maf_any/{wildcards.pheno}"

ruleorder:
    test_rare_variants_association_nofilter_maf_any >
    test_rare_variants_association_nofilter_maf_threshold >
    test_rare_variants_association_filtered_maf_any >
    test_rare_variants_association_filtered_maf_threshold

rule test_rare_variants_association_filtered_maf_any:
    input:
        vcf="output/filtered-annotated-variants/{filter}.vcf.gz",
        tbi="output/filtered-annotated-variants/{filter}.vcf.gz.tbi",
        ped="output/variant-analysis/input/samples.ped",
        ucsc_gene_file="vendor/rvtests/refFlat_hg38.txt.gz"
    output: multiext("output/variant-analysis/rva.{filter}.maf_any/{pheno}", \
                     ".CMC.assoc", ".Skat.assoc", ".SkatO.assoc", ".VariableThresholdPrice.assoc")
    shell: "rvtest --noweb --inVcf {input.vcf:q} \
                   --pheno {input.ped:q} --pheno-name {wildcards.pheno} \
                   --geneFile {input.ucsc_gene_file:q} \
                   --burden cmc --vt price --kernel skat,skato \
                   --out output/variant-analysis/rva.{wildcards.filter}.maf_any/{wildcards.pheno}"

rule test_rare_variants_association_nofilter_maf_threshold:
    input:
        vcf="output/clinvar-annotated-variants/annotated.vcf.gz",
        tbi="output/clinvar-annotated-variants/annotated.vcf.gz.tbi",
        ped="output/variant-analysis/input/samples.ped",
        ucsc_gene_file="vendor/rvtests/refFlat_hg38.txt.gz"
    output: multiext("output/variant-analysis/rva.nofilter.maf_{threshold,[0-9.]+}/{pheno}", \
                     ".CMC.assoc", ".Skat.assoc", ".SkatO.assoc", ".VariableThresholdPrice.assoc")
    shell: "rvtest --noweb --inVcf {input.vcf:q} \
                   --pheno {input.ped:q} --pheno-name {wildcards.pheno} \
                   --freqUpper {wildcards.threshold} \
                   --geneFile {input.ucsc_gene_file:q} \
                   --burden cmc --vt price --kernel skat,skato \
                   --out output/variant-analysis/rva.nofilter.maf_{wildcards.threshold}/{wildcards.pheno}"

rule test_rare_variants_association_filtered_maf_threshold:
    input:
        vcf="output/filtered-annotated-variants/{filter}.vcf.gz",
        tbi="output/filtered-annotated-variants/{filter}.vcf.gz.tbi",
        ped="output/variant-analysis/input/samples.ped",
        ucsc_gene_file="vendor/rvtests/refFlat_hg38.txt.gz"
    output: multiext("output/variant-analysis/rva.{filter}.maf_{threshold,[0-9.]+}/{pheno}", \
                     ".CMC.assoc", ".Skat.assoc", ".SkatO.assoc", ".VariableThresholdPrice.assoc")
    shell: "rvtest --noweb --inVcf {input.vcf:q} \
                   --pheno {input.ped:q} --pheno-name {wildcards.pheno} \
                   --freqUpper {wildcards.threshold} \
                   --geneFile {input.ucsc_gene_file:q} \
                   --burden cmc --vt price --kernel skat,skato \
                   --out output/variant-analysis/rva.{wildcards.filter}.maf_{wildcards.threshold}/{wildcards.pheno}"

rule summarize_association_test_results_cva:
    input: "output/variant-analysis/cva/{path}"
    output: "output/variant-analysis/cva/{path}.summary"
    shell: "echo -e 'NAME\tP' > {output:q} && cut -d $'\t' -f 3,16 <(tail -n +2 {input:q} | grep ADD) >> {output:q}"

rule summarize_association_test_results_rva_cmc:
    input: "output/variant-analysis/rva.{path}.CMC.assoc"
    output: "output/variant-analysis/rva.{path}.CMC.assoc.summary"
    shell: "echo -e 'NAME\tP' > {output:q} && cut -d $'\t' -f 1,7 <(tail -n +2 {input:q}) >> {output:q}"

rule summarize_association_test_results_rva_skat:
    input: "output/variant-analysis/rva.{path}.Skat.assoc"
    output: "output/variant-analysis/rva.{path}.Skat.assoc.summary"
    shell: "echo -e 'NAME\tP' > {output:q} && cut -d $'\t' -f 1,13 <(tail -n +2 {input:q}) >> {output:q}"

rule summarize_association_test_results_rva_skato:
    input: "output/variant-analysis/rva.{path}.SkatO.assoc"
    output: "output/variant-analysis/rva.{path}.SkatO.assoc.summary"
    shell: "echo -e 'NAME\tP' > {output:q} && cut -d $'\t' -f 1,8 <(tail -n +2 {input:q}) >> {output:q}"

rule summarize_association_test_results_rva_vt:
    input: "output/variant-analysis/rva.{path}.VariableThresholdPrice.assoc"
    output: "output/variant-analysis/rva.{path}.VariableThresholdPrice.assoc.summary"
    shell: "echo -e 'NAME\tP' > {output:q} && cut -d $'\t' -f 1,14 <(tail -n +2 {input:q}) >> {output:q}"

rule make_volcano_plot_from_association_test_results:
    input:
        rscript="src/make-volcano-plot.r",
        results="output/variant-analysis/{path}"
    output: "output/variant-analysis/{path}.volcano.png"
    shell: "Rscript {input.rscript:q} {input.results:q} {output:q}"

rule make_manhattan_plot_from_association_test_results:
    input:
        rscript="src/make-manhattan-plot.r",
        results="output/variant-analysis/{path}"
    output: "output/variant-analysis/{path}.manhattan.png"
    shell: "Rscript {input.rscript:q} {input.results:q} {output:q}"

rule make_qq_plot_from_association_test_results:
    input:
        rscript="src/make-qq-plot.r",
        res_summary="output/variant-analysis/{path}.summary"
    output: "output/variant-analysis/{path}.qq.png"
    shell: "Rscript {input.rscript:q} {input.res_summary:q} {output:q}"

rule make_bar_plot_from_association_test_results:
    input:
        rscript="src/make-bar-plot.r",
        res_summary="output/variant-analysis/{path}.summary"
    output: "output/variant-analysis/{path}.bar.png"
    shell: "Rscript {input.rscript:q} {input.res_summary:q} {output:q}"