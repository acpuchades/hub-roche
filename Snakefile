import re

SPLICEAI_GENOME_SCORES_DIR = "vendor/spliceai/genome_scores_v1.3-194103939/genome_scores_v1.3-ds.20a701bc58ab45b59de2576db79ac8d0"

samples, = glob_wildcards("data/laboratorio/vcf/{sample}.vcf.gz")
fastq_reads, = glob_wildcards("data/laboratorio/fastq/{read}.fastq.gz")
fastq_names = list(set(re.sub('_R[12]_001', '', x) for x in fastq_reads))

auto_chrN = list(map(str, range(1, 23)))
nuclear_chrN = auto_chrN + ['X', 'Y']
all_chrN = nuclear_chrN + ['M']

rvtests = ["Skat", "VariableThresholdPrice", "CMC", "Kbac"]

rule all:
    input:
        "output/variant-analysis/controls.hardy",
        "output/variant-analysis/cva.ALS.glm.logistic.hybrid.volcano.png",
        "output/variant-analysis/cva.ALS.glm.logistic.hybrid.manhattan.png",
        "output/variant-analysis/cva.MS.glm.logistic.hybrid.volcano.png",
        "output/variant-analysis/cva.MS.glm.logistic.hybrid.manhattan.png",
        expand("output/variant-analysis/rva.ALS.{t}.assoc", t=rvtests),
        expand("output/variant-analysis/rva.MS.{t}.assoc", t=rvtests)

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

rule sort_vcf_file:
    input: "data/laboratorio/vcf/{sample}.vcf.gz"
    output: "output/sorted-variants/{sample}.vcf.gz"
    shell: "bcftools sort -Oz {input:q} -o {output:q}"

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

rule remove_format_pl_from_vcf_file:
    input: "output/normalized-variants/{sample}.vcf.gz"
    output: "output/noformatpl-variants/{sample}.vcf.gz"
    shell: "bcftools annotate -x FORMAT/PL -Oz {input:q} -o {output:q}"

rule merge_vcf_files:
    input:
        vcf=expand("output/noformatpl-variants/{sample}.vcf.gz", sample=samples),
        tbi=expand("output/noformatpl-variants/{sample}.vcf.gz.tbi", sample=samples)
    output: "output/merged-variants/merged.vcf.gz"
    shell: "bcftools merge -m none {input.vcf:q} -Oz -o {output:q}"

rule calculate_maf_from_vcf_file:
    input: "output/merged-variants/merged.vcf.gz"
    output: "output/tagged-variants/variants.bcf"
    shell: "bcftools +fill-tags {input:q} -Ob -o {output:q}"

rule extract_variant_file_sample_names:
    input: "output/tagged-variants/variants.bcf"
    output: "output/renamed-variants/samples.original.txt"
    shell: "bcftools query --list-samples {input:q} > {output:q}"

rule normalize_sample_names:
    input: "output/renamed-variants/samples.original.txt"
    output: "output/renamed-variants/samples.normalized.txt"
    shell: "sed -e 's/-ambBED$//' -e 's/-//g' \
        -e 's/ELA1912ADN1/ELA1912-2/' -e 's/ELA1923ADN1/ELA1923-2/' \
        -e 's/ELA1924ADN1/ELA1924-2/' -e 's/ADN[12]$//' {input:q} > {output:q}"

rule normalize_vcf_sample_headers:
    input:
        vcf="output/tagged-variants/variants.bcf",
        normalized_names="output/renamed-variants/samples.normalized.txt"
    output: "output/renamed-variants/normalized.vcf"
    shell: "bcftools reheader -s {input.normalized_names} {input.vcf:q} | gunzip > {output:q}"

rule sort_sample_names:
    input: "output/renamed-variants/samples.normalized.txt"
    output: "output/renamed-variants/samples.sorted.txt"
    shell: "sort -k1,1 {input:q} > {output:q}"

rule sort_vcf_sample_headers:
    input:
        vcf="output/renamed-variants/normalized.vcf",
        sorted_names="output/renamed-variants/samples.sorted.txt"
    output: "output/renamed-variants/renamed.vcf"
    shell: "bcftools view -S {input.sorted_names} -Ov -o {output:q} {input.vcf:q}"

rule filter_vcf_files_in:
    input: "output/renamed-variants/renamed.vcf"
    output: "output/filtered-variants/included.vcf"
    shell: "bcftools view -f PASS {input:q} -Ou -o {output:q}"

rule filter_vcf_files_out:
    input: "output/renamed-variants/renamed.vcf"
    output: "output/filtered-variants/excluded.vcf"
    shell: "bcftools view -f FAIL {input:q} -Ou -o {output:q}"

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
        vcf="output/filtered-variants/included.vcf"
    output: "output/vep-annotated-variants/annotated.vcf"
    shell: "vep --cache --dir vendor/vep -i {input.vcf:q} --fasta {input.vep_refseq:q} --everything \
        --plugin dbNSFP,{input.dbnsfp_db:q},SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,REVEL_score,REVEL_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore \
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

rule export_variant_stats:
    input: "output/clinvar-annotated-variants/annotated.bcf"
    output: "output/analysis-report/variant-stats.tsv"
    shell: "bcftools +split-vep -H -d \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%INFO/gno_id\t%INFO/dbsnp_id\t%INFO/clinvar_id\t%INFO/gno_af_all\t%CSQ\n' -A tab \
        {input:q} | sort | uniq > {output:q}"

rule filter_annotated_variants_maf_lt:
    input:  "output/clinvar-annotated-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/maf_lt_{threshold}.bcf"
    shell: "bcftools view -i 'MAF < {wildcards.threshold}' -Ob -o {output:q} {input:q}"

rule filter_annotated_variants_maf_gt:
    input: "output/clinvar-annotated-variants/annotated.bcf"
    output: "output/filtered-annotated-variants/maf_gt_{threshold}.bcf"
    shell: "bcftools view -i 'MAF > {wildcards.threshold}' -Ob -o {output:q} {input:q}"

rule generate_sample_information_file:
    input: "output/renamed-variants/samples.normalized.txt"
    output: "output/phenotype-info/samples.tsv"
    shell: "Rscript src/make-sample-info.r > {output:q}"

rule extract_control_samples:
    input: "output/phenotype-info/samples.tsv"
    output: "output/phenotype-info/control.samples"
    shell: "gawk '$3 == \"CONTROL\" {{ print $1 }}' {input:q} | tail -n +2 > {output:q}"

rule generate_psam_file:
    input: "output/phenotype-info/samples.tsv"
    output: "output/phenotype-info/samples.psam"
    shell: "tail -n +2 {input:q} | gawk ' \
                BEGIN {{ \
                    OFS=\"\\t\"; \
                    print \"#IID\", \"SEX\", \"ALS\", \"MS\" \
                }} \
                {{ \
                    print $1, $2, \
                        $3 == \"CONTROL\" ? 1 : $3 == \"ALS\" ? 2 : -9, \
                        $3 == \"CONTROL\" ? 1 : $3 == \"MS\"  ? 2 : -9  \
                }} \
            ' > {output:q}"

rule generate_ped_file:
    input: "output/phenotype-info/samples.tsv"
    output: "output/phenotype-info/samples.ped"
    shell: "tail -n +2 {input:q} | gawk ' \
                BEGIN {{ \
                    OFS=\"\\t\"; \
                    print \"FID\", \"IID\", \"FATID\", \"MATID\", \"SEX\", \"ALS\", \"MS\" \
                }} \
                {{ \
                    print $1, $1, 0, 0, \
                        $2 == \"M\" ? 1 : $2 == \"F\" ? 2 : 0, \
                        $3 == \"CONTROL\" ? 1 : $3 == \"ALS\" ? 2 : 0, \
                        $3 == \"CONTROL\" ? 1 : $3 == \"MS\"  ? 2 : 0  \
                }} \
            ' > {output:q}"

rule convert_unfiltered_variants_to_plink:
    input:
        bcf="output/clinvar-annotated-variants/annotated.bcf",
        psam="output/phenotype-info/samples.psam"
    output: multiext("output/variant-analysis/unfiltered", ".pgen", ".pvar", ".psam")
    shell: "plink2 --make-pgen --bcf {input.bcf:q} \
                   --psam {input.psam:q} --double-id  \
                   --out output/variant-analysis/unfiltered"

rule test_hwe_on_unfiltered_variants_in_controls:
    input: multiext("output/variant-analysis/unfiltered", ".pgen", ".pvar", ".psam")
    output: "output/variant-analysis/controls.hardy"
    shell: "plink2 --pfile output/variant-analysis/unfiltered \
                   --hardy --snps-only --autosome \
                   --keep-if ALS == 1 --keep-founders \
                   --out output/variant-analysis/controls"

rule test_common_variants_association:
    input: multiext("output/variant-analysis/unfiltered", ".pgen", ".pvar", ".psam")
    output: expand("output/variant-analysis/cva.{c}.glm.logistic.hybrid", c=["ALS", "MS"])
    shell: "plink2 --pfile output/variant-analysis/unfiltered --maf 0.05 \
                   --glm allow-no-covars --out output/variant-analysis/cva"

rule download_ucsc_ref_flat_genes_hg19:
    output: "vendor/rvtests/refFlat_hg19.txt.gz"
    shell: "wget -O - 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz' > {output:q}"

rule test_rare_variants_association:
    input:
        vcf="output/clinvar-annotated-variants/annotated.vcf.gz",
        tbi="output/clinvar-annotated-variants/annotated.vcf.gz.tbi",
        ped="output/phenotype-info/samples.ped",
        ucsc_gene_file="vendor/rvtests/refFlat_hg19.txt.gz"
    output: multiext("output/variant-analysis/rva.{pheno}", \
                     ".CMC.assoc", ".Kbac.assoc", ".Skat.assoc", \
                     ".VariableThresholdPrice.assoc")
    shell: "rvtest --noweb --inVcf {input.vcf:q} \
                   --pheno {input.ped:q} --pheno-name {wildcards.pheno} \
                   --freqUpper 0.05 --burden cmc --vt price --kernel skat,kbac \
                   --geneFile {input.ucsc_gene_file:q} \
                   --out output/variant-analysis/rva.{wildcards.pheno}"

rule make_volcano_plot_from_association_test_results:
    input:
        rscript="src/make-volcano-plot.r",
        results="output/variant-analysis/{results}"
    output: "output/variant-analysis/{results}.volcano.png"
    shell: "Rscript {input.rscript:q} {input.results:q} {output:q}"

rule make_manhattan_plot_from_association_test_results:
    input:
        rscript="src/make-manhattan-plot.r",
        results="output/variant-analysis/{results}"
    output: "output/variant-analysis/{results}.manhattan.png"
    shell: "Rscript {input.rscript:q} {input.results:q} {output:q}"