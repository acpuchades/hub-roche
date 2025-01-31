import re

SPLICEAI_GENOME_SCORES_DIR = "vendor/spliceai/genome_scores_v1.3-194103939/genome_scores_v1.3-ds.20a701bc58ab45b59de2576db79ac8d0"

samples, = glob_wildcards("data/laboratorio/vcf/{sample}.vcf.gz")
fastq_reads, = glob_wildcards("data/laboratorio/fastq/{read}.fastq.gz")
fastq_names = list(set(re.sub('_R[12]_001', '', x) for x in fastq_reads))

auto_chrN = list(map(str, range(1, 23)))
nuclear_chrN = auto_chrN + ['X', 'Y']
all_chrN = nuclear_chrN + ['M']

sys.path.insert(1, './src')

rule all:
    input:
        "output/merged-variants/merged.vcf.gz",
        "output/filtered-variants/included.vcf",
        "output/filtered-variants/excluded.vcf",
        "output/vep-annotated-variants/merged.vcf",
        "output/gnomad-annotated-variants/annotated-nsv.vcf",
        "output/clinvar-annotated-variants/annotated.vcf",
        "output/analysis-report/pipeline-summary.png",
        "output/analysis-report/variants-coverage.png",
        "output/analysis-report/variants-coverage.png",
        "output/analysis-report/variants-quality.png",
        "output/analysis-report/variants-reported.png",
        "output/analysis-report/gnomad-af.png",
        "output/common-variants/plink.eigenvec",
        "output/common-variants/plink.eigenval",
        "output/common-variants/variants-pca.html",
        expand("output/fastq-qc-plots/{read}_fastqc.html", read=fastq_reads),
        expand("output/trimmed-fastq-reads/{read}_{rn}{up}.fastq.gz", read=fastq_names, rn=[1,2], up=["U", "P"]),
        expand("output/trimmed-fastq-qc-plots/{read}_{rn}P_fastqc.html", read=fastq_names, rn=[1,2])

rule clean:
    shell: "rm -rf output"

rule compress_output_file_bgzip:
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

rule unpack_annovar:
    input: "vendor/annovar.latest.tar.gz"
    output: directory("vendor/annovar")
    shell: "tar -xvf {input} -C vendor"

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
        vcf = expand("output/noformatpl-variants/{sample}.vcf.gz", sample=samples),
        tbi = expand("output/noformatpl-variants/{sample}.vcf.gz.tbi", sample=samples)
    output: "output/merged-variants/merged.vcf.gz"
    shell: "bcftools merge -m none {input.vcf:q} -Oz -o {output:q}"

rule extract_variant_file_sample_names:
    input: "output/merged-variants/merged.vcf.gz"
    output: "output/renamed-variants/samples.original.txt"
    shell: "bcftools query --list-samples {input:q} > {output:q}"

rule sort_vcf_sample_names:
    input: "output/renamed-variants/samples.original.txt"
    output: "output/renamed-variants/samples.sorted.txt"
    shell: "sort -k1,1 {input:q} > {output:q}"

rule sort_renamed_variant_file_samples:
    input:
        vcf="output/merged-variants/merged.vcf.gz",
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
        vcf="output/renamed-variants/renamed.vcf"
    output: multiext("output/annovar-annotated-variants/annotated.hg38_multianno", ".vcf", ".txt")
    shell: "{input.av_tableav} -buildver hg38 {input.vcf:q} \
            -out output/annovar-annotated-variants/annotated \
            -protocol refGene,cytoBand,dbnsfp30a \
            -operation g,r,f -nastring . -vcfinput -polish \
            vendor/annovar/humandb"

rule download_vep_cache_files:
    output: "vendor/vep/homo_sapiens_vep_{ver}_{genome}.tar.gz"
    shell: "curl https://ftp.ensembl.org/pub/release-{wildcards.ver}/variation/indexed_vep_cache/homo_sapiens_vep_{wildcards.ver}_{wildcards.genome}.tar.gz -o {output:q}"

rule download_vep_refseq_files:
    output: "vendor/vep/homo_sapiens_vep_{ver}_GRCh38.dna.primary_assembly.fa.gz"
    shell: "curl https://ftp.ensembl.org/pub/release-{wildcards.ver}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -o {output:q}"

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
        vep_cache="vendor/vep/homo_sapiens/112_GRCh38",
        vep_refseq="vendor/vep/homo_sapiens_vep_112_GRCh38.dna.primary_assembly.fa.bgz",
        dbnsfp_db="vendor/vep/dbnsfp/dbNSFP4.8a_GRCh38.gz",
        dbnsfp_tbi="vendor/vep/dbnsfp/dbNSFP4.8a_GRCh38.gz.tbi",
        spliceai_snv=f"{SPLICEAI_GENOME_SCORES_DIR}/spliceai_scores.raw.snv.hg38.vcf.gz",
        spliceai_indel=f"{SPLICEAI_GENOME_SCORES_DIR}/spliceai_scores.raw.indel.hg38.vcf.gz",
        vcf="output/filtered-variants/included.vcf"
    output: "output/vep-annotated-variants/merged.vcf"
    shell: "vep --cache --dir vendor/vep -i {input.vcf:q} --fasta {input.vep_refseq:q} --everything \
        --plugin dbNSFP,{input.dbnsfp_db:q},SIFT_score,SIFT_converted_rankscore,SIFT_pred,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,REVEL_score,REVEL_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore \
        --plugin SpliceAI,snv={input.spliceai_snv:q},indel={input.spliceai_indel:q} \
        --fork 4 --offline --vcf -o {output:q}"

rule split_vep_annotated_variant_files:
    input:
        vcf="output/vep-annotated-variants/merged.vcf.gz",
        tbi="output/vep-annotated-variants/merged.vcf.gz.tbi"
    output: "output/vep-annotated-variants/split.{chrN}.vcf.gz"
    shell: "bcftools view -Oz {input.vcf:q} -o {output:q} --regions {wildcards.chrN}"

rule download_gnomad_exome_files:
    output: "vendor/gnomad/{version}/exomes/{file}"
    shell: "curl 'https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.version}/vcf/exomes/{wildcards.file}' -o {output:q}"

rule download_gnomad_genome_files:
    output: "vendor/gnomad/{version}/genomes/{file}"
    shell: "curl 'https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.version}/vcf/genomes/{wildcards.file}' -o {output:q}"

rule download_gnomad_sv_files:
   output: "vendor/gnomad/{version}/genome_sv/{file}"
   shell: "curl 'https://storage.googleapis.com/gcp-public-data--gnomad/release/{wildcards.version}/genome_sv/{wildcards.file}' -o {output:q}"

rule download_clinvar_files:
    output: "vendor/clinvar/{build}/{file}"
    shell: "curl 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{wildcards.build}/{wildcards.file}' -o {output:q}"

rule download_dbsnp_files:
    output: "vendor/dbsnp/GRCh38/{file}"
    shell: "curl 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/{wildcards.file}' -o {output:q}"

rule interpolate_vcfanno_gnomad_file_mitochondrial:
    input:
        config="config/gnomad-chrM.conf.in",
        gnomad_vcf="vendor/gnomad/3.1/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
    output: "output/gnomad-annotated-variants/vcfanno.chrM.conf"
    shell: "GNOMAD_VCF_FILE={input.gnomad_vcf:q} envsubst < {input.config:q} > {output:q}"

rule interpolate_vcfanno_gnomad_file_nuclear:
    input:
        config="config/gnomad.conf.in",
        gnomad_vcf="vendor/gnomad/4.1/genomes/gnomad.genomes.v4.1.sites.{chrN}.vcf.bgz"
    output: "output/gnomad-annotated-variants/vcfanno.{chrN}.conf"
    shell: "GNOMAD_VCF_FILE={input.gnomad_vcf:q} envsubst < {input.config:q} > {output:q}"

rule interpolate_vcfanno_gnomad_file_sv:
   input:
       config="config/gnomad-sv.conf.in",
       gnomad_vcf="vendor/gnomad/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz"
   output: "output/gnomad-annotated-variants/vcfanno.genome_sv.conf"
   shell: "GNOMAD_VCF_FILE={input.gnomad_vcf:q} envsubst < {input.config:q} > {output:q}"

rule gnomad_annotate_mitochondrial_variants:
    input:
        gnomad_vcf="vendor/gnomad/3.1/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz",
        gnomad_tbi="vendor/gnomad/3.1/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi",
        vcfanno_config="output/gnomad-annotated-variants/vcfanno.chrM.conf",
        vcf="output/vep-annotated-variants/split.chrM.vcf.gz"
    output: "output/gnomad-annotated-variants/annotated.chrM.vcf"
    shell: "vcfanno -p $(nproc --all) -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule gnomad_annotate_nuclear_variants:
    input:
        gnomad_vcf="vendor/gnomad/4.1/genomes/gnomad.genomes.v4.1.sites.{chrN}.vcf.bgz",
        gnomad_tbi="vendor/gnomad/4.1/genomes/gnomad.genomes.v4.1.sites.{chrN}.vcf.bgz.tbi",
        vcfanno_config="output/gnomad-annotated-variants/vcfanno.{chrN}.conf",
        vcf="output/vep-annotated-variants/split.{chrN}.vcf.gz"
    output: "output/gnomad-annotated-variants/annotated.{chrN}.vcf"
    shell: "vcfanno -p $(nproc --all) -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule concat_gnomad_annotated_variant_files:
    input:
        vcf=expand("output/gnomad-annotated-variants/annotated.chr{N}.vcf.gz", N=all_chrN),
        tbi=expand("output/gnomad-annotated-variants/annotated.chr{N}.vcf.gz.tbi", N=all_chrN)
    output: "output/gnomad-annotated-variants/annotated-nsv.vcf"
    shell: "bcftools concat {input.vcf:q} -Ov -o {output:q}"

rule gnomad_annotate_structural_variants:
    input:
       gnomad_sv_vcf="vendor/gnomad/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz",
       gnomad_sv_tbi="vendor/gnomad/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi",
       vcfanno_config="output/gnomad-annotated-variants/vcfanno.genome_sv.conf",
       vcf="output/gnomad-annotated-variants/annotated-nsv.vcf"
    output: "output/gnomad-annotated-variants/annotated.vcf"
    shell: "vcfanno -p $(nproc --all) -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule generate_chr_conversion_file:
    output: "output/gnomad-annotated-variants/chr-conv.txt"
    run:
        with open(output[0], "w") as out:
            for n in all_chrN:
                if n == "M":
                    print(f'chr{n}\tMT', file=out)
                else:
                    print(f'chr{n}\t{n}', file=out)

rule rename_chrs_in_gnomad_annotate_variants:
    input:
        vcf="output/gnomad-annotated-variants/annotated.vcf",
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
    shell: "vcfanno -p $(nproc --all) -lua config/custom.lua {input.vcfanno_config:q} {input.vcf:q} > {output:q}"

rule export_variant_stats:
    input: "output/clinvar-annotated-variants/annotated.vcf"
    output: "output/analysis-report/variant-stats.tsv"
    shell: "bcftools +split-vep -H -d \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%INFO/gno_id\t%INFO/dbsnp_id\t%INFO/clinvar_id\t%INFO/gno_af_all\t%CSQ\n' -A tab \
        {input:q} | sort | uniq > {output:q}"

rule export_pipeline_graph:
    output: "output/analysis-report/pipeline-summary.png"
    shell: "snakemake --rulegraph --forceall | dot -Tpng > {output:q}"

rule generate_variant_stats_plots:
    input: "output/analysis-report/variant-stats.tsv"
    output: "output/analysis-report/variants-coverage.png",
        "output/analysis-report/variants-quality.png",
        "output/analysis-report/variants-reported.png",
        "output/analysis-report/gnomad-af.png"
    shell: "Rscript src/summary-stats.r"

rule filter_common_variants_by_maf:
    input: "output/clinvar-annotated-variants/annotated.vcf.gz"
    output: "output/common-variants/maf_filtered.ge.0.05.bcf"
    shell: "bcftools view -i 'gno_af_all > 0.05' -Ob -o {output:q} {input:q}"

rule filter_annotated_variants_with_maf_lesser_than:
    input: "output/clinvar-annotated-variants/annotated.vcf.gz"
    output: "output/rare-variants/maf_filtered.lt.{threshold}.bcf"
    shell: "bcftools view -i 'gno_af_all < {wildcards.threshold}' -Ob -o {output:q} {input:q}"

rule convert_common_variants_to_plink:
    input: "output/common-variants/maf_filtered.ge.0.05.bcf"
    output: multiext("output/common-variants/plink", ".bed", ".bim", ".fam")
    shell: "plink --bcf {input:q} --allow-extra-chr --out output/common-variants/plink"

rule extract_principal_components_from_variants:
    input: "output/common-variants/plink.bed"
    output: multiext("output/common-variants/plink", ".eigenvec", ".eigenval")
    shell: "plink --bfile output/common-variants/plink --pca header tabs --out output/common-variants/plink"

rule generate_variants_pca_plot:
    input: "output/common-variants/plink.eigenvec"
    output: "output/common-variants/variants-pca.html"
    shell: "RSTUDIO_PANDOC='$(which pandoc)' Rscript src/pca-plots.r"

rule generate_quality_plots_from_fastq_files:
    input: "data/laboratorio/fastq/{path}.fastq.gz"
    output: "output/fastq-qc-plots/{path}_fastqc.html"
    shell: "fastqc {input:q} -o $(dirname {output:q})"

rule trim_fastq_read_files:
    input:
        fastq_r1="data/laboratorio/fastq/{read}_R1_001.fastq.gz",
        fastq_r2="data/laboratorio/fastq/{read}_R2_001.fastq.gz"
    output:
        trim_1p="output/trimmed-fastq-reads/{read}_1P.fastq.gz",
        trim_1u="output/trimmed-fastq-reads/{read}_1U.fastq.gz",
        trim_2p="output/trimmed-fastq-reads/{read}_2P.fastq.gz",
        trim_2u="output/trimmed-fastq-reads/{read}_2U.fastq.gz"
    shell: "trimmomatic PE -phred33 {input.fastq_r1:q} {input.fastq_r2:q} \
                -baseout 'output/trimmed-fastq-reads/{wildcards.read}.fastq.gz' \
                ILLUMINACLIP:config/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:50"

rule generate_quality_plots_from_trimmed_fastq_files:
    input: "output/trimmed-fastq-reads/{prefix}P.fastq.gz"
    output: "output/trimmed-fastq-qc-plots/{prefix}P_fastqc.html"
    shell: "fastqc {input:q} -o $(dirname {output:q})"

rule generate_phenotype_file:
    output: "output/phenotype-info/samples.pheno"
    shell: "Rscript src/make-pheno.r > {output:q}"