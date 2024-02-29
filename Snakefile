
from os.path import join, exists
import sys
import snakemake
import time


# Import config
# Load config file
#if not config:
 #   raise SystemExit(
  #      "Config not found, are you sure your profile contains a path to an "
   #     "existing configfile?"
    #)

#with open("config.yaml", "r") as config_file:
#    config = yaml.safe_load(config_file)


configfile: "config.yaml"
in_dir=config["input_dir"]
SAMPLES=glob_wildcards(in_dir+"/{sample}.fastq")
classifier=config["classifier"]
print(SAMPLES)
type="kraken"

rule all:
    input:
         expand(config["out_dir"]+"/classification_{classifier}/{sample}.kraken.report",sample=SAMPLES.sample,classifier=classifier),
         expand(config["out_dir"]+"/classification_{classifier}/{sample}.kraken",sample=SAMPLES.sample,classifier=classifier),
         expand(config["out_dir"]+"/classification_{classifier}/{sample}.kraken.bracken_species.report",sample=SAMPLES.sample,classifier=classifier) if config["bracken_options"]["run_bracken"] else "run.txt",
         #expand("/filtered_"+in_dir+"/{sample}.fastq", sample=SAMPLES.sample) if config["quality"]["perform"] else "run.txt",
         #expand("/aligned_"+in_dir+"/{sample}.fastq", sample=SAMPLES.sample) if config["run_align"] else "run.txt"

# Rule: Generate MultiQC report
rule multiqc:
    input:
        fastqcs=expand(in_dir+"/{sample}.fastq", sample=SAMPLES.sample)
    output:
        output_dir=directory("/multiqc_out"),
        multiqc_output="./multiqc_out/multiqc_report.html"
    shell:
        """
         fastqc -o {output.output_dir} {input.fastqcs} 
         multiqc -o {output.multiqc_output} {output.output_dir}"""


# Rule: Quality filtering with NanoFilt
rule quality_filter:
    input:
        in_dir+"/{sample}.fastq"
    output:
        directory("/filtered_"+in_dir+"/{sample}.fastq")
    params:
        quality_threshold=config["quality"].get("threshold",0),
        headcrop=config["quality"].get("head_crop",0),
        tailcrop=config["quality"].get("tail_crop",0),
        max_length=config["quality"].get("max_length",0)
    run:
       shell_cmd = "NanoFilt"	
       if params.quality_threshold:
            shell_cmd += " -q {params.quality_threshold}"
       if params.headcrop:
            shell_cmd += " --headcrop {params.headcrop}"
       if params.tailcrop:
            shell_cmd += " --tailcrop {params.tailcrop}"
       if params.max_length:
           shell_cmd += " --maxlength {params.max_length}"
       shell_cmd += " {input} > {output}; ln -s {in_dir} {output}"
       shell(shell_cmd)


# Rule: Align to human reference and exclude human reads + we filter resulting SAM file creating BAM file containing only the reads that did not align
# to the human genome reference
rule align_to_host:
    input:
        fastq=in_dir+"/{sample}.fastq",
        host_reference=config.get("host_reference")
    output:
        directory("/aligned_"+in_dir+"/{sample}.fastq")
    shell:
        """
        if [ -n "{input[host_reference]}" ]; then
            minimap2 -ax map-ont {input[host_reference]} {input[fastq]} | 
            samtools view -b -f 4 - > {output}/{wildcards.sample}_human_reads.bam && 
            samtools fastq -f 12 -1 {output}/{wildcards.sample}.fastq {output}/{wildcards.sample}_human_reads.bam
            ln -s {input[fastq]} {output}/{wildcards.sample}.fastq
            """

rule download_kraken_database:
    output:
        hash=join(config["kraken_options"].get("db","kraken_dbs"),"hash.k2d"),
        opts=join(config["kraken_options"].get("db","kraken_dbs"),"opts.k2d"),
        taxo=join(config["kraken_options"].get("db","kraken_dbs"),"taxo.k2d")
    params:
        db_root=config["kraken_build_options"],

        add_genome=config["kraken_build_options"].get("add_genome",0),
        method=config["kraken_build_options"]["method"],
        references=config["kraken_build_options"]["references"]
    container: "singularity_env.sif"
    shell:
        """
        bash ./scripts/download_kraken2_db.sh --method {params.db_root} --references {params.references} --output_dir {params.db_root}
        """






# Rule: Index reference for Kraken2
#rule index_kraken:
#    output:
#        hash="{kraken_db}/hash.k2d",
#        opts= "{kraken_db}/opts.k2d",
#        taxo="{kraken_db}/taxo.k2d"
#
#    params:
#        db_path=config["kraken_options"].get("db", "kraken_dbs"),
#        collection=config["kraken_build_options"]["collection"]
#    shell:
#        """
#            list=("greengenes" "silva" "rdp")
#            if [[ " {list[@]} " =~ " {params.collection} " ]]; then
#                kraken-build --special {params.collection} --db {params.db_path}
#            else
#                kraken-build --download-library {params.collection} --db {params.db_path}
#            fi
#
#        """



rule index_centrifuge_db:
    output:
        ex1 = join(config["centrifuge_build_options"]["centrifuge_db"],config["centrifuge_build_options"]["reference_name"],".1.cf"),
        ex2 = join(config["centrifuge_build_options"]["centrifuge_db"],config["centrifuge_build_options"]["reference_name"],".2.cf"),
        ex3 = join(config["centrifuge_build_options"]["centrifuge_db"],config["centrifuge_build_options"]["reference_name"],".3.cf")

    input:
        conv=config["centrifuge_build_options"]["centrifuge_db"]+"/seqid2taxid.map",
        tree=config["centrifuge_build_options"]["centrifuge_db"]+"/taxonomy/nodes.dmp",
        name_table=config["centrifuge_build_options"]["centrifuge_db"]+"/taxonomy/names.dmp",
        fa_path=config["centrifuge_build_options"]["centrifuge_db"]+"/library/"
    params:
        threads=config["centrifuge_build_options"].get("threads","1"),
        reference_name=config["centrifuge_build_options"].get("reference_name","ex")
    container: "singularity_env.sif"
    shell:
        """
        cat {input.fa_path}/*/*.fna |
        centrifuge-build -p {params.threads} --conversion-table {input[conv]} \
                         --taxonomy-tree {input[tree]} --name-table {input[name_table]} \
                          {params[reference_name]}
        """


# Rule: Kraken2 classification_{classifier}
rule kraken2:
    input:
        fastq=in_dir+"/{sample}.fastq",
        kraken_db=config["kraken_options"].get("db","kraken_dbs"),
        hash=config["kraken_options"].get("db","kraken_dbs")+"/hash.k2d",
        opts=config["kraken_options"].get("db","kraken_dbs")+"/opts.k2d",
        taxo=config["kraken_options"].get("db","kraken_dbs")+"/taxo.k2d"
    output:
        krak = join(config["out_dir"], "classification_kraken2/{sample}.kraken"),
        krak_report = join(config["out_dir"], "classification_kraken2/{sample}.kraken.report")
    container: "singularity_env.sif"
    shell:
        """
           time kraken2 --db {input[kraken_db]} --threads {config[threads]} --output {output.krak} --report {output.krak_report} {input.fastq} --use-names
        """
rule centrifuge:
    input:
        fastq=in_dir+"/{sample}.fastq",
        centrifuge_db=join(config["centrifuge_build_options"]["centrifuge_db"],config["centrifuge_build_options"]["reference_name"]),
        ex1 = join(config["centrifuge_build_options"]["centrifuge_db"],config["centrifuge_build_options"]["reference_name"],".1.cf"),
        ex2 = join(config["centrifuge_build_options"]["centrifuge_db"],config["centrifuge_build_options"]["reference_name"],".2.cf"),
        ex3 = join(config["centrifuge_build_options"]["centrifuge_db"],config["centrifuge_build_options"]["reference_name"],".3.cf"),
    output:
        report_file=join(config["out_dir"], "classification_centrifuge/{sample}.report"),
        stdout=join(config["out_dir"], "classification_centrifuge/{sample}.centrifuge"),
        kraken_report=join(config["out_dir"],"classification_centrifuge/{sample}.kraken.report")
    params:
        threads=config["centrifuge_options"].get("threads",1)
    container: "singularity_env.sif"
    shell:"""
    centrifuge -x {input.centrifuge_db} -S {output.stdout} --report-file {output.report_file} -f {input.fastq} 
    centrifuge-kreport -x {input.centrifuge_db} {output.report_file}  > {output.kraken_report}

    """

rule bracken:
    input:
     krak_report = join(config["out_dir"], "classification_{classifier}/{sample}.kraken.report"),
     krak = join(config["out_dir"], "classification_{classifier}/{sample}.kraken")
    output:
     join(config["out_dir"], "classification_{classifier}/{sample}.kraken.bracken_{level}.report")
    params:
     kraken_db = config["kraken_options"].get("db", "kraken_dbs"),
     readlen = config["bracken_options"]["read_length"],
     threshold = config["bracken_options"].get("threshold",10),
     level = config["bracken_options"].get("taxonomic_level","S"),
     out = join(config["out_dir"],"classification_{classifier}/{sample}.bracken.report")
    container: "singularity_env.sif"
    shell: """
       bracken -d {params.kraken_db} -i {input.krak_report} -o {params.out} -w {output} -l {params.level} -t {params.threshold} -l {params.readlen}
    """






#Downstream analysis with R
rule combine_kreport:
    input:
        expand(config["out_dir"]+"/classification_{classifier}/{sample}.kraken.report",sample=SAMPLES.sample,classifier=classifier)
    output:
        config["out_dir"]+f"/classification_{classifier}/combined.kraken.report"
    shell:
        """
        python combine_kreports.py -r {input} -o {output} --display-headers
        """

rule kreport2krona:
    input:
        config["out_dir"]+f"/classification_{classifier}"+"/{sample}.kraken.report"
    output:
        join(config["out_dir"],"krona_results/{sample}.krona")
    shell:
        """
        python kreport2krona.py -r {input} -o {output}
        """

rule filter_bracken:
    input:
        join(config["out_dir"], f"classification_{classifier}","/{sample}.bracken_{level}.report")
    output:
        join(config["out_dir"],"filtered_bracken/{sample}.bracken_{level}.report")
    params:
        to_=config["filter_bracken"]["filter"],
        list=config["filter_bracken"]["exclude"]
    shell:
        """
        python filter_bracken_out.py -i {input} -o {output} --{params.to_} {params.list}
        """

rule make_biom:
    output:
        join(config["out_dir"],f"analysis/{type}_table.biom")
    input:
        expand(config["out_dir"]+f"classification_{classifier}"+"/{sample}."+f"{type}.report",sample=SAMPLES.sample, classifier=classifier)

    shell:
        """
        kraken-biom {input} -o {output}
        """

rule make_phyloseq:
    input:
       biom_file=join(config["out_dir"],"analysis/{type}_table.biom")
    output:
        join(config["out_dir"],"analysis/phyloseq_{type}.Rdata")
    singularity: "docker://Zagh05/MetaLung:metalung"
    script:
        'scripts/make_phyloseq.R'

rule alpha_diversity:
    input:
        phylo_obj=join(config["out_dir"],f"classification_{classifier}"+"/phyloseq_{type}.Rdata")
    output:
        join(config["out_dir"], "analysis/alpha_diversity_{type}.png")
    params:
        measures = config["alpha_diversity"]["measures"],
        title = config["alpha_diversity"]["title"],
        color = config["alpha_diversity"]["color"]
    singularity: "docker://Zagh05/MetaLung:metalung"
    script:
        'scripts/alpha_diversity.R'

rule beta_diversity:
    input:
        phylo_obj=join(config["out_dir"],f"classification_{classifier}"+"/phyloseq_{type}.Rdata")
    output:
        join(config["out_dir"], "analysis/beta_diversity_{type}.png")
    params:
        method = config["beta_diversity"]["method"],
        title = config["beta_diversity"]["title"],
        color = config["beta_diversity"].get("color","NA"),
        shape = config["beta_diversity"].get("shape","NA"),
        type = config["beta_diversity"]["type"],
        wrap = config["beta_diversity"].get("wrap","NA")

    singularity: "docker://Zagh05/MetaLung:metalung"
    script:
        'scripts/beta_diversity.R'


rule taxonomic_barplots:
    input:
        phylo_obj=join(config["out_dir"],f"classification_{classifier}"+"/phyloseq_{type}.Rdata")
    output:
        join(config["out_dir"],"analysis/taxonomic_barplots_{type}.pdf")
    params:
        title = config["taxonomic_barplots"]["title"],
        groups = config["taxonomic_barplots"]["groups"],
        tax_ranks = config["taxonomic_barplots"]["tax_ranks"],
        abundance = config["taxonomic_barplots"]["abundance"],
        abundance_threshold = config["taxonomic_barplots"]["abundance_threshold"]
    singularity: "docker://Zagh05/MetaLung:metalung"
    script:
        'scripts/taxonomic_barplots.R'

rule abund_heatmap:
    input:
        phylo_obj=join(config["out_dir"], f"classification_{classifier}"+"/phyloseq_{type}.Rdata")
    output:
        join(config["out_dir"],"analysis/abund_heatmap_{type}.png")
    params:
        subset  =   config["abund_heatmap"]["subset"],
        method =   config["abund_heatmap"]["method"],
        distance =   config["abund_heatmap"]["distance"],
        sample_label =  config["abund_heatmap"]["sample_label"],
        taxa_label =   config["abund_heatmap"]["taxa_label"],
        wrap =   config["abund_heatmap"].get("wrap",0)
    singularity: "docker://Zagh05/MetaLung:metalung"
    script:
        'scripts/abund_heatmap.R'

rule differential_abundance:
    input:
        phylo_obj=join(config["out_dir"],f"classification_{classifier}"+"/phyloseq_{type}.Rdata")
    output:
        join(config["out_dir"],"analysis/differential_abundance/"+"compositional_PCA_plot_{type}.pdf")

    params:
        subset = config["diff_abund"].get("subset",0),
        group = config["diff_abund"]["group"],
        output_dir = join(config["out_dir"],"differential_abundance")
    singularity: "docker://Zagh05/MetaLung:metalung"
    script:
        'scripts/differential_abundance.R'