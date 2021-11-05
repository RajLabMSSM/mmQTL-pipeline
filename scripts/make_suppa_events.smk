shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; module purge; conda activate snakemake; ml R/4.0.3;")

GTF = config["GTF"]
prefix = config["prefix"]

rule all:
    input: 
        events = expand(prefix + ".events_{event_type}_strict.ioe", event_type = ["SE", "MX","RI","AF", "AL", "A3", "A5"]),
        total = prefix + ".all_suppa_events.ioe"

rule SUPPA_events:
    input:
        config["GTF"]
    output:
        events = expand(prefix + ".events_{event_type}_strict.ioe", event_type = ["SE", "MX","RI","AF", "AL", "A3", "A5"]),
        total = prefix + ".all_suppa_events.ioe"
    params:
        prefix = prefix + ".events"
    shell:
        "conda activate isoseq-pipeline;"
        "suppa.py generateEvents -i {input} -o {params.prefix} -e SE SS MX RI FL -f ioe --pool-genes;"
        " awk 'FNR==1 && NR!=1 {{ while (/^<header>/) getline; }} 1 {{print}}' {output.events} > {output.total}"

