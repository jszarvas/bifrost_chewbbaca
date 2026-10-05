#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback

from pathlib import Path
from bifrostlib import common
from bifrostlib.datahandling import SampleReference
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import ComponentReference
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from snakemake.io import directory
import datetime
from hashlib import md5
os.umask(0o2)

def md5_suffix(value: str, length: int = 8) -> str:
    return md5(value.encode("utf-8")).hexdigest()[:length]

print(config)

try:
    component_ref = ComponentReference(name=config['component_name'])
    component: Component = Component.load(reference=component_ref)
    if component is None:
        raise Exception("invalid component passed")

    # Get species
    detected_species = config["species"]

    # Mapping dictionary from component config
    schema_mapping = component["options"]["chewbbaca_species_mapping"]["schema"]

    # Resolve schema name
    SCHEMA_NAME = schema_mapping[detected_species]

    # Resolve schema directory
    SCHEMA_DIR = os.path.join(os.environ["BIFROST_CG_MLST_DIR"], "schemes", SCHEMA_NAME)

    # for calling alleles with chewBBACA 
    # Mapping dictionary from component config
    schema_call = component["options"]["chewbbaca_species_call"]["schema"]

    # Resolve schema name
    SCHEMA_NAME_CALL = schema_call[detected_species]

    # Resolve schema directory
    SCHEMA_DIR_CALL = os.path.join(os.environ["BIFROST_CG_MLST_DIR"], "schemes", SCHEMA_NAME_CALL)

    run_outputdir = config['rundir']

    # make samples+schema specific workdir
    if config['samples']:
        samples_str = "\n".join(config['samples'])
        config_str = f"{samples_str}\n{SCHEMA_NAME_CALL}"
        suffix = md5_suffix(config_str)
        chewbbaca_workdir = os.path.join(run_outputdir, f"chewbbaca_{suffix}")
    else:
        raise ValueError(f"No samples found in config")

    print(f"[INFO] {chewbbaca_workdir}")

    samples = {}
    samplecomponents = {}
    long_name_lookup = {}
    short_names = {}
    for sample_name in config['samples']:
        sn = sample_name.split("___")[1]
        long_name_lookup[sn] = sample_name
        short_names[sample_name] = sn
        sample_ref = SampleReference(_id=config.get('sample_id', None), name=sample_name)
        samples[sample_name]: Sample = Sample.load(sample_ref)
        if samples[sample_name] is None:
            raise Exception("invalid sample passed")

        samplecomponent_ref = SampleComponentReference(
            name=SampleComponentReference.name_generator(samples[sample_name].to_reference(), component.to_reference())
        )
        samplecomponents[sample_name] = SampleComponent.load(samplecomponent_ref)
        if samplecomponents[sample_name] is None:
            samplecomponents[sample_name] = SampleComponent(
                sample_reference=samples[sample_name].to_reference(),
                component_reference=component.to_reference()
            )

        common.set_status_and_save(samples[sample_name], samplecomponents[sample_name], "Running")
    
except Exception as error:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

print(f"[INFO] {run_outputdir}")

envvars:
    "BIFROST_INSTALL_DIR",
    "BIFROST_CG_MLST_DIR",
    "CONDA_PREFIX",
    "BIFROST_CPUS_BIG"

JOB_CPUS = 4

rule all:
    input:
        expand(f"{run_outputdir}/{{sample}}/{component['name']}/datadump_complete", sample = config['samples'])
    params:
        input_samples = config['samples']
    run:
        for sample_name in params.input_samples:
            common.set_status_and_save(samples[sample_name], samplecomponents[sample_name], "Success")

rule set_time_start:
    output:
        start_file = f"{run_outputdir}/{{sample}}/{component['name']}/time_start.txt"
    run:
        import time
        with open(output.start_file, "w") as fh:
            fh.write(str(time.time()))

rule setup:
    input:
        expand(f"{run_outputdir}/{{sample}}/{component['name']}/time_start.txt", sample = config['samples'])
    output:
        init_file = touch(f"{run_outputdir}/{{sample}}/{component['name']}/initialized")
    run:
        samplecomponents[wildcards.sample]['path'] = os.path.dirname(output.init_file)
        samplecomponents[wildcards.sample].save()

#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************

rule_name = "blast_locus_call"
rule blast_locus_call:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.out.log",
        err_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{run_outputdir}/{{sample}}/{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.set_time_start.output.start_file,
        genome = lambda wc: f"{samples[wc.sample]['categories']['contigs']['summary']['data']}"
    params:
        samplecomponent_ref_json = lambda wc: samplecomponents[wc.sample].to_reference().json,
        chewbbaca_schemes = f"{os.environ['BIFROST_CG_MLST_DIR']}/schemes/"
    threads: JOB_CPUS
    output:
        locus_call_results = directory(f"{run_outputdir}/{{sample}}/{component['name']}/blast_locus_call_results"),
        locus_calls = f"{run_outputdir}/{{sample}}/{component['name']}/blast_locus_call_results/locus_calls.fa",
        locus_call_done = f"{run_outputdir}/{{sample}}/{component['name']}/blast_locus_call_done"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "rule__blast_locuscall.py")

rule set_blast_time_end:
    input:
        rules.blast_locus_call.output.locus_call_done
    output:
        blast_end_file = f"{run_outputdir}/{{sample}}/{component['name']}/blast_time_end.txt"
    run:
        import time
        with open(output.blast_end_file, "w") as fh:
            fh.write(str(time.time()))

rule set_chewbbaca_time_start:
    input:
        rules.blast_locus_call.output.locus_call_done,
    output:
        chewbbaca_start_file = f"{run_outputdir}/{{sample}}/{component['name']}/chewbbaca_time_start.txt"
    run:
        import time
        with open(output.chewbbaca_start_file, "w") as fh:
            fh.write(str(time.time()))

rule_name = "gather_loci_for_chewbbaca"
rule gather_loci_for_chewbbaca:
    input:
        genome = lambda wc: f"{run_outputdir}/{long_name_lookup[wc.short]}/{component['name']}/blast_locus_call_results/locus_calls.fa"
    output:
        chewbbaca_input = f"{chewbbaca_workdir}/input/{{short}}.fa"
    threads: JOB_CPUS
    shell:
        r"""
        ln -sf $(realpath {input.genome}) {output.chewbbaca_input}
        """

CHEWB_SCRIPT = os.path.realpath(
    os.path.join(os.path.dirname(workflow.snakefile), "../chewBBACA/CHEWBBACA/chewBBACA.py")
)

rule_name = "run_chewbbaca_on_batch"
rule run_chewbbaca_on_batch:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{chewbbaca_workdir}/log/{rule_name}.out.log",
        err_file = f"{chewbbaca_workdir}/log/{rule_name}.err.log",
    benchmark:
        f"{chewbbaca_workdir}/benchmarks/{rule_name}.benchmark"
    input:
        expand(f"{run_outputdir}/{{sample}}/{component['name']}/chewbbaca_time_start.txt", sample = config['samples']),
        genomes = expand(f"{chewbbaca_workdir}/input/{{short}}.fa", short = list(long_name_lookup.keys()))
    output:
        chewbbaca_results = directory(f"{chewbbaca_workdir}/chewbbaca_results"),
        chewbbaca_done = f"{chewbbaca_workdir}/chewbbaca_done"
    params:
        chewbbaca_workdir = chewbbaca_workdir,
        chewbbaca_script = CHEWB_SCRIPT,
        schema_name = SCHEMA_NAME_CALL,
        schema_dir = SCHEMA_DIR_CALL
    threads: 16
    shell:
        r"""
        set -euo pipefail

        echo "DEBUG: schema_name = {params.schema_name}"
        echo "DEBUG: schema_dir  = {params.schema_dir}"

        mkdir -p {output.chewbbaca_results}
        echo "{params.schema_name}" > {output.chewbbaca_results}/schema

        {params.chewbbaca_script} AlleleCall \
            -i {params.chewbbaca_workdir}/input \
            -g "{params.schema_dir}" \
            -o {output.chewbbaca_results}/output \
            --cpu {threads} \
            --cds --wait-time 120 --lock-stale 2 \
            1>> {log.out_file} \
            2>> {log.err_file}

        touch {output.chewbbaca_done}
        """

rule_name = "separate_chewbbaca_batch"
rule separate_chewbbaca_batch:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.out.log",
        err_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.err.log",
    input:
        rules.run_chewbbaca_on_batch.output.chewbbaca_done
    output:
        chewbbaca_results = directory(f"{run_outputdir}/{{sample}}/{component['name']}/chewbbaca_results"),
        chewbbaca_alleles = f"{run_outputdir}/{{sample}}/{component['name']}/chewbbaca_results/output/results_alleles.tsv",
        chewbbaca_stats = f"{run_outputdir}/{{sample}}/{component['name']}/chewbbaca_results/output/results_statistics.tsv",
        chewbbaca_schema = f"{run_outputdir}/{{sample}}/{component['name']}/chewbbaca_results/schema",
        chewbbaca_done = f"{run_outputdir}/{{sample}}/{component['name']}/chewbbaca_results/chewbbaca_done"
    params:
        chewbbaca_workdir = chewbbaca_workdir,
        schema_name = SCHEMA_NAME_CALL,
        short_name = lambda wc: short_names[wc.sample]
    threads: 1
    shell:
        r"""
        set -euo pipefail

        echo "DEBUG: chewbbaca_workdir = {params.chewbbaca_workdir}"

        mkdir -p {output.chewbbaca_results}
        ln -s {params.chewbbaca_workdir} {output.chewbbaca_results}/chewbbaca_batch
        echo "{params.schema_name}" > {output.chewbbaca_schema}

        grep -P "FILE|{params.short_name}" {params.chewbbaca_workdir}/chewbbaca_results/output/results_alleles.tsv > {output.chewbbaca_alleles}
        grep -P "FILE|{params.short_name}" {params.chewbbaca_workdir}/chewbbaca_results/output/results_statistics.tsv > {output.chewbbaca_stats}

        touch {output.chewbbaca_done}
        """


#* Dynamic section: end ****************************************************************************

rule set_time_end:
    input:
        rules.run_chewbbaca_on_batch.output.chewbbaca_done
    output:
        end_file = f"{run_outputdir}/{{sample}}/{component['name']}/time_end.txt"
    run:
        import time
        with open(output.end_file, "w") as fh:
            fh.write(str(time.time()))

rule_name = "git_version"
rule git_version:
    log:
        out_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.out.log",
        err_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.err.log",
    input:
        rules.setup.output.init_file
    output:
        git_hash = f"{run_outputdir}/{{sample}}/{component['name']}/git_hash.txt"
    run:
        import subprocess, os

        snake_dir = os.path.dirname(workflow.snakefile)

        # Best effort: get commit hash; if not a git repo, write "-"
        try:
            git_hash = subprocess.check_output(
                ["git", "-C", snake_dir, "rev-parse", "HEAD"],
                stderr=subprocess.STDOUT,
                text=True
            ).strip()
        except Exception as e:
            git_hash = "-"
            os.makedirs(os.path.dirname(log.err_file), exist_ok=True)
            with open(log.err_file, "a") as fh:
                fh.write(f"[git_version] Could not determine git hash from {snake_dir}: {e}\n")

        with open(output.git_hash, "w") as fh:
            fh.write(str(git_hash))

rule dump_info:
    input:
        start_file = rules.set_time_start.output.start_file,
        end_file = rules.set_time_end.output.end_file,
        git_hash = rules.git_version.output.git_hash,
        blast_end = rules.set_blast_time_end.output.blast_end_file,
        chewbbaca_start = rules.set_chewbbaca_time_start.output.chewbbaca_start_file,
    output:
        runtime_flag = touch(f"{run_outputdir}/{{sample}}/{component['name']}/runtime_set")
    run:
        import time
        from bifrostlib.datahandling import SampleComponent

        with open(input.start_file) as fh:
            t_start = float(fh.read().strip())
        with open(input.end_file) as fh:
            t_end = float(fh.read().strip())
        with open(input.git_hash) as fh:
            git_hash = str(fh.read().strip())

        with open(input.blast_end) as fh: t_b_end = float(fh.read().strip())

        with open(input.chewbbaca_start) as fh: t_c_start = float(fh.read().strip())

        runtime_minutes = (t_end - t_start) / 60.0
        print(f"runtime in minutes {runtime_minutes}")

        sc = SampleComponent.load(samplecomponents[wildcards.sample].to_reference())
        sc["time_start"] = datetime.datetime.fromtimestamp(t_start).strftime("%Y-%m-%d %H:%M:%S")
        sc["time_end"] = datetime.datetime.fromtimestamp(t_end).strftime("%Y-%m-%d %H:%M:%S")
        sc["time_running"] = round(runtime_minutes, 3)
        sc["git_hash"] = git_hash

        sc["blast_time_running"] = round(((t_b_end - t_start) / 60.0), 3)
        sc["chewbbaca_time_running"] = round(((t_end - t_c_start) / 60.0), 3)

        sc.save()


#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.out.log",
        err_file = f"{run_outputdir}/{{sample}}/{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{run_outputdir}/{{sample}}/{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.separate_chewbbaca_batch.output.chewbbaca_results,
        rules.dump_info.output.runtime_flag
    output:
        f"{run_outputdir}/{{sample}}/{component['name']}/datadump_complete"
    params:
        samplecomponent_ref_json = lambda wc: samplecomponents[wc.sample].to_reference().json
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------







