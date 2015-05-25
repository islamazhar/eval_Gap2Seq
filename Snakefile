"""
Submit this job on uppmax as:
    snakemake --debug --keep-going -j 999 --wait-for-files 30 --cluster "sbatch -A {params.account} -p {params.partition} -n {params.n}  -t {params.runtime} -C {params.memsize} -J {params.jobname} --mail-type={params.mail_type} --mail-user={params.mail}"
"""
shell.prefix("set -o pipefail; ")
configfile: "config.json"



####################################################
########## standard python functions ###############
####################################################

import re
import os

def parse_quast(quast_report):
    lines = open(quast_report, 'r').readlines()
    for line in lines:
        if re.match(r"Total length \(\>\= 0 bp\)",line):
            assembly_size = int(line.strip().split()[6])
        if re.match(r"NG50",line):
            NG50 = int(line.strip().split()[2])
        if re.match(r"NGA50",line):
            NGA50 = int(line.strip().split()[2]) 
        if re.match(r"# misassemblies",line):
            large_misassm = int(line.strip().split()[3]) 
        if re.match(r"# local misassemblies",line):
            local_misassm = int(line.strip().split()[4]) 
        if re.match(r"ESIZE_ASSEMBLY",line):
            ESIZE_ASSEMBLY = float(line.strip().split()[2])
        if re.match(r"ESIZE_GENOME",line):
            ESIZE_GENOME = float(line.strip().split()[2]) 
        if re.match(r"CORR_ESIZE_ASSEMBLY",line):
            CORR_ESIZE_ASSEMBLY = float(line.strip().split()[2])
        if re.match(r"CORR_ESIZE_GENOME",line):
            CORR_ESIZE_GENOME = float(line.strip().split()[2]) 

    try:
        NGA50
    except UnboundLocalError:
        NGA50 = "."
    try:
        NG50
    except UnboundLocalError:
        NG50 = "."

    misassmblies = large_misassm #+ local_misassm
    return(misassmblies, NG50, NGA50, assembly_size, ESIZE_ASSEMBLY, ESIZE_GENOME, CORR_ESIZE_ASSEMBLY, CORR_ESIZE_GENOME)

def  parse_gnu_time(stderr_file):
    lines = open(stderr_file, 'r').readlines()

    for l in lines:
        usertime_match =  re.search('User time \(seconds\): [\d.]+', l)
        wct_match = re.search('Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): [\d.:]+', l) 
        mem_match = re.search('Maximum resident set size \(kbytes\): [\d.:]+', l) 
        if usertime_match:
            usertime = float(usertime_match.group().split(':')[1].strip())
        if wct_match:
            wallclocktime = wct_match.group().split()[7]
        if mem_match:
            mem_tmp = int(mem_match.group().split()[5])
            memory_gb = mem_tmp / 4000000.0 

    vals = list(map(lambda x: float(x), wallclocktime.split(":") ))
    if len(vals) == 3:
        h,m,s = vals
        tot_wallclock_secs = h*3600.0 + m*60.0 + s
    elif len(vals) == 2:
        m,s = vals
        tot_wallclock_secs = m*60.0 + s

    return usertime, tot_wallclock_secs, memory_gb


def performance_input(wildcards):
  input_list_to_performance_latex_table = []
  
  if wildcards.dataset == "staph":
    for gapfiller in config["GAPFILLERS"]:
      for assembly in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
        fname = config["OUTBASE"]+"staph/{0}_{1}_time_and_mem.txt".format(gapfiller, assembly)
        if os.path.isfile(fname):
          input_list_to_performance_latex_table.append(fname)

  if wildcards.dataset == "rhodo":
    for gapfiller in config["GAPFILLERS"]:
      for assembly in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
        fname = config["OUTBASE"]+"rhodo/{0}_{1}_time_and_mem.txt".format(gapfiller, assembly)
        if os.path.isfile(fname):
          input_list_to_performance_latex_table.append(fname)


  if wildcards.dataset == "hs14":
    for gapfiller in config["GAPFILLERS"]:
      for assembly in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
        fname = config["OUTBASE"]+"hs14/{0}_{1}_time_and_mem.txt".format(gapfiller, assembly)
        if os.path.isfile(fname):
          input_list_to_performance_latex_table.append(fname)

  return input_list_to_performance_latex_table


def quality_input(wildcards):
  """
    Allows for missing files to input as some jobs may fail due to time constraints etc.
  """
  input_list_to_quality_latex_table = []

  if wildcards.dataset == "staph":
    for gapfiller in config["GAPFILLERS"]:
      for assembly in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
        fname = config["OUTBASE"]+"staph/quast_{0}_{1}.csv".format(gapfiller, assembly)
        if os.path.isfile(fname):
          input_list_to_quality_latex_table.append(fname)

  if wildcards.dataset == "rhodo":
    for gapfiller in config["GAPFILLERS"]:
      for assembly in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
        fname = config["OUTBASE"]+"rhodo/quast_{0}_{1}.csv".format(gapfiller, assembly)
        if os.path.isfile(fname):
          input_list_to_quality_latex_table.append(fname)

  if wildcards.dataset == "hs14":
    for gapfiller in config["GAPFILLERS"]:
      for assembly in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
        fname = config["OUTBASE"]+"hs14/quast_{0}_{1}.csv".format(gapfiller, assembly)
        if os.path.isfile(fname):
          input_list_to_quality_latex_table.append(fname)

  return input_list_to_quality_latex_table


def scaffolded_fasta_files(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    if dataset == "staph":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Allpaths-LG", "Bambus2", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          input_.append(config["OUTBASE"]+"{0}/{1}_{2}.fa".format(dataset, gapfiller, assembler))

    if dataset == "rhodo":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Allpaths-LG", "Bambus2", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          input_.append(config["OUTBASE"]+"{0}/{1}_{2}.fa".format(dataset, gapfiller, assembler))

    if dataset == "hs14":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Allpaths-LG", "Bambus2", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          input_.append(config["OUTBASE"]+"{0}/{1}_{2}.fa".format(dataset, gapfiller, assembler))

  return input_


def eval_csv_files(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    if dataset == "staph":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Allpaths-LG", "Bambus2", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          infile_name = config["OUTBASE"]+"{0}/{1}_{2}.fa".format(dataset, gapfiller, assembler)
          if os.path.isfile(infile_name):
            outfile_name = config["OUTBASE"]+"{0}/quast_{1}_{2}.csv".format(dataset, gapfiller, assembler)
            input_.append(outfile_name)

    if dataset == "rhodo":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          infile_name = config["OUTBASE"]+"{0}/{1}_{2}.fa".format(dataset, gapfiller, assembler)
          if os.path.isfile(infile_name):
            outfile_name = config["OUTBASE"]+"{0}/quast_{1}_{2}.csv".format(dataset, gapfiller, assembler)
            input_.append(outfile_name)


    if dataset == "hs14":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Bambus2", "Allpaths-LG", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          infile_name = config["OUTBASE"]+"{0}/{1}_{2}.fa".format(dataset, gapfiller, assembler)
          if os.path.isfile(infile_name):
            outfile_name = config["OUTBASE"]+"{0}/quast_{1}_{2}.csv".format(dataset, gapfiller, assembler)
            input_.append(outfile_name)

  return input_

def performance_csv_files(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    if dataset == "staph":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Allpaths-LG", "Bambus2", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          infile_name = config["OUTBASE"]+"{0}/{1}_{2}.stderr".format(dataset, gapfiller, assembler)
          if os.path.isfile(infile_name):
            outfile_name = config["OUTBASE"]+"{0}/{1}_{2}_time_and_mem.txt".format(dataset, gapfiller, assembler)
            input_.append(outfile_name)

    if dataset == "rhodo":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Allpaths-LG", "Bambus2", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          infile_name = config["OUTBASE"]+"{0}/{1}_{2}.stderr".format(dataset, gapfiller, assembler)
          if os.path.isfile(infile_name):
            outfile_name = config["OUTBASE"]+"{0}/{1}_{2}_time_and_mem.txt".format(dataset, gapfiller, assembler)
            input_.append(outfile_name)

    if dataset == "hs14":
      for gapfiller in config["GAPFILLERS"]:
        for assembler in ["ABySS", "ABySS2", "Allpaths-LG", "Bambus2", "CABOG", "MSR-CA", "SGA", "SOAPdenovo", "Velvet"]:
          infile_name = config["OUTBASE"]+"{0}/{1}_{2}.stderr".format(dataset, gapfiller, assembler)
          if os.path.isfile(infile_name):
            outfile_name = config["OUTBASE"]+"{0}/{1}_{2}_time_and_mem.txt".format(dataset, gapfiller, assembler)
            input_.append(outfile_name)

  return input_


def latex_tables(wildcards):
  input_= []

  for dataset in config["DATASETS"]:
    input_.append(config["OUTBASE"]+"performance_table_{0}.tex".format(dataset))
    input_.append(config["OUTBASE"]+"quality_table_{0}.tex".format(dataset))

      # if experiment == "rhodo":
      #   input_.append(config["OUTBASE"]+"performance_table_{0}_{1}.tex".format(experiment, contamine))
      #   input_.append(config["OUTBASE"]+"quality_table_{0}_{1}.tex".format(experiment, contamine))

      # if experiment == "hs14":
      #   input_.append(config["OUTBASE"]+"performance_table_{0}_{1}.tex".format(experiment, contamine))
      #   input_.append(config["OUTBASE"]+"quality_table_{0}_{1}.tex".format(experiment, contamine))

      # if experiment == "sim":
      #   input_.append(config["OUTBASE"]+"performance_table_{0}_{1}.tex".format(experiment, contamine))
      #   input_.append(config["OUTBASE"]+"quality_table_{0}_{1}.tex".format(experiment, contamine))

      # if experiment == "assemblathon3k":
      #   input_.append(config["OUTBASE"]+"performance_table_{0}_{1}.tex".format(experiment, contamine))
      #   input_.append(config["OUTBASE"]+"quality_table_{0}_{1}.tex".format(experiment, contamine))
  return input_

###########################################################
###########################################################


############## TARGET RULES ####################
rule all:
    input: latex_tables
    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

# rule test:
#     input: config["OUTBASE"]+"quality_table_sim_30.tex",
#            config["OUTBASE"]+"percentage_quality_table_sim_30.tex",
#            expand(config["OUTBASE"]+"sim/30/quast_{scaffolder}_NA.csv",scaffolder=config["SCAFFOLDERS"])
#     params: 
#         runtime="15:00",
#         memsize = "mem128GB",
#         partition = "core",
#         n = "1",
#         jobname="all",
#         account=config["SBATCH"]["ACCOUNT"],
#         mail=config["SBATCH"]["MAIL"],
#         mail_type=config["SBATCH"]["MAIL_TYPE"]

rule eval_gapfillers:
    input:  eval_csv_files,
            performance_csv_files
            #expand(config["OUTBASE"]+"{experiment}/performance_table.tex",  experiment=config["EXPERIMENTS"]),
            #expand(config["OUTBASE"]+"{experiment}/quality_table.tex",  experiment=config["EXPERIMENTS"]), 

    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

rule run_gapfillers:
    input: scaffolded_fasta_files  
    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]



############## TARGET RULES END ####################


rule Gap2Seq:
    input:  scaffolds = lambda wildcards: config[wildcards.dataset][wildcards.assembler]["SCAFFOLDS"],
            reads = lambda wildcards: config[wildcards.dataset]["MERGED_READS"]
    output: fasta= config["OUTBASE"]+"{dataset}/GAP2SEQ_{assembler}.fa",
            stderr= config["OUTBASE"]+"{dataset}/GAP2SEQ_{assembler}.stderr"
    log: config["OUTBASE"]+"{dataset}/GAP2SEQ_{assembler}.log"
    params: 
        runtime = lambda wildcards: config["SBATCH"][wildcards.dataset]["GAP2SEQ_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname = "{dataset}_{assembler}"+"_GAP2SEQ",
        account = config["SBATCH"]["ACCOUNT"],
        mail = config["SBATCH"]["MAIL"],
        mail_type = config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        stdout = config["OUTBASE"]+"{0}/GAP2SEQ_{1}.stdout".format(wildcards.dataset, wildcards.assembler)
        shell("{time} Gap2Seq -filled {output.fasta} -scaffolds {input.scaffolds} -reads {input.reads} 1>  {stdout} 2> {output.stderr} ")
        shell("cat {stdout} {output.stderr} > {log}")

rule GapCloser:
    input:  scaffolds = lambda wildcards: config[wildcards.dataset][wildcards.assembler]["SCAFFOLDS"],
            config = lambda wildcards: config[wildcards.dataset]["GAPCLOSER_CONFIG"]
    output: fasta= config["OUTBASE"]+"{dataset}/GAPCLOSER_{assembler}.fa",
            stderr= config["OUTBASE"]+"{dataset}/GAPCLOSER_{assembler}.stderr"
    log: config["OUTBASE"]+"{dataset}/GAPCLOSER_{assembler}.log"
    params: 
        runtime = lambda wildcards: config["SBATCH"][wildcards.dataset]["GAPCLOSER_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["n"],
        jobname = "{dataset}_{assembler}"+"_GAPCLOSER",
        account = config["SBATCH"]["ACCOUNT"],
        mail = config["SBATCH"]["MAIL"],
        mail_type = config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        stdout = config["OUTBASE"]+"{0}/GAPCLOSER_{1}.stdout".format(wildcards.dataset, wildcards.assembler)
        shell("{time} GapCloser -a {input.scaffolds} -b {input.config} -o {output.fasta} 1> {stdout} 2> {output.stderr} ")
        shell("cat {stdout} {output.stderr} > {log}")

rule GapFiller_bwa:
    input:  scaffolds = lambda wildcards: config[wildcards.dataset][wildcards.assembler]["SCAFFOLDS"],
            config = lambda wildcards: config[wildcards.dataset]["GAPFILLER_BWA_CONFIG"],
    output: fasta= config["OUTBASE"]+"{dataset}/GAPFILLER_BWA_{assembler}.fa",
            stderr= config["OUTBASE"]+"{dataset}/GAPFILLER_BWA_{assembler}.stderr"
    log: config["OUTBASE"]+"{dataset}/GAPFILLER_BWA_{assembler}.log"
    params: 
        runtime = lambda wildcards: config["SBATCH"][wildcards.dataset]["GAPFILLER_BWA_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_n"],
        jobname = "{dataset}_{assembler}"+"_GAPFILLER_BWA",
        account = config["SBATCH"]["ACCOUNT"],
        mail = config["SBATCH"]["MAIL"],
        mail_type = config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        stdout = config["OUTBASE"]+"{0}/GAPFILLER_BWA_{1}.stdout".format(wildcards.dataset, wildcards.assembler)
        shell("{time} GapFiller  -s {input.scaffolds} -l {input.config} -b bwa_{wildcards.dataset}_{wildcards.assembler} > {stdout} 2> {output.stderr} ")
        shell("cat {stdout} {output.stderr} > {log}")
        shell("mv bwa_{wildcards.dataset}_{wildcards.assembler}/bwa_{wildcards.dataset}_{wildcards.assembler}.gapfilled.final.fa {output.fasta}")
        shell("rm -r bwa_{wildcards.dataset}_{wildcards.assembler}")

rule GapFiller_bowtie:
    input:  scaffolds = lambda wildcards: config[wildcards.dataset][wildcards.assembler]["SCAFFOLDS"],
            config = lambda wildcards: config[wildcards.dataset]["GAPFILLER_BOWTIE_CONFIG"],
    output: fasta= config["OUTBASE"]+"{dataset}/GAPFILLER_BOWTIE_{assembler}.fa",
            stderr= config["OUTBASE"]+"{dataset}/GAPFILLER_BOWTIE_{assembler}.stderr"
    log: config["OUTBASE"]+"{dataset}/GAPFILLER_BOWTIE_{assembler}.log"
    params: 
        runtime = lambda wildcards: config["SBATCH"][wildcards.dataset]["GAPFILLER_BOWTIE_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_n"],
        jobname = "{dataset}_{assembler}"+"_GAPFILLER_BOWTIE",
        account = config["SBATCH"]["ACCOUNT"],
        mail = config["SBATCH"]["MAIL"],
        mail_type = config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        stdout = config["OUTBASE"]+"{0}/GAPFILLER_BOWTIE_{1}.stdout".format(wildcards.dataset, wildcards.assembler)
        shell("{time} GapFiller  -s {input.scaffolds} -l {input.config} -b bowtie_{wildcards.dataset}_{wildcards.assembler} > {stdout} 2> {output.stderr} ")
        shell("cat {stdout} {output.stderr} > {log}")
        shell("mv bowtie_{wildcards.dataset}_{wildcards.assembler}/bowtie_{wildcards.dataset}_{wildcards.assembler}.gapfilled.final.fa {output.fasta}")
        shell("rm -r bowtie_{wildcards.dataset}_{wildcards.assembler}")

rule ORIGINAL:
    input:  scaffolds = lambda wildcards: config[wildcards.dataset][wildcards.assembler]["SCAFFOLDS"],
    output: fasta= config["OUTBASE"]+"{dataset}/ORIGINAL_{assembler}.fa"
    params: 
        runtime="10:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname = "{dataset}_{assembler}_{contamine}"+"_ORIGINAL",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        shell("{time} cp {input.scaffolds} {output.fasta} ")

rule QUAST:
   input: scaffolds=config["OUTBASE"]+"{dataset}/{scaffolder}_{assembler}.fa",
            ref = lambda wildcards: config[wildcards.dataset]["REF"]
   output: nice_format=config["OUTBASE"]+"{dataset}/quast_{scaffolder}_{assembler}.csv" 
   params: 
       runtime=lambda wildcards: config["SBATCH"][wildcards.dataset]["quast_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_n"],
       jobname="quast_{dataset}_{scaffolder}",
       account=config["SBATCH"]["ACCOUNT"],
       mail=config["SBATCH"]["MAIL"],
       mail_type=config["SBATCH"]["MAIL_TYPE"]
   run:
       env = config["LOAD_PYTHON_ENV"]
       shell("{env}")
       python = config["PYTHON2"]
       #out=config['OUTBASE']
       path=config["quast_rules"]["path"]
       min_contig =  config["quast_rules"]["min_contig"]
       outpath="/tmp/{0}_{1}_{2}/QUAST/".format(wildcards.dataset, wildcards.scaffolder, wildcards.assembler)
       shell(" {python} {path}quast.py -R  {input.ref} -o {outpath} -s --no-plots {input.scaffolds} ") 
       misassmblies, N50, NA50, tot_length, ESIZE_ASSEMBLY, ESIZE_GENOME, CORR_ESIZE_ASSEMBLY, CORR_ESIZE_GENOME = parse_quast(outpath+"report.txt")
       #e_size = get_esize(input.scaffolds)
       print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(wildcards.assembler, wildcards.scaffolder, tot_length, N50, misassmblies,  NA50, ESIZE_ASSEMBLY, ESIZE_GENOME, CORR_ESIZE_ASSEMBLY, CORR_ESIZE_GENOME), file=open(output.nice_format, 'w'))    



rule time_and_mem:
   input:  stderr=config["OUTBASE"]+"{dataset}/{scaffolder}_{assembler}.stderr"
   output: outfile=config["OUTBASE"]+"{dataset}/{scaffolder}_{assembler}_time_and_mem.txt"
   params: 
       runtime="15:00",
       memsize = "'mem128GB|mem256GB|mem512GB'",
       partition = "core",
       n = "1",
       jobname="{assembler}"+"_time_and_mem",
       account=config["SBATCH"]["ACCOUNT"],
       mail=config["SBATCH"]["MAIL"],
       mail_type=config["SBATCH"]["MAIL_TYPE"]
   run:
       usertime, wallclocktime, memory_gb =  parse_gnu_time(input.stderr)
       print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(wildcards.dataset, wildcards.assembler, wildcards.scaffolder, usertime, wallclocktime, memory_gb), file=open(output.outfile, 'w') )
       

rule performace_latex_table:
   input: files=performance_input
   output: table=config["OUTBASE"]+"performance_table_{dataset}.tex"
   params: 
       runtime="15:00",
       memsize = "'mem128GB|mem256GB|mem512GB'",
       partition = "core",
       n = "1",
       jobname="performace_latex_table_{dataset}",
       account=config["SBATCH"]["ACCOUNT"],
       mail=config["SBATCH"]["MAIL"],
       mail_type=config["SBATCH"]["MAIL_TYPE"]
   run:
        table_file = open(output.table, 'w')
        print("{0} & {1} & {2} & {3} & {4} & {5}  \\\ \hline".format('organism', 'assembly', 'tool', 'user time', 'wall clock time', 'peak memory (Gb)'), file=table_file)
        # for file_ in input.files:
        #     line=open(file_,'r').readlines()[0]
        #     print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*line.strip().split()), file=table_file)


        prev_scaffolder = -1

        for file_ in input.files:
            line=open(file_,'r').readlines()[0]  
            vals = line.strip().split()
            #scaffolder = vals[2]
            experiment, assembler, scaffolder, user_time, wct, peak_mem = vals
            if prev_scaffolder == -1:
                sum_values = [experiment, "TOTAL", scaffolder, 0, 0, 0]
                completed_experiment_count = 0

            if scaffolder != prev_scaffolder and prev_scaffolder != -1:
                average_values = list(sum_values)
                average_values[0] = "AVERAGE"
                for i, item in enumerate(sum_values): 
                    try:
                        value = float(item)
                        average_values[i] = round((value/completed_experiment_count),1)
                        sum_values[i] = round(item,1)
                    except ValueError:
                        pass

                print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*map(lambda x: x, sum_values)), file=table_file)
                print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*map(lambda x: x, average_values)), file=table_file)
                experiment, assembler, scaffolder, user_time, wct, peak_mem = vals


                sum_values = [experiment, "TOTAL", scaffolder, 0, 0, 0]
                completed_experiment_count = 1
                for i, item in enumerate(vals): 
                    try:
                        value = float(item)
                        sum_values[i] = value
                    except ValueError:
                        pass

            else:
                completed_experiment_count += 1
                for i, item in enumerate(vals): 
                    try:
                        value = float(item)
                        sum_values[i] += value
                    except ValueError:
                        pass

            prev_scaffolder = scaffolder

            print("{0} & {1} & {2} & {3} & {4} & {5} \\\ ".format(*line.strip().split()), file=table_file)

        average_values = list(sum_values)
        average_values[0] = "AVERAGE"
        for i, item in enumerate(sum_values): 
            try:
                value = float(item)
                average_values[i] = round((value/completed_experiment_count),1)
                sum_values[i] = round(item,1)
            except ValueError:
                pass

        print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*map(lambda x: x, sum_values)), file=table_file)
        print("{0} & {1} & {2} & {3} & {4} & {5} \\\ \hline".format(*map(lambda x: x, average_values)), file=table_file)


rule quality_latex_table:
   input: files=quality_input 
   output: table=config["OUTBASE"]+"quality_table_{dataset}.tex"
   params: 
       runtime="15:00",
       memsize = "'mem128GB|mem256GB|mem512GB'",
       partition = "core",
       n = "1",
       jobname="quality_latex_table_{dataset}",
       account=config["SBATCH"]["ACCOUNT"],
       mail=config["SBATCH"]["MAIL"],
       mail_type=config["SBATCH"]["MAIL_TYPE"]
   run:
        table_file = open(output.table, 'w')
        print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9}  \\\ \hline".format('assembly', 'tool', 'assm-size', 'NG50', 'misassm', 'NGA50', '$E_{assm}$', '$E_{genome}$', '$EA_{assm}$' , '$EA_{genome}$'), file=table_file)

        prev_scaffolder = -1

        for file_ in input.files:
            line=open(file_,'r').readlines()[0]  
            vals = line.strip().split()
            #scaffolder = vals[2]
            assembler, scaffolder, totsize, NG50, misassm, NGA50, E_a, E_g, EA_a, EA_g = vals
            if prev_scaffolder == -1:
                sum_values = ["TOTAL", scaffolder, 0, 0, 0, 0, 0, 0, 0, 0]
                completed_experiment_count = 0

            if scaffolder != prev_scaffolder and prev_scaffolder != -1:
                average_values = list(sum_values)
                average_values[0] = "AVERAGE"
                for i, item in enumerate(sum_values): 
                    try:
                        value = float(item)
                        average_values[i] = round((value/completed_experiment_count),1)
                        sum_values[i] = round(item,0) 
                    except ValueError:
                        pass

                print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9}  \\\ \hline".format(*map(lambda x: x,sum_values)), file=table_file)
                print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9}  \\\ \hline".format(*map(lambda x: x, average_values)), file=table_file)

                assembler, scaffolder, totsize, NG50, misassm, NGA50, E_a, E_g, EA_a, EA_g = vals
                sum_values = ["TOTAL", scaffolder, 0, 0, 0, 0, 0, 0, 0, 0]
                completed_experiment_count = 1
                for i, item in enumerate(vals): 
                    try:
                        value = float(item)
                        sum_values[i] = value
                    except ValueError:
                        pass

            else:
                completed_experiment_count += 1
                for i, item in enumerate(vals): 
                    try:
                        value = float(item)
                        sum_values[i] += value
                    except ValueError:
                        pass

            prev_scaffolder = scaffolder

            print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9} \\\ ".format(*line.strip().split()), file=table_file)

        average_values = list(sum_values)
        average_values[0] = "AVERAGE"
        for i, item in enumerate(sum_values): 
            try:
                value = float(item)
                average_values[i] = round((value/completed_experiment_count),1)
                sum_values[i] = round(item,0) 

            except ValueError:
                pass

        print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9} \\\ \hline".format(*map(lambda x: x,sum_values)), file=table_file)
        print("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} & {8} & {9} \\\ \hline".format(*map(lambda x: x, average_values)), file=table_file)




