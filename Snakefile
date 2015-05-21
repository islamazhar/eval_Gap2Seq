"""
Submit this job on uppmax as:
    snakemake --debug --keep-going -j 999 --cluster "sbatch -A {params.account} -p {params.partition} -n {params.n}  -t {params.runtime} -C {params.memsize} -J {params.jobname} --mail-type={params.mail_type} --mail-user={params.mail}"
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



def latex_tables(wildcards):
  input_= []

  for experiment in config["DATASETS"]:
    input_.append(config["OUTBASE"]+"performance_table_{0}_{1}.tex".format(experiment, contamine))
    input_.append(config["OUTBASE"]+"percentage_quality_table_{0}_{1}.tex".format(experiment, contamine))
    input_.append(config["OUTBASE"]+"quality_table_{0}_{1}.tex".format(experiment, contamine))

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
    input: latex_tables,
            #expand(config["OUTBASE"]+"{experiment}/performance_table.tex",  experiment=config["EXPERIMENTS"]),
            #expand(config["OUTBASE"]+"{experiment}/quality_table.tex",  experiment=config["EXPERIMENTS"]), 
            #lambda wildcards: expand(config["INBASE"]+"{realdata}/{assembler}_mapped.bam", realdata=config["REALDATA"], assembler=config["REALDATA"]["ASSEMBLERS"])
            contaminated_mapped_files,
            reversed_read_files
    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

rule test:
    input: config["OUTBASE"]+"quality_table_sim_30.tex",
           config["OUTBASE"]+"percentage_quality_table_sim_30.tex",
           expand(config["OUTBASE"]+"sim/30/quast_{scaffolder}_NA.csv",scaffolder=config["SCAFFOLDERS"])
    params: 
        runtime="15:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname="all",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]

rule eval_scaffolders:
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

rule run_scaffolders:
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

rule create_experiment_indata:
    input: contaminated_mapped_files,
            reversed_read_files
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
            reads = lambda wildcards: config[wildcards.dataset]["MERGED_READS"],
    output: fasta= config["OUTBASE"]+"{dataset}/GAP2SEQ_{assembler}.fa",
            stderr= config["OUTBASE"]+"{dataset}/GAP2SEQ_{assembler}.stderr"
    log: config["OUTBASE"]+"{dataset}/GAP2SEQ.log
    params: 
        runtime = lambda wildcards: config["SBATCH"][wildcards.dataset]["gap2seq_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.dataset]["memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.dataset]["small_n"],
        path = config["gap2seq_rules"]["path"],
        python=config["PYTHON2"],
        jobname = "{dataset}_{contamine}_{assembler}"+"_GAP2SEQ",
        account = config["SBATCH"]["ACCOUNT"],
        mail = config["SBATCH"]["MAIL"],
        mail_type = config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        stdout = config["OUTBASE"]+"{0}/GAP2SEQ_{1}.stdout".format(wildcards.dataset, wildcards.assembler)
        shell("{time} Gap2Seq -filled {output.fasta} -scaffolds {input.scaffolds} -reads {input.reads} 1>  {stdout} 2> {output.stderr} ")
        shell("cat {stdout} {output.stderr} > {log}")


rule ORIGINAL:
    input:  ctgs = lambda wildcards: config[wildcards.experiment][wildcards.contamine][wildcards.assembler]["CONTIGS"],
    output: fasta= config["OUTBASE"]+"{experiment}/{contamine}/ORIGINAL_{assembler}.fa"
    params: 
        runtime="10:00",
        memsize = "mem128GB",
        partition = "core",
        n = "1",
        jobname = "{experiment}_{assembler}_{contamine}"+"_ORIGINAL",
        account=config["SBATCH"]["ACCOUNT"],
        mail=config["SBATCH"]["MAIL"],
        mail_type=config["SBATCH"]["MAIL_TYPE"]
    run:
        time = config["GNUTIME"]
        shell("{time} cp {input.ctgs} {output.fasta} ")

rule QUAST:
   input: scaffolds=config["OUTBASE"]+"{experiment}/{contamine}/{scaffolder}_{assembler}.fa",
            ref = lambda wildcards: config[wildcards.experiment]["REF"]
   output: nice_format=config["OUTBASE"]+"{experiment}/{contamine}/quast_{scaffolder}_{assembler}.csv" 
   params: 
       runtime=lambda wildcards: config["SBATCH"][wildcards.experiment]["quast_time"],
        memsize = lambda wildcards: config["SBATCH"][wildcards.experiment]["small_memsize"],
        partition = lambda wildcards: config["SBATCH"][wildcards.experiment]["small_partition"],
        n = lambda wildcards: config["SBATCH"][wildcards.experiment]["small_n"],
       jobname="quast_{experiment}_{contamine}_{scaffolder}",
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
       outpath="/tmp/{0}_{1}_{2}_{3}/QUAST/".format(wildcards.experiment, wildcards.contamine, wildcards.scaffolder, wildcards.assembler)
       shell(" {python} {path}quast.py -R  {input.ref} -o {outpath} -s --no-plots {input.scaffolds} ") 
       misassmblies, N50, NA50, tot_length, ESIZE_ASSEMBLY, ESIZE_GENOME, CORR_ESIZE_ASSEMBLY, CORR_ESIZE_GENOME = parse_quast(outpath+"report.txt")
       #e_size = get_esize(input.scaffolds)
       print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}".format(wildcards.assembler, wildcards.scaffolder, tot_length, N50, misassmblies,  NA50, ESIZE_ASSEMBLY, ESIZE_GENOME, CORR_ESIZE_ASSEMBLY, CORR_ESIZE_GENOME), file=open(output.nice_format, 'w'))    



rule time_and_mem:
   input:  stderr=config["OUTBASE"]+"{experiment}/{contamine}/{scaffolder}_{assembler}.stderr"
   output: outfile=config["OUTBASE"]+"{experiment}/{contamine}/{scaffolder}_{assembler}_time_and_mem.txt"
   params: 
       runtime="15:00",
       memsize = "'mem128GB|mem256GB|mem512GB'",
       partition = "core",
       n = "1",
       jobname="{experiment}_{assembler}_{contamine}"+"_time_and_mem",
       account=config["SBATCH"]["ACCOUNT"],
       mail=config["SBATCH"]["MAIL"],
       mail_type=config["SBATCH"]["MAIL_TYPE"]
   run:
       usertime, wallclocktime, memory_gb =  parse_gnu_time(input.stderr)
       print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(wildcards.experiment, wildcards.assembler, wildcards.scaffolder, usertime, wallclocktime, memory_gb), file=open(output.outfile, 'w') )
       

rule performace_latex_table:
   input: files=performance_input
   output: table=config["OUTBASE"]+"performance_table_{experiment}_{contamine}.tex"
   params: 
       runtime="15:00",
       memsize = "'mem128GB|mem256GB|mem512GB'",
       partition = "core",
       n = "1",
       jobname="performace_latex_table_{experiment}_{contamine}",
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
   output: table=config["OUTBASE"]+"quality_table_{experiment}_{contamine}.tex"
   params: 
       runtime="15:00",
       memsize = "'mem128GB|mem256GB|mem512GB'",
       partition = "core",
       n = "1",
       jobname="quality_latex_table_{experiment}_{contamine}",
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



rule quality_latex_table_percentages:
   input: files=quality_input 
   output: table=config["OUTBASE"]+"percentage_quality_table_{experiment}_{contamine}.tex"
   params: 
       runtime="15:00",
       memsize = "'mem128GB|mem256GB|mem512GB'",
       partition = "core",
       n = "1",
       jobname="percenatge_quality_latex_table_{experiment}_{contamine}",
       account=config["SBATCH"]["ACCOUNT"],
       mail=config["SBATCH"]["MAIL"],
       mail_type=config["SBATCH"]["MAIL_TYPE"]
   run:
        # get original assembly metrics to compute perceantages
        assembly_dict = {}
        for file_ in input.files:
            line=open(file_,'r').readlines()[0]  
            vals = line.strip().split()
            assembler, scaffolder, totsize, NG50, misassm, NGA50, E_a, E_g, EA_a, EA_g = vals
            totsize, NG50, misassm, NGA50, E_a, E_g, EA_a, EA_g = map(lambda x : float(x), [totsize, NG50, misassm, NGA50, E_a, E_g, EA_a, EA_g] )
            if scaffolder == "ORIGINAL":
                assembly_dict[assembler] = {'assm-size' : totsize , 'NG50' : NG50 ,
                                         'misassm' :  misassm, 'NGA50' : NGA50, 'E_a' : E_a ,
                                          'E_g' : E_g, 'EA_A':EA_a, 'EA_g' : EA_g }

        table_file = open(output.table, 'w')
        print("{0} & {1} & {2} & {3} & {4}  \\\ \hline".format('assembly', 'tool', 'size-diff', 'misassm', '$EA_s/EA_c$'), file=table_file)

        prev_scaffolder = -1

        for file_ in input.files:
            line=open(file_,'r').readlines()[0]  
            vals = line.strip().split()
            assembler, scaffolder, totsize, NG50, misassm, NGA50, E_a, E_g, EA_a, EA_g = vals

            totsize_perc = round(int(totsize) / assembly_dict[assembler]['assm-size'] , 3)
            misassm_diff = int(misassm) - int(assembly_dict[assembler]['misassm'])
            EA_g_perc = round(float(EA_g) / assembly_dict[assembler]['EA_g'] , 1)


            if prev_scaffolder == -1:
                sum_values = ["TOTAL", scaffolder, 0, 0, 0]
                completed_experiment_count = 0

            if scaffolder != prev_scaffolder and prev_scaffolder != -1:
                average_values = list(sum_values)
                average_values[0] = "AVERAGE"
                for i, item in enumerate(sum_values): 
                    try:
                        value = float(item)
                        average_values[i] = round((value/completed_experiment_count),3)
                        sum_values[i] = round(item,3) 
                    except ValueError:
                        pass

                print("{0} & {1} & {2} & {3} & {4}  \\\ \hline".format(*map(lambda x: x,sum_values)), file=table_file)
                print("{0} & {1} & {2} & {3} & {4}  \\\ \hline".format(*map(lambda x: x, average_values)), file=table_file)

                assembler, scaffolder, totsize, NG50, misassm, NGA50, E_a, E_g, EA_a, EA_g = vals
                sum_values = ["TOTAL", scaffolder, 0, 0, 0]
                completed_experiment_count = 1
                for i, item in enumerate([assembler, scaffolder, totsize_perc, misassm_diff, EA_g_perc]): 
                    try:
                        value = float(item)
                        sum_values[i] = value
                    except ValueError:
                        pass

            else:
                completed_experiment_count += 1
                for i, item in enumerate([assembler, scaffolder, totsize_perc, misassm_diff, EA_g_perc]): 
                    try:
                        value = float(item)
                        sum_values[i] += value
                    except ValueError:
                        pass

            prev_scaffolder = scaffolder


            print("{0} & {1} & {2} & {3} & {4} \\\ ".format(assembler, scaffolder, totsize_perc, misassm_diff,EA_g_perc), file=table_file)

        average_values = list(sum_values)
        average_values[0] = "AVERAGE"
        for i, item in enumerate(sum_values): 
            try:
                value = float(item)
                average_values[i] = round(value/completed_experiment_count, 3)
                sum_values[i] = round(item,3) 

            except ValueError:
                pass

        print("{0} & {1} & {2} & {3} & {4} \\\ \hline".format(*map(lambda x: x,sum_values)), file=table_file)
        print("{0} & {1} & {2} & {3} & {4} \\\ \hline".format(*map(lambda x: x, average_values)), file=table_file)




rule CONTAMINE:
  input: reads1=lambda wildcards: config[wildcards.realdata]["ORIGINAL_READ1"],
          reads2=lambda wildcards: config[wildcards.realdata]["ORIGINAL_READ2"]
  output: file1=config["INBASE"]+"{realdata}/{contamine}/frag1.fastq",
	        file2=config["INBASE"]+"{realdata}/{contamine}/frag2.fastq"
  params:
     runtime="45:00",
     memsize = "'mem128GB|mem256GB|mem512GB'",
     partition = "core",
     n = "2",
     jobname="{realdata}"+"_contamine_{contamine}",
     account=config["SBATCH"]["ACCOUNT"],
     mail=config["SBATCH"]["MAIL"],
     mail_type=config["SBATCH"]["MAIL_TYPE"],
      ref=lambda wildcards: config[wildcards.realdata]["REF"],
      length=lambda wildcards: config[wildcards.realdata]["READ_LENGTH"]
  run:
    if wildcards.contamine == "15":
        mean = 300
        sd=30
        c= 0.15
    elif wildcards.contamine == "40":
        mean = 400
        sd=40
        c=0.40
    elif wildcards.contamine == "0":
        mean = 300
        sd=30
        c=0.00001       
    """    
        Receipt for simulate contaimination reads for real libraries:
        1. count the number of reads in the original file = n
        2. n_c =  n_r*c_ratio / (1- c_ratio) where c_ratio is the contamine ratio, n_c is the number of contemned reads
        3. simulate k contaminated reads from 300,30 or 400,40. With the 
            same average length as original reads in fastq.zg format
        4. cat original_reads.gz sim_reads.gz > reads.fastq.gz
    """    
    nr_readpairs = int(list(shell("zcat {input.reads1} | wc -l ", iterable=True))[0]) /4
    nr_contamine_reads = int( (nr_readpairs*c)/(1 - c) )
    sim_path = "/tmp/{0}.frag".format(wildcards.realdata)
    shell("generate_library {params.ref} {nr_contamine_reads} {params.length} {sim_path} {mean} {sd}")
    shell("zcat -f {input.reads1} {sim_path}1.fastq > {output.file1} ")
    shell("zcat -f {input.reads2} {sim_path}2.fastq > {output.file2} ")

rule MAP:
  input: file1=config["INBASE"]+"{realdata}/{contamine}/frag1.fastq",
          file2=config["INBASE"]+"{realdata}/{contamine}/frag2.fastq",
          contigs=map_input
  output: bam=config["INBASE"]+"{realdata}/{contamine}/{assembly}_mapped.bam"
  params:
     runtime=lambda wildcards: config["SBATCH"][wildcards.realdata]["besst_time"],
     memsize = "'mem128GB|mem256GB|mem512GB'",
     partition = lambda wildcards: config["SBATCH"][wildcards.realdata]["partition"],
     n = lambda wildcards: config["SBATCH"][wildcards.realdata]["n"],
     prefix=config["INBASE"]+"{realdata}/{contamine}/{assembly}_mapped",
     jobname="{realdata}_{contamine}_map_{assembly}",
     account=config["SBATCH"]["ACCOUNT"],
     mail=config["SBATCH"]["MAIL"],
     mail_type=config["SBATCH"]["MAIL_TYPE"]
  run:
    if wildcards.realdata == "staph":
      shell("map_reads -sort --bowtie2 {input.file1} {input.file2} {input.contigs} {params.prefix}")
    else:
      shell("map_reads -sort --bowtie2 --local {input.file1} {input.file2} {input.contigs} {params.prefix}")

    print("map_reads -sort --bowtie2 {0} {1} {2} {3}".format(input.file1, input.file2, input.contigs, params.prefix))

rule REVCOMP:
    input:  reads1 = lambda wildcards: config[wildcards.experiment][wildcards.contamine]["READ1"],
            reads2 = lambda wildcards: config[wildcards.experiment][wildcards.contamine]["READ2"]
    output: reads1=config["INBASE"]+"{experiment}/{contamine}/frag1_fr.fastq",
            reads2=config["INBASE"]+"{experiment}/{contamine}/frag2_fr.fastq"
    params:
     runtime="20:00",
     memsize = "'mem128GB|mem256GB|mem512GB'",
     partition = "core",
     n = "1",
     jobname="revcomp_{experiment}"+"_contamine_{contamine}",
     account=config["SBATCH"]["ACCOUNT"],
     mail=config["SBATCH"]["MAIL"],
     mail_type=config["SBATCH"]["MAIL_TYPE"]

    run:
        shell("revcompfastx {input.reads1} > {output.reads1} ")
        shell("revcompfastx {input.reads2} > {output.reads2} ")

# rule clean:
#     input:
#     output:
#     run:

# rule test:
#     input:
#     output:
#     run:
