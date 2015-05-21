#!/bin/bash -l

#SBATCH -A b2013169
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 1-10:00:00
#SBATCH -J gap_filling_evaluation

#module add bioinfo-tools
#module add bwa
#module add bowtie

for genome in Staphylococcus_aureus Rhodobacter_sphaeroides 
do
	for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
	do

		scaf_path='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Assembly/'$assembly'/genome.scf.fasta'

		for tool in Salmela # GapFiller_bwa GapFiller_bowtie #GapCloser  
		do

			echo 'Running' $tool 'on' $assembly 'on' $genome 

			tool_cfg='/proj/b2013169/private/kristoffer/gap_filling/'$tool'/'$genome'.config'
			tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/'
			filled_scf_name='filled_scfs.fasta'
			
			if [ ! -d $tool_out ]; then
			    mkdir -p $tool_out
			else
				rm -r $tool_out
				mkdir -p $tool_out
			fi

			if [ ! -f $scaf_path ]; then
			    echo "$scaf_path not found!"
			    echo "Skipping gap filling with $assembly on $genome"
			    continue
			fi


			if [ "$tool" = "Salmela" ]; then
				for n in 1 2 3 4 5 6 7 8 9 10
				do  
					echo "Running experiment $n" 
					tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/'
					if [ ! -d $tool_out ]; then
					    mkdir -p $tool_out
					else
						rm -r $tool_out
						mkdir -p $tool_out
					fi
					lib_path1='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/frag_merged.fastq.gz'				
					lib_path2='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/shortjump_merged.fastq.gz'
					
					/usr/bin/time -v GapTerminator -k 35 -filled $tool_out$filled_scf_name -scaffolds  $scaf_path -reads $lib_path1','$lib_path2 1> $tool_out'outsahlin.stdout' 2> $tool_out'outsahlin.stderr'
				done

			elif [ "$tool" = "GapCloser" ]; then
				/usr/bin/time -v GapCloser -a $scaf_path -b $tool_cfg -o $tool_out$filled_scf_name 1> $tool_out'outsahlin.stdout' 2> $tool_out'outsahlin.stderr'

			else 
				if [ ! -d $tool ]; then
					mkdir $tool
				else
					rm -r $tool
					mkdir $tool
				fi
				/usr/bin/time -v GapFiller  -s $scaf_path -l $tool_cfg -b $tool 1> $tool'/outsahlin.stdout' 2> $tool'/outsahlin.stderr'
				
				mv $tool $tool_out
				mv $tool_out$tool'/'$tool'.gapfilled.final.fa' $tool_out$filled_scf_name
				mv $tool_out$tool'/outsahlin.stdout' $tool_out
				mv $tool_out$tool'/outsahlin.stderr' $tool_out

			fi 
		done
	done
done


echo "FINISHED GAP CLOSING RUNS -- STARTING WITH EVALUATION"

for genome in Staphylococcus_aureus Rhodobacter_sphaeroides 
do
	genome_path='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/genome.fasta'
	for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
	do
		#ORIGINAL ASSEMBLY
		echo 'Evaluation of original assembly for' $assembly 'on' $genome

		original_scaffolds='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Assembly/'$assembly'/genome.scf.fasta'
		original_out='/proj/b2013169/private/data/gap_filling/original/'$genome'/'$assembly'/QUAST/'
		if [ ! -f $original_scaffolds ]; then
		    echo "$original_scaffolds not found!"
		    continue
		fi

		QUAST --strict-NA -s -o $original_out -R $genome_path $original_scaffolds 1> $original_out'outsahlin.stdout' 2> $original_out'outsahlin.stderr'


		#THE TOOLS

		filled_scf_name='filled_scfs.fasta'

		for tool in  Salmela #GapFiller_bwa GapFiller_bowtie GapCloser 
		do
			echo 'Evaluation of' $tool 'on' $assembly 'on' $genome

			tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/'
			tool_scaffolds='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/'$filled_scf_name		
			
			if [ "$tool" = "Salmela" ]; then
				for n in 1 2 3 4 5 6 7 8 9 10
				do  
					echo "Evaluating experiment $n" 
					tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/'
					tool_scaffolds='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/'$filled_scf_name		

					tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/QUAST/'
					if [ ! -d $tool_out ]; then
						mkdir $tool_out					    
					fi
					QUAST --strict-NA -s -o $tool_out -R $genome_path $tool_scaffolds 1> $tool_out'outsahlin.stdout' 2> $tool_out'outsahlin.stderr'
				done
			else
				
				QUAST --strict-NA -s -o $tool_out -R $genome_path $tool_scaffolds 1> $tool_out'outsahlin.stdout' 2> $tool_out'outsahlin.stderr'

			fi			

		done 
	done
done



echo "FINISHED EVALUATION OF GAP CLOSING -- PARSING EVALUATION RESULTS"


for genome in Staphylococcus_aureus Rhodobacter_sphaeroides 
do
	quast_csv='/home/kris/Work/gap_filling/main_output/quast_table_'$genome'.csv'
	performance_csv='/home/kris/Work/gap_filling/main_output/performance_table_'$genome'.csv'
	echo 'tool,assembly,misassemblies,seq-error-len,Unaligned\_length,NGA50,nr\_gaps,tot\_gap\_length' > $quast_csv
	echo 'User\_time,Wall\_clock\_time,memory' > $performance_csv

	for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
	do

		for tool in original GapCloser Salmela GapFiller_bwa GapFiller_bowtie
		do
			
			#parse results individually into intermediate files
			if [ "$tool" = "original" ]; then
				quast_ctg_rep_stdout='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/contigs_reports/contigs_report_genome.scf.stdout'
				quast_misass='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/contigs_reports/misassemblies_report.txt'
				quast_report='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/report.txt'
				gnu_time_mem='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/outsahlin.stderr'
				outpath='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/'
				python length_sum_misassembles_quast.py  $quast_ctg_rep_stdout	$quast_misass $quast_report	$gnu_time_mem $outpath --N 4000

				scaf_file='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Assembly/'$assembly'/genome.scf.fasta'
				gap_stats="$(python gap_statistics.py $scaf_file )"

				quality_results="$(cat $outpath'gap_filling_quality_eval.csv')"
				performance_results="$(cat $outpath'gap_filling_performance_eval.csv')"
				echo $tool','$assembly','$quality_results','$gap_stats >> $quast_csv
				echo $tool','$assembly','$performance_results >> $performance_csv
			
			elif [ "$tool" = "Salmela" ]; then
				for n in 1 2 3 4 5 6 7 8 9 10
				do  
					echo "Parsing results experiment $n" 
					quast_ctg_rep_stdout='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/QUAST/contigs_reports/contigs_report_filled_scfs.stdout'
					quast_misass='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/QUAST/contigs_reports/misassemblies_report.txt'
					quast_report='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/QUAST/report.txt'
					gnu_time_mem='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/outsahlin.stderr'
					outpath='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-'$n'/'
					python length_sum_misassembles_quast.py  $quast_ctg_rep_stdout	$quast_misass $quast_report	$gnu_time_mem $outpath --N 4000
	 			done
	 				
	 				outpath='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/'
	 				python average_results.py $outpath 10 

	 				scaf_file='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/experiment-1/filled_scfs.fasta'
					gap_stats="$(python gap_statistics.py $scaf_file )"
					#aggregate results into final file

					quality_results="$(cat $outpath'gap_filling_quality_eval.csv')"
					performance_results="$(cat $outpath'gap_filling_performance_eval.csv')"
					echo $tool','$assembly','$quality_results','$gap_stats >> $quast_csv
					echo $tool','$assembly','$performance_results >> $performance_csv
					#tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/-experiment-'$n'/QUAST'
				
			else
				quast_ctg_rep_stdout='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/contigs_reports/contigs_report_filled_scfs.stdout'
				quast_misass='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/contigs_reports/misassemblies_report.txt'
				quast_report='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/report.txt'
				gnu_time_mem='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/outsahlin.stderr'
				outpath='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/'
				python length_sum_misassembles_quast.py  $quast_ctg_rep_stdout	$quast_misass $quast_report	$gnu_time_mem $outpath --N 4000
 				
 				scaf_file='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/filled_scfs.fasta'
				gap_stats="$(python gap_statistics.py $scaf_file )"
				#aggregate results into final file

				quality_results="$(cat $outpath'gap_filling_quality_eval.csv')"
				performance_results="$(cat $outpath'gap_filling_performance_eval.csv')"
				echo $tool','$assembly','$quality_results','$gap_stats >> $quast_csv
				echo $tool','$assembly','$performance_results >> $performance_csv
			fi

			#quast_ctg_rep_stdout='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/contigs_reports/contigs_report_filled_scfs.stdout'


	
			#python format_results.py  $quast_ctg_rep_stdout	$quast_misass $quast_report	$gnu_time_mem $outpath

			# if [ "$tool" = "original" ]; then
			# 	scaf_file='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Assembly/'$assembly'/genome.scf.fasta'
			# else
 		# 		scaf_file='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/filled_scfs.fasta'
			# fi
			# gap_stats="$(python gap_statistics.py $scaf_file )"
			# #aggregate results into final file

			# quality_results="$(cat $outpath'gap_filling_quality_eval.csv')"
			# performance_results="$(cat $outpath'gap_filling_performance_eval.csv')"
			# echo $tool','$assembly','$quality_results','$gap_stats >> $quast_csv
			# echo $tool','$assembly','$performance_results >> $performance_csv
		done
	done
done

python make_csv_for_latex.py -q '/home/kris/Work/gap_filling/main_output/quast_table_Staphylococcus_aureus.csv' -o '/home/kris/Work/gap_filling/main_output/saur_quast.csv'
python make_csv_for_latex.py -q '/home/kris/Work/gap_filling/main_output/quast_table_Rhodobacter_sphaeroides.csv' -o '/home/kris/Work/gap_filling/main_output/rhodo_quast.csv'
python make_csv_for_latex.py -p '/home/kris/Work/gap_filling/main_output/performance_table_Staphylococcus_aureus.csv' -o '/home/kris/Work/gap_filling/main_output/saur_perf.csv'
python make_csv_for_latex.py -p '/home/kris/Work/gap_filling/main_output/performance_table_Rhodobacter_sphaeroides.csv' -o '/home/kris/Work/gap_filling/main_output/rhodo_perf.csv'

