echo "FINISHED EVALUATION OF GAP CLOSING -- PARSING EVALUATION RESULTS"


for genome in Staphylococcus_aureus Rhodobacter_sphaeroides 
do
	quast_csv='/home/kris/Work/gap_filling/main_output/quast_table_'$genome'.csv'
	performance_csv='/home/kris/Work/gap_filling/main_output/performance_table_'$genome'.csv'
	echo 'tool,assembly,#mismatches,#indels,#short_indels,#long_indels,#misassemblies,#local_misassemblies,#Ns_per_100_kbp,#mismatches_per_100_kbp,#indels_per_100_kbp,nr_gaps,tot_gap_length' > $quast_csv
	echo 'User_time,Wall_clock_time,memory' > $performance_csv

	for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
	do

		for tool in original GapCloser Salmela #GapFiller_bwa GapFiller_bowtie
		do
			
			#parse results individually into intermediate files
			
			quast_ctg_rep_stdout='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/contigs_reports/contigs_report_filled_scfs.stdout'
			quast_misass='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/contigs_reports/misassemblies_report.txt'
			quast_report='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST/report.txt'
			gnu_time_mem='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/outsahlin.stderr'
			outpath='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/'

			#assembly_quality='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/gap_filling_quality_eval.csv'
			#assembly_performance='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/gap_filling_performance_eval.csv'
			python format_results.py  $quast_ctg_rep_stdout	$quast_misass $quast_report	$gnu_time_mem $outpath

			if [ "$tool" = "original" ]; then
				scaf_file='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Assembly/'$assembly'/genome.scf.fasta'
			else
 				scaf_file='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/filled_scfs.fasta'
			fi
			gap_stats="$(python gap_statistics.py $scaf_file )"
			#aggregate results into final file

			quality_results="$(cat $outpath'gap_filling_quality_eval.csv')"
			performance_results="$(cat $outpath'gap_filling_performance_eval.csv')"
			echo $tool','$assembly','$quality_results','$gap_stats >> $quast_csv
			echo $tool','$assembly','$performance_results >> $performance_csv
		done
	done
done
