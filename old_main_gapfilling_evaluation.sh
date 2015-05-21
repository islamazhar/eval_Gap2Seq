#!/bin/bash -l

#SBATCH -A b2013169
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J rhody_scaf_evaluation

module add bioinfo-tools
module add bwa
module add bowtie

for genome in Staphylococcus_aureus Rhodobacter_sphaeroides 
do
	for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
	do

		scaf_path='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Assembly/'$assembly'/genome.scf.fasta'

		lib1_path1='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/frag_1.fastq.gz'
		lib2_path2='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/frag_2.fastq.gz'
		lib1_path1='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/shortjump_1.fastq.gz'
		lib2_path2='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/shortjump_2.fastq.gz'

		lib_salmela_path='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/frag_merged.fastq.gz'

		GapFiller_cfg_bowtie='/proj/b2013169/private/kristoffer/gap_filling/GapFiller/'$genome'_bowtie.config'
		GapFiller_cfg_bwa='/proj/b2013169/private/kristoffer/gap_filling/GapFiller/'$genome'_bwa.config'
		GapCloser_cfg='/proj/b2013169/private/kristoffer/gap_filling/GapCloser/'$genome'.config'

		Gapfiller_out_folder='/proj/b2013169/private/data/gap_filling/GapFiller/'$genome'/'$assembly'/'
		GapCloser_out_folder='/proj/b2013169/private/data/gap_filling/GapCloser/'$genome'/'$assembly'/'
		Salmela_out_folder='/proj/b2013169/private/data/gap_filling/Salmela/'$genome'/'$assembly'/'

		 
		
		if [ ! -f $scaf_path ]; then
		    echo "$scaf_path not found!"
		    continue
		fi

		if [ ! -d 'bwa_out' ]; then
		    mkdir 'bwa_out'
		fi
		if [ ! -d 'bowtie_out' ]; then
		    mkdir 'bowtie_out'
		fi

		if [ ! -d $Gapfiller_out_folder ]; then
		    mkdir $Gapfiller_out_folder
		fi
		if [ ! -d $GapCloser_out_folder ]; then
		    mkdir $GapCloser_out_folder
		fi
		if [ ! -d $Salmela_out_folder ]; then
		    mkdir $Salmela_out_folder
		fi

		echo 'Running GapCloser on' $assembly 'on' $genome 
		/usr/bin/time -v GapCloser -a $scaf_path -b $GapCloser_cfg -o $GapCloser_out_folder'out' 1> $GapCloser_out_folder'outsahlin.stdout' 2> $GapCloser_out_folder'outsahlin.stderr'

		echo 'Running GapFiller on' $assembly 'on' $genome 'with aligner bwa'
		/usr/bin/time -v GapFiller  -s $scaf_path -l $GapFiller_cfg_bwa -b 'bwa_out' 1> 'bwa_out/outsahlin.stdout' 2> 'bwa_out/outsahlin.stderr'
		mv -f 'bwa_out' $Gapfiller_out_folder'bwa_out'

		echo 'Running GapFiller on' $assembly 'on' $genome 'with aligner bowtie'
		/usr/bin/time -v GapFiller  -s $scaf_path -l $GapFiller_cfg_bowtie -b 'bowtie_out' 1> 'bowtie_out/outsahlin.stdout' 2> 'bowtie_out/outsahlin.stderr'
		mv -f 'bowtie_out' $Gapfiller_out_folder'bowtie_out'

		echo 'Running Salmela on' $assembly 'on' $genome 
		/usr/bin/time -v GapTerminator -filled $Salmela_out_folder'out' -scaffolds  $scaf_path -reads $lib_salmela_path 1> $Salmela_out_folder'outsahlin.stdout' 2> $Salmela_out_folder'outsahlin.stderr'
	done
done

echo "FINISHED GAP CLOSING RUNS -- STARTING WITH EVALUATION"

for genome in Staphylococcus_aureus Rhodobacter_sphaeroides 
do
	genome_path='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Data.original/genome.fasta'
	for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
	do
		
		original_scaffolds='/proj/b2013169/private/data/gage.cbcb.umd.edu/data/'$genome'/Assembly/'$assembly'/genome.scf.fasta'
		original_out='/proj/b2013169/private/data/gap_filling/original/'$genome'/'$assembly'/QUAST'
		if [ ! -f $original_scaffolds ]; then
		    echo "$original_scaffolds not found!"
		    continue
		fi

		QUAST -s -o $original_out -R $genome_path $original_scaffolds
		
		for tool in GapFiller_bwa GapFiller_bowtie GapCloser Salmela
		do
			if [ "$tool" = "GapCloser" ]; then
				tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST'
				tool_scaffolds='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/out'
			elif  [ "$tool" = "GapFiller_bwa" ]; then
				tool_out='/proj/b2013169/private/data/gap_filling/GapFiller/'$genome'/'$assembly'/bwa/QUAST'
				tool_scaffolds='/proj/b2013169/private/data/gap_filling/GapFiller/'$genome'/'$assembly'/bwa_out/bwa_out.gapfilled.final.fa'
			elif  [ "$tool" = "GapFiller_bowtie" ]; then
				tool_out='/proj/b2013169/private/data/gap_filling/GapFiller/'$genome'/'$assembly'/bowtie/QUAST'
				tool_scaffolds='/proj/b2013169/private/data/gap_filling/GapFiller/'$genome'/'$assembly'/bowtie_out/bowtie_out.gapfilled.final.fa'
			elif [ "$tool" = "Salmela" ]; then
				tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/QUAST'
				tool_scaffolds='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/out'
			fi
			
			QUAST -s -o $tool_out -R $genome_path $tool_scaffolds
			
		done 
	done
done

echo "FINISHED EVALUATION OF GAP CLOSING -- PARSING QUAST RESULTS"



# echo "#Indels,#Inversions,#Translocations,#Relocations,#Gaps>=1000,MEANGAPERROR,#Contigs>N50,CorrectedN50,CorrectedE-size" > /proj/b2013169/private/data/BESST_paper/BESST/rhody_path_est/results.txt
# for scaffolder in BESST  #SSPACE SOPRA Opera
# do
# #echo $scaffolder
# for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
# do
# echo $assembly
# if [ "$scaffolder" = "BESST" ]; then
#         scaff_file="BESST_output/pass1/Scaffolds_pass1.fa"
# #elif  [ "$scaffolder" = "SSPACE" ]; then
# #       scaff_file="sspace_out/sspace_out.final.scaffolds.fasta"
# #elif [ "$scaffolder" = "SOPRA" ]; then
# #       scaff_file="scaffolds_h2.2_L150_w4.fasta"
# #elif [ "$scaffolder" = "Opera" ]; then
# #       scaff_file="scaffoldSeq.fasta"
# fi

# gage-validate "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody_no_score/$assembly/gage_eval/" "/proj/b2013169/private/data/bacgenomes/references/modified_references/rhody.reference.fasta" "/proj/b2013169/private/data/gage.cbcb.umd.edu/data/Rhodobacter_sphaeroides/Assembly/$assembly/genome.ctg.fasta"  "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody_no_score/$assembly/$scaff_file"
# echo -n "$scaffolder,$assembly," >> /proj/b2013169/private/data/BESST_paper/BESST/rhody_no_score/results.txt
# tail -n 1 "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody_path_est/$assembly/gage_eval/gage-results.eval"  >> /proj/b2013169/private/data/BESST_paper/BESST/rhody_no_score/results.txt
# done
# done



# #echo 'Supplementary file'

# #echo "Scaffolder,Assembly,Scaffolds,NG50,CorrNG50,CorrE-size,MisAssemblies,time" > supplemental_results.txt

# #for scaffolder in BESST SSPACE SOPRA Opera
# #do
# #for assembly in ABySS Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
# #do
# #echo -n "$scaffolder,$assembly," >> supplemental_results.txt
# #grep  "Total units:" "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody/$assembly/gage_eval/gage-results.eval" | tail -n 1 | cut -f 3 -d ' ' | tr -d '\n' >> supplemental_results.txt

# #echo -n ',' >> supplemental_results.txt

# #grep  "N50:" "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody/$assembly/gage_eval/gage-results.eval" | tail -n 1 | cut -f 2 -d ' ' | tr -d '\n' >> supplemental_results.txt

# #echo -n ',' >> supplemental_results.txt

# #tail -n 1 "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody/$assembly/gage_eval/gage-results.eval" | cut -f 8-9 -d ',' | tr -d '\n' >> supplemental_results.txt

# #echo -n ',' >> supplemental_results.txt

# #tail -n 1 "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody/$assembly/gage_eval/gage-results.eval" | cut -f 1-4 -d ',' | awk -F ',' '{s+=$1+$2+$3+$4} END {print s}' | tr -d '\n' >> supplemental_results.txt

# #echo -n ',' >> supplemental_results.txt
# #grep 'real'  "/proj/b2013169/private/data/BESST_paper/$scaffolder/rhody/$assembly/gage_eval/runtime.txt" | cut -f 2 >> supplemental_results.txt

# #done
# #done