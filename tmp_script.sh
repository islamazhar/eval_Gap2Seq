for genome in Staphylococcus_aureus Rhodobacter_sphaeroides 
do
	for assembly in ABySS ABySS2 Allpaths-LG Bambus2 CABOG MSR-CA SGA SOAPdenovo Velvet
	do

		for tool in GapFiller_bwa GapFiller_bowtie #GapCloser Salmela 
		do

			echo 'Running' $tool 'on' $assembly 'on' $genome 
			tool_out='/proj/b2013169/private/data/gap_filling/'$tool'/'$genome'/'$assembly'/'

			if [ ! -f $scaf_path ]; then
			    echo "$scaf_path not found!"
			    echo "Skipping gap filling with $assembly on $genome"
			    continue
			fi

				mv $tool_out$tool'/outsahlin.stdout' $tool_out
				mv $tool_out$tool'/outsahlin.stderr' $tool_out

		 
		done
	done
done
