- Place `Gap2Seq`, `GapCloser`, `GapFiller`, and `Quast` in the `$HOME`. 
For gapfiller see this https://github.com/islamazhar/eval_Gap2Seq#preliminaries.

Add the following lines to bashrc file

    export Gap2Seq=$HOME/Gap-filling/eval-insertfill/Gap2Seq/build/Gap2Seq
    export GapCloser=$HOME/GapCloser-src-v1.12-r6/v1.12-r6/bin/GapCloser
    export GapFiller=$HOME/gapFiler/GapFiller
    
    
   Reload the bashrc file by running `source ~/.bashrc`
- cd to the `frag` and `short` reads of Staph from GAUGE data folder. In my case `cd /home/islamazhar/GAGE_Data/Staph/Data/original`
You need to merge the paired reads in a single file. 
    - `cat frag_1.fastq frag_2.fastq > frag_merged.fastq` ;  `gzip frag_merged.fastq`
    -  `cat shortjump_1.fastq shortjump_2.fastq > frag_merged.fastq` ;  `gzip shortjump_merged.fastq`
- Go to `config.jason` inside `eval_Gap2Seq` folder. Change the variables. I have kept only `Staph` in `"DATASETS"` variable.
- Go to `tools config` folder inside `eval_Gap2Seq`. Edit the `Staph.config` file inside `GapCloser`,`GapFiller_bowtie`, `gapFiller_bwa`, and `Salemla`.   
`Salmela` is actually `Gap2Seq`. `Salmela` is the last name of the first author of the paper. 
- Use my `Snakefile` from there https://github.com/islamazhar/eval_Gap2Seq/blob/master/Snakefile
- Run `snakemake run_gapfillers`. It should run and place the results inside `tools_configs` folders such as the result of `GapCloser` is in 
    `OUTBASE` folder
-  Run `pyenv install 2.7.6` and then `snakemake eval_gapfillers`. This should give us the results. However this commands fails as discussed in the following.

## Problems
- I could not run Gap2Seq I am getting the following error. https://paste.ubuntu.com/p/gJ7xwXDcBb/ (pool allocation failing line no 24).
That is way I have removed Gap2Seq from line 142 of `config.json`.

Solution:  
    - One way to run `Gap2Seq` I think is to run it on high configuration PC (UIU server). But I not too sure.
I run it on a short genome.. the one we tested our method on and it ran fine.
- I cound not run QUAST 2.3 because of some werid errors. I replaced Quast 2.3 with 5.0.2 and it worked. 
For this reason I tried to change the `correct_quast.py` file as well since it was written to 
parse quast result of 2.3 not 5.0.2. And since 2.3 and 5.0.2 have different reporting format
`correct_quast.py` is not being able to purse quast-5.0.2 reports. But I have been able to change it properly errors keep showing up.
 Hence the second snakemake `snakemake eval_gapfillers` command fails. 
    
    Solution:
    - If we can run quast-2.3 our problem will be solved   
    - Otherwise we need to change correct_quast.py :-(
