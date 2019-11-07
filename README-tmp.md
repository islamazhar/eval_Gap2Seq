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
- Go to `tools config` folder inside `eval_Gap2Seq`. Edit the `Staph.config` file inside `GapCloser`,`GapFiller_bowtie`, `gapFiller_bwa`.
- Run `$ pyenv shell 3.4.1` and then `snakemake run_gapfillers`. Use my `Snakefile` from here https://gist.github.com/islamazhar/a88b983cd4b4f46fd3cead1e3a9255a2  
. The ouptputs are place inside`OUTBASE` folder
-  Run `pyenv install 2.7.6` and then `snakemake eval_gapfillers`. This should give us the results.

## Problems
- I could not run Gap2Seq I am getting the following error. https://paste.ubuntu.com/p/gJ7xwXDcBb/ (pool allocation failing line no 24).
That is way I have removed Gap2Seq from line 142 of `config.json`.
