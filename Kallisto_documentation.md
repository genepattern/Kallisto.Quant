# Kallisto (v2.0)

This is a bare-bones implementation of Kallisto, intended to be used for RNA quantitation against a human index, Gencode release 37 (GRCh38.p13), or a mouse index, Gencode release M26 (GRCm39).

Module Author: Edwin Juarez

Contact: https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help

Algorithm Version: Kallisto 0.46.1

<!-- ## Summary
*To be added* -->

## References
**Used verion 0.46.1, which is the most recent as of 2021-06-08**
- https://pachterlab.github.io/kallisto/manual

<!-- ### Functionality yet to be implemented:
*To be added*
*To be added* -->

### Technical notes:
- Docker container for this module: [genepattern/kallisto:5.46.1](https://hub.docker.com/layers/genepattern/kallisto/5.46.1/images/sha256-3161a3bdd346f8ebcd677ba6dec9f829ac8e38740db639676bd94a7c2f6e6895?context=repo)
- GitHub repo for this module: https://github.com/genepattern/Kallisto/releases/tag/v2
- This module uses Sleuth (version 0.30.0) for the gene aggregation step (but does not compute differential expression): https://anaconda.org/bioconda/r-sleuth & https://pachterlab.github.io/sleuth/manual

## Parameters

#### fastq_files
- Default: [blank]
- Required: Yes
- Description: The FASTQ (or fastq.gz) to be quantitated and translated into RNA seq transcript and gene counts.

#### transcriptome
- Default: Human
- Required: Yes
- Description: Which transcriptome to use for the pseudoalignment (Human or Mouse).

#### output_filename
- Default: 
- Required: Yes
- Description: 

#### extra_commands
- Default: [blank]
- Required: No
- Description: Write here what other parameters you want to pass to Kallisto, e.g., for single-end data, you can write "--single --fragment-length=200 --sd=30". This string will be passed to kallisto when it is called.


## Output Files
- module_log.txt: A list of the commands ran by the module.
- RNASeq_quant/<output_filename>_abundance.h5: transcript-level abundance. Use this file if you inted to use tximport.
- RNASeq_quant/<output_filename>_abundance.tsv: transcript-level abundance (same as the h5 file).
- RNASeq_quant/<output_filename>_normalized_gene_level.csv: normalized gene-level abundance.
- RNASeq_quant/<output_filename>_raw_estimated_counts.csv: raw gene-level estimated counts.
- RNASeq_quant/<output_filename>_transcript_level.csv: transcript-level abundance.
- RNASeq_quant/run_info.json: kallisto's parameters used.
- stdout.txt: A list of non-essential messages printed by Kallisto, this may be helpful to debug should any errors occur.
- stderr.txt: If there were errors in completing the job, they'd appear here. Note that the mere presence of a file with this name does not imply that there were errors. Here is an example of the contents of the stderr.txt file from our test run:
```
[quant] fragment length distribution is truncated gaussian with mean = 200, sd = 30
[index] k-mer length: 31
[index] number of targets: 234,485
[index] number of k-mers: 143,310,200
[index] number of equivalence classes: 1,002,767
[quant] running in single-end mode
[quant] will process file 1: /opt/gpbeta_2/gp_home/users/.cache/uploads/cache/datasets.genepattern.org/data/module_support_files/Kallisto/test_data/SRR1515119_50k.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 50,000 reads, 192 reads pseudoaligned
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 501 rounds
[bstrp] running EM for the bootstrap: 1
[bstrp] running EM for the bootstrap: 2

'gene_mode' is TRUE. Sleuth will do counts aggregation at the gene level for downstream normalization, transformation, and modeling steps, as well as for plotting and results.
reading in kallisto results
dropping unused factor levels
.
normalizing est_counts
6 targets passed the filter
normalizing tpm
merging in metadata
aggregating by column: gene_name
6 genes passed the filter
summarizing bootstraps
.
Warning messages:
1: In sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column = "gene_name",  :
  There is only one sample present, but you also provided a model. The model will be set to NULL to prevent downstream errors.
The sample can be viewed using sleuth_live after preparation, but you need more than one sample to run the other aspects of Sleuth.
2: `select_()` is deprecated as of dplyr 0.7.0.
Please use `select()` instead.
This warning is displayed once every 8 hours.
Call `lifecycle::last_warnings()` to see where this warning was generated. 
3: In data.table::melt(bs_df, id.vars = "bootstrap_num", variable.name = "target_id",  :
  The melt generic in data.table has been passed a data.frame and will attempt to redirect to the relevant reshape2 method; please note that reshape2 is deprecated, and this redirection is now deprecated as well. To continue using melt methods from reshape2 while both libraries are attached, e.g. melt.list, you can prepend the namespace like reshape2::melt(bs_df). In the next version, this warning will become an error.
4: In data.table::dcast(scaled_bs, sample ~ target_id, value.var = "scaled_reads_per_base") :
  The dcast generic in data.table has been passed a data.frame and will attempt to redirect to the reshape2::dcast; please note that reshape2 is deprecated, and this redirection is now deprecated as well. Please do this redirection yourself like reshape2::dcast(scaled_bs). In the next version, this warning will become an error.
Warning message:
In check_quant_mode(obj, which_units) :
  your sleuth object is in gene mode, but you selected 'est_counts'. Selecting 'scaled_reads_per_base'...
```

## License

Kallisto itself and this GenePattern module are distributed under a modified BSD license. This module's license is available at https://raw.githubusercontent.com/genepattern/Kallisto/develop/LICENSE


Version Comments
Version	Release Date	Description
2	2021-06-08	Adding Mouse index and upgrading Human index
1	2019-04-22	Initial release of Kallisto
