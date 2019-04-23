# Kallisto (v1.0)

This is a bare-bones implementation of Kallisto, intended to be used for RNA quantitation against a human index (Homo_sapiens.GRCh38.95_kalllisto_index).

Author: Edwin Juarez

Contact: https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help

Algorithm Version: Kallisto 0.45

## Summary
*To be added*

## References
**Used verion 0.45, which is the most recent as of 2019-04-22**
- https://pachterlab.github.io/kallisto/manual
- https://anaconda.org/bioconda/kallisto
*More to be added*

### Functionality yet to be implemented:
*To be added*

### Technical notes:
*To be added*

## Parameters

#### extra_commands
- Default: [blank]
- Required: No
- Description: Write here what other parameters you want to pass to Kallisto, e.g., for single-end data, you can write "--single --fragment-length=200 --sd=30". This string will be passed to kallisto when it is called.

#### fastq_files
- Default: [blank]
- Required: Yes
- Description: The FASTQ (or fastq.gz) to be quantitated and translated into RNA seq transcript and gene counts.

## Output Files
- module_log.txt: A list of the commands ran by the module.
- RNASeq_quant/abundance.h5: transcript-level abundance.
- RNASeq_quant/abundance.tsv: transcript-level abundance (same as the h5 file).
- RNASeq_quant/gene_expression.csv: gene-level abundance.
- RNASeq_quant/run_info.json: kallisto's parameters used.


## License

Kallisto itself and this GenePattern module are distributed under a modified BSD license. This module's license is available at https://raw.githubusercontent.com/genepattern/Kallisto/develop/LICENSE

## Platform Dependencies
Task Type: RNASeq quantitation
CPU Type:
any

Operating System:
any

Language:
C++

Version Comments
Version	Release Date	Description
1	2019-04-22	Initial release of Kallisto
