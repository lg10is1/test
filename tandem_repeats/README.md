# Tandem repeats

Tandem-repeat logistic association test.

## Requirements

- R with `data.table` and `logistf`

## Usage

The TR genotype matrix is split into chunks. Run one chunk per invocation; the numeric argument is the chunk index:

```bash
Rscript tr_logistic_gwas.R 1
```

Expected inputs in the working directory:

- `tmp/dscz/<chunk>`: case haplotype TR dosage chunk
- `tmp/dcontrol/<chunk>`: control haplotype TR dosage chunk
- `lord.RData`: list of TR identifiers per chunk
- `sample.info`: sample table; column 2 is recoded to `SCZ`

For each TR, a logistic regression of disease status on dosage is fit. The same code was used for the SRS and LRS TR matrices.
