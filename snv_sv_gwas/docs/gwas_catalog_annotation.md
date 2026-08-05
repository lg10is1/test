# Local GWAS Catalog annotation

`08_annotate_gwas_catalog_scz_related.py` annotates an input table containing `Gene_refGene` by matching the local GWAS Catalog fields `MAPPED_GENE`, `MAPPED_TRAIT`, and `DISEASE/TRAIT`.

Supply an existing local associations zip with `--gwas-zip PATH` or `GWAS_CATALOG_ASSOC_ZIP`. The public script does not download data. Record the selected catalog release and checksum.

```bash
python3 scripts/08_annotate_gwas_catalog_scz_related.py \
  --input /path/to/lead_annotations.tsv \
  --gwas-zip /path/to/gwas_catalog_associations.zip \
  --out-dir /path/to/output
```

The script writes an annotated table, a gene summary, and the matching association rows. Trait-keyword defaults are defined in the script and should be reviewed against the study methods.
