# Path placeholders

The full workflow scripts intentionally use placeholder paths instead of private project paths.

Common placeholders:

| Placeholder | Meaning |
|---|---|
| `/path/to/EOSCZ_PROJECT` | Main project root containing GWAS, Figure outputs, TGS callsets, sample metadata |
| `/path/to/SCZ` | Methylation project root used by meQTL scripts |
| `/path/to/shared_resources` | Shared reference files, BAM directories, genome FASTA/GTF resources |
| `/path/to/local/figure` | Local directory for manually prepared annotation/input tables |
| `/path/to/annovar` | ANNOVAR installation directory |
| `/path/to/igv.sh` | IGV executable used to render optional IGV snapshots |

PanGenie public example:

```text
/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas/set00/gwas/SCZ.mlm.ngspc.fastGWA
/path/to/EOSCZ_PROJECT/TGS_callset/Pangenie_v3/06.gwas/set00/NGS.QCsite.QCind.{bed,bim,fam}
```

DeepVariant/Paragraph examples:

```text
/path/to/EOSCZ_PROJECT/GWAS/Deepvariant_paragraph/deepvar_gwas/deepvar/03_gwas/SCZ.deepvar.mlm.geno0.1.maf0.01.fastGWA
/path/to/EOSCZ_PROJECT/GWAS/Deepvariant_paragraph/deepvar_gwas/paragraph_test/03_gwas/SCZ.paragraph_test.MODEL.mlm.geno0.1.maf0.01.fastGWA
/path/to/EOSCZ_PROJECT/GWAS/Deepvariant_paragraph/chr_all2.strict_step2_genimi.common_samples.merged.{bed,bim,fam}
```

