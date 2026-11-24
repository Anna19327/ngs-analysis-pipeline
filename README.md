# NGS Analysis Pipeline

Scripts and workflows for Next-Generation Sequencing (NGS) variant filtering, annotation, and parsing.

### `scripts/annotation_and_parsing/filter_tp53_canonical.py`

Python script for parsing VCF files annotated with **Ensembl VEP** (`CSQ` in the `INFO` field), filtering variants for the **`TP53`** gene, and exporting results to a tab-separated table (`tp53_consequences.tsv`).

#### Features
* Parses `CSQ` field definitions from the `##INFO=<ID=CSQ,...>` header line and maps each annotation value to its corresponding field name.
* Filters annotations to `SYMBOL == "TP53"`.
* Retains canonical transcripts (`CANONICAL == "YES"`) by default, or includes all *TP53* transcripts when `--all-transcripts` is passed.
* Prints warnings to `stderr` if the `CSQ` header definition is missing or `CANONICAL` is absent while canonical filtering is active.
* Writes filtered records to `tp53_consequences.tsv` (`CHROM`, `POS`, `ID`, `REF`, `ALT`, `CANONICAL`, `Consequence`, `IMPACT`, `Feature`, `BIOTYPE`, `EXON`, `HGVSc`, `HGVSp`, `Protein_position`, `Amino_acids`, `Codons`).

#### Usage
```bash
python scripts/annotation_and_parsing/filter_tp53_canonical.py path/to/annotated.vcf
python scripts/annotation_and_parsing/filter_tp53_canonical.py path/to/annotated.vcf --all-transcripts
```
