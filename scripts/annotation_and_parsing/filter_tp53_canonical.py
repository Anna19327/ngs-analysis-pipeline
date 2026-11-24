import csv
import sys

vcf_file = sys.argv[1] if len(sys.argv) > 1 else "tp53_vep.vcf"
canonical_only = "--all-transcripts" not in sys.argv[2:]

csq_fields = None
rows = []

with open(vcf_file, encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("##INFO=<ID=CSQ,"):
            description = line.split('Description="', 1)[1]
            description = description.split('">', 1)[0]

            fmt = description.split("Format: ", 1)[1]
            csq_fields = fmt.split("|")
            continue

        if line.startswith("#"):
            continue

        columns = line.rstrip("\n").split("\t")

        chrom, pos, variant_id, ref, alt = columns[:5]
        info = columns[7]

        csq_value = next(
            (
                item[4:] for item in info.split(";")
                if item.startswith("CSQ=")
            ),
            None,
        )

        if csq_value is None:
            continue

        for annotation in csq_value.split(","):
            values = annotation.split("|")

            record = dict(zip(csq_fields, values))

            if record.get("SYMBOL") != "TP53":
                continue

            if canonical_only and record.get("CANONICAL") != "YES":
                continue

            rows.append({
                "CHROM": chrom,
                "POS": pos,
                "ID": variant_id,
                "REF": ref,
                "ALT": alt,
                "CANONICAL": record.get("CANONICAL", ""),
                "Consequence": record.get("Consequence", ""),
                "IMPACT": record.get("IMPACT", ""),
                "Feature": record.get("Feature", ""),
                "BIOTYPE": record.get("BIOTYPE", ""),
                "EXON": record.get("EXON", ""),
                "HGVSc": record.get("HGVSc", ""),
                "HGVSp": record.get("HGVSp", ""),
                "Protein_position": record.get(
                    "Protein_position", ""
                ),
                "Amino_acids": record.get("Amino_acids", ""),
                "Codons": record.get("Codons", ""),
            })

if csq_fields is None:
    raise SystemExit("ERROR: CSQ definition not found in VCF header")

if canonical_only and "CANONICAL" not in csq_fields:
    print(
        "WARNING: 'CANONICAL' field is missing in CSQ header. "
        "Re-run VEP with --canonical or pass --all-transcripts.",
        file=sys.stderr,
    )

output_file = "tp53_consequences.tsv"

if rows:
    with open(output_file, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=rows[0].keys(),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

print("CSQ fields:", " | ".join(csq_fields))
print("TP53 canonical annotations:", len(rows))
print("Output:", output_file)
