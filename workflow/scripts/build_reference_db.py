#!/usr/bin/env python3
# =============================================================================
# build_reference_db.py: build reference FASTA + taxonomy tables for BLAST.
#
# Two modes:
#   public  : reformat MIDORI2 (COX1, lrRNA) or SILVA (18S) into
#             - a FASTA with short IDs (pub0000001 ...), because makeblastdb
#               -parse_seqids rejects long IDs
#             - a taxonomy table: id <tab> lineage <tab> original header ID
#             SILVA sequences (RNA) are converted to DNA and only eukaryotes
#             are kept.
#   augment : add the curated Sanger library (IDs "SAN-...") to a public
#             reference, with the lineage Eukaryota;Arthropoda;Thecostraca;Facetotecta
#
# Usage:
#   build_reference_db.py public --input MIDORI.fasta --format midori \
#       --out-fasta X.fasta --out-tax X.tax.tsv
#   build_reference_db.py augment --input X.fasta --tax X.tax.tsv \
#       --library sanger.fasta --out-fasta Y.fasta --out-tax Y.tax.tsv
# =============================================================================

import argparse
import re
import sys

from yfutils import read_fasta, write_fasta

# Lineage given to every sequence of our Sanger library.
FACET_LINEAGE = "Eukaryota;Arthropoda;Thecostraca;Facetotecta"


def clean_rank(name):
    """MIDORI ranks look like 'Arthropoda_6656' (name_taxid). Drop the taxid."""
    return re.sub(r"_\d+$", "", name.strip())


def parse_midori(desc):
    """Return (original_id, lineage) from a MIDORI2 header.
    MIDORI2 headers: <accession and coordinates>###<rank_taxid>;<rank_taxid>;..."""
    if "###" in desc:
        acc, lineage = desc.split("###", 1)
    else:
        parts = desc.split(None, 1)
        acc, lineage = parts[0], parts[1] if len(parts) > 1 else ""
    ranks = [clean_rank(r) for r in lineage.split(";") if r.strip()]
    # Remove the artificial "root" rank.
    ranks = [r for r in ranks if r.lower() != "root"]
    return acc.split()[0], ";".join(ranks)


def parse_silva(desc):
    """Return (original_id, lineage) from a SILVA header:
    <accession.start.stop> <rank>;<rank>;...;<species>"""
    parts = desc.split(None, 1)
    return parts[0], parts[1].strip() if len(parts) > 1 else ""


def mode_public(args):
    """Reformat a public reference."""
    n_in, n_out = 0, 0
    records = []
    with open(args.out_tax, "w") as tax:
        tax.write("id\tlineage\toriginal_id\n")
        for _, desc, seq in read_fasta(args.input):
            n_in += 1
            if args.format == "midori":
                orig, lineage = parse_midori(desc)
            else:
                orig, lineage = parse_silva(desc)
                # Keep only eukaryotic 18S sequences.
                if not lineage.startswith("Eukaryota"):
                    continue
                # SILVA stores RNA: convert U to T.
                seq = seq.upper().replace("U", "T")
            seq = seq.upper().replace("-", "").replace(".", "")
            if len(seq) < 50:
                continue  # skip fragments too short to be useful
            n_out += 1
            new_id = f"pub{n_out:07d}"
            records.append((new_id, seq))
            tax.write(f"{new_id}\t{lineage}\t{orig}\n")
    write_fasta(records, args.out_fasta)
    print(f"read {n_in} sequences, wrote {n_out}", file=sys.stderr)


def mode_augment(args):
    """Append the Sanger library to a public reference."""
    records = [(n, s) for n, _, s in read_fasta(args.input)]
    lib = [(n, s.replace("-", "")) for n, _, s in read_fasta(args.library)]
    # Copy the public taxonomy, then add the library.
    with open(args.tax) as fin, open(args.out_tax, "w") as fout:
        fout.write(fin.read())
        for n, _ in lib:
            fout.write(f"{n}\t{FACET_LINEAGE}\t{n}\n")
    write_fasta(records + lib, args.out_fasta)
    print(f"public {len(records)} + library {len(lib)} sequences", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="mode", required=True)
    p = sub.add_parser("public")
    p.add_argument("--input", required=True)
    p.add_argument("--format", choices=["midori", "silva"], required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-tax", required=True)
    a = sub.add_parser("augment")
    a.add_argument("--input", required=True)
    a.add_argument("--tax", required=True)
    a.add_argument("--library", required=True)
    a.add_argument("--out-fasta", required=True)
    a.add_argument("--out-tax", required=True)
    args = ap.parse_args()
    if args.mode == "public":
        mode_public(args)
    else:
        mode_augment(args)


if __name__ == "__main__":
    main()
