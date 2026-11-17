#!/usr/bin/env python3
# =============================================================================
# export_dwc.py: Darwin Core tables for publishing the public-data detections
# on GBIF/OBIS (core: Occurrence; extension: DNA-derived data, Abarenkov et al.
# 2023, "Publishing DNA-derived data through biodiversity data platforms").
#
# One occurrence per confirmed ASV. Coordinates are the representative
# coordinates of the dataset given in public_datasets.tsv: replace them with
# per-sample coordinates where the original study provides them.
# Review every field before publishing.
# =============================================================================

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--detections", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--out-occurrence", required=True)
    ap.add_argument("--out-dna", required=True)
    args = ap.parse_args()

    det = pd.read_csv(args.detections, sep="\t", dtype=str)
    ds = pd.read_csv(args.datasets, sep="\t", comment="#", dtype=str).set_index("dataset")

    occ, dna = [], []
    for r in det.itertuples():
        oid = f"{r.dataset}:{r.asv}"
        occ.append({
            "occurrenceID": oid,
            "basisOfRecord": "MaterialSample",
            "scientificName": "Facetotecta",
            "kingdom": "Animalia",
            "phylum": "Arthropoda",
            "occurrenceStatus": "present",
            "decimalLatitude": r.lat,
            "decimalLongitude": r.lon,
            "locality": r.region,
            "associatedSequences": ds.loc[r.dataset, "runs"],
            "identificationRemarks": (f"BLAST top hit Facetotecta ({r.identity}% identity) against "
                                      f"MIDORI2/SILVA plus the voucher-linked Facetotecta library; "
                                      f"placed in Facetotecta (EPA-ng); MOTU {r.motu}"),
            "bibliographicCitation": ds.loc[r.dataset, "citation"],
        })
        dna.append({
            "occurrenceID": oid,
            "DNA_sequence": r.sequence,
            "target_gene": r.marker,
            "pcr_primer_forward": ds.loc[r.dataset, "fwd"],
            "pcr_primer_reverse": ds.loc[r.dataset, "rev"],
            "otu_db": "Facetotecta reference library (this study) + MIDORI2/SILVA",
            "otu_seq_comp_appr": "BLASTN top hit + EPA-ng placement + phylogenetic filter",
        })
    pd.DataFrame(occ).to_csv(args.out_occurrence, sep="\t", index=False)
    pd.DataFrame(dna).to_csv(args.out_dna, sep="\t", index=False)


if __name__ == "__main__":
    main()
