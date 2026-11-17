#!/usr/bin/env python3
# =============================================================================
# yfutils.py: small helpers shared by the pipeline's Python scripts.
#
# - FASTA/FASTQ reading and writing (plain or gzip), without extra packages,
#   so every conda environment can use it
# - platform tags: each sequence ID starts with a tag and a dash:
#     SAN- Sanger reference, PB- PacBio, ILL- Illumina, PUB- public data,
#     OUT- outgroup
# - writing "metrics" files: key <tab> value, collected by summarize_numbers.py
# =============================================================================

import gzip
import os


def _open(path, mode="rt"):
    """Open a file, using gzip when the name ends in .gz."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path):
    """Yield (id, description, sequence) for each record of a FASTA file.
    id is the first word of the header; description is the full header."""
    name, desc, chunks = None, None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                # Return the previous record before starting a new one.
                if name is not None:
                    yield name, desc, "".join(chunks)
                desc = line[1:].strip()
                name = desc.split()[0] if desc else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, desc, "".join(chunks)


def read_fasta_dict(path):
    """Read a FASTA file into an ordered dict {id: sequence}."""
    return {n: s for n, _, s in read_fasta(path)}


def write_fasta(records, path, width=0):
    """Write (id, sequence) pairs to FASTA. width=0 writes each sequence on one line."""
    with _open(path, "wt") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            if width and width > 0:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


def read_fastq(path):
    """Yield (id, header, sequence, quality) for each FASTQ record."""
    with _open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n\r")
            fh.readline()  # the "+" line
            qual = fh.readline().rstrip("\n\r")
            header = header.rstrip("\n\r")[1:]
            yield header.split()[0], header, seq, qual


def write_fastq_record(fh, name, seq, qual):
    """Write one FASTQ record to an open file handle."""
    fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


# Reverse complement, including IUPAC ambiguity codes.
_COMP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn-?",
                      "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn-?")


def revcomp(seq):
    """Reverse complement of a DNA sequence."""
    return seq.translate(_COMP)[::-1]


def platform_of(seq_id):
    """Platform tag of a sequence ID (text before the first dash)."""
    tag = seq_id.split("-", 1)[0]
    return tag if tag in {"SAN", "PB", "ILL", "PUB", "OUT"} else "UNK"


def write_metrics(path, metrics):
    """Write a dict of metrics as key <tab> value lines."""
    if not path:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for k, v in metrics.items():
            fh.write(f"{k}\t{v}\n")


def read_membership(path):
    """Read a membership table (sequence <tab> MOTU) into a dict."""
    members = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "sequence":
                continue  # header line
            members[parts[0]] = parts[1]
    return members


def read_taxonomy(path):
    """Read a taxonomy table (id <tab> lineage [<tab> original id])."""
    tax = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] != "id":
                tax[parts[0]] = parts[1]
    return tax
