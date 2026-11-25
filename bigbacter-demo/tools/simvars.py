#!/usr/bin/env python3
"""
Simulate SNP variants of an assembly and generate paired FASTQ files.

Pipeline, built on existing tools:

  1. CFSAN snp-mutator (`snpmutator`) introduces an exact number of random
     substitutions into the input assembly, once per variant. The SNP count for
     each variant is drawn uniformly from --min-snps..--max-snps, so the set
     spans a range of distances from the reference.
  2. wgsim (default) or InSilicoSeq (`iss`) simulates paired Illumina reads
     from each variant genome at --depth fold coverage.
  3. Realized depth is measured back from the FASTQ files and reported, so the
     target depth is verified rather than assumed.

For each variant it keeps the mutated assembly (gzipped) alongside its reads,
plus a truth-set TSV of every introduced SNP, per-variant VCFs, and a run
manifest listing the SNP count, seed and realized depth for each variant.

Sample names never contain an underscore. Underscore is the field delimiter in
Illumina FASTQ filenames (sample_S1_L001_R1_001.fastq.gz), so a sample name
carrying one is split apart by anything that parses those names.

This is a standalone simulator: it knows nothing about any downstream pipeline
or sample metadata. The manifest carries the sample names and FASTQ paths if
you need to build an input sheet for something else.

Dependencies
    pip install snp-mutator            # snpmutator
    pip install InSilicoSeq            # iss, only for --simulator iss
    wgsim: conda install -c bioconda wgsim, or
           git clone https://github.com/lh3/wgsim && cc -O2 -o wgsim wgsim.c -lz -lm

Note on contigs: snp-mutator concatenates a multi-contig assembly into a
single sequence, so the reference written here (reference.fasta) is that
concatenation, and reads may span artificial contig junctions. For a draft
assembly with many short contigs, consider the junction artefact before
treating these reads as a faithful resequencing simulation.
"""

import argparse
import csv
import glob
import gzip
import os
import random
import shutil
import subprocess
import sys

# read length baked into each bundled InSilicoSeq error model
ISS_MODEL_READ_LENGTH = {
    "hiseq": 126,
    "novaseq": 151,
    "miseq": 301,
    "miseq-20": 301,
    "miseq-24": 301,
    "miseq-28": 301,
    "miseq-32": 301,
    "miseq-36": 301,
    "nextseq": 301,
}


def die(msg):
    sys.exit(f"error: {msg}")


def run(cmd, label):
    """Run a subprocess, surfacing its output only on failure."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(f"\n{label} failed: {' '.join(cmd)}\n")
        sys.stderr.write(proc.stdout[-4000:])
        sys.stderr.write(proc.stderr[-4000:])
        die(f"{label} exited {proc.returncode}")
    return proc


def read_fasta(path):
    """Minimal FASTA reader: [(id, description, sequence)]."""
    recs, name, desc, chunks = [], None, "", []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, desc, "".join(chunks)))
                head = line[1:].strip()
                parts = head.split(None, 1)
                name = parts[0] if parts else "seq"
                desc = parts[1] if len(parts) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        recs.append((name, desc, "".join(chunks)))
    return recs


def fastq_bases(path):
    """Total sequenced bases in a FASTQ, for measuring realized depth."""
    op = gzip.open if path.endswith(".gz") else open
    total = 0
    with op(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                total += len(line.strip())
    return total


def compress(path, keep=False):
    """gzip a file, preferring pigz/gzip binaries over the gzip module.

    The timestamp is suppressed (-n / mtime=0) so that two runs with the same
    seed produce byte-identical .gz files, not just identical contents."""
    out = path + ".gz"
    tool = shutil.which("pigz") or shutil.which("gzip")
    if tool:
        with open(out, "wb") as fh:
            proc = subprocess.run([tool, "-n", "-c", path], stdout=fh)
        if proc.returncode == 0:
            if not keep:
                os.remove(path)
            return out
        if os.path.exists(out):
            os.remove(out)
    with open(path, "rb") as src, gzip.GzipFile(out, "wb", mtime=0) as dst:
        shutil.copyfileobj(src, dst, length=1 << 20)
    if not keep:
        os.remove(path)
    return out


def parse_vcf(path):
    """Yield (chrom, pos, ref, alt) from a snp-mutator VCF."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                yield f[0], f[1], f[3], f[4]


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("assembly", help="input assembly FASTA")

    var = ap.add_argument_group("variant generation (snp-mutator)")
    var.add_argument("-n", "--variants", type=int, default=6,
                     help="number of variant genomes to create (default: 6)")
    var.add_argument("--min-snps", type=int, default=10,
                     help="low end of the per-variant SNP range (default: 10)")
    var.add_argument("--max-snps", type=int, default=200,
                     help="high end of the per-variant SNP range (default: 200)")
    var.add_argument("--insertions", type=int, default=0,
                     help="insertions per variant (default: 0, SNPs only)")
    var.add_argument("--deletions", type=int, default=0,
                     help="deletions per variant (default: 0, SNPs only)")
    var.add_argument("--seqid", default=None,
                     help="sequence ID for the mutated FASTA and VCF CHROM column. "
                          "Defaults to the first contig's ID; a value with spaces would "
                          "produce a malformed VCF.")

    rd = ap.add_argument_group("read simulation")
    rd.add_argument("--depth", type=float, default=50.0,
                    help="target fold coverage per variant (default: 50)")
    rd.add_argument("--simulator", choices=["wgsim", "iss"], default="wgsim",
                    help="wgsim is ~100x faster; iss gives empirical Illumina error "
                         "profiles (default: wgsim)")
    rd.add_argument("--read-length", type=int, default=150,
                    help="read length for wgsim (default: 150). With --simulator iss "
                         "the model's own read length is used instead.")
    rd.add_argument("--fragment-length", type=int, default=350, help="mean fragment size (default: 350)")
    rd.add_argument("--fragment-sd", type=int, default=50, help="fragment size SD (default: 50)")
    rd.add_argument("--error-rate", type=float, default=0.002,
                    help="wgsim per-base error rate (default: 0.002)")
    rd.add_argument("--iss-model", default="novaseq",
                    choices=sorted(ISS_MODEL_READ_LENGTH), help="iss error model (default: novaseq)")
    rd.add_argument("--iss-mode", choices=["kde", "basic", "perfect"], default="kde",
                    help="iss error mode; kde uses --iss-model (default: kde)")
    rd.add_argument("--threads", type=int, default=2, help="threads for iss (default: 2)")

    out = ap.add_argument_group("outputs")
    out.add_argument("-o", "--outdir", default="simulated", help="output directory (default: simulated)")
    out.add_argument("--prefix", default="sim", help="sample name prefix (default: sim)")
    out.add_argument("--seed", type=int, default=42, help="master RNG seed (default: 42)")
    out.add_argument("--no-gzip", action="store_true",
                     help="leave FASTQ files and mutated assemblies uncompressed")
    args = ap.parse_args()

    # ---- validate -------------------------------------------------------- #
    if not os.path.isfile(args.assembly):
        die(f"assembly not found: {args.assembly}")
    if args.min_snps < 0 or args.max_snps < args.min_snps:
        die("require 0 <= --min-snps <= --max-snps")
    if args.variants < 1:
        die("--variants must be at least 1")
    if args.depth <= 0:
        die("--depth must be positive")
    # Underscore separates fields in Illumina FASTQ filenames, so a sample name
    # containing one cannot be recovered by anything that parses those names.
    if "_" in args.prefix:
        die(f"--prefix must not contain an underscore: {args.prefix!r}")
    if any(c.isspace() for c in args.prefix) or "/" in args.prefix or "\\" in args.prefix:
        die(f"--prefix must not contain whitespace or path separators: {args.prefix!r}")

    if not shutil.which("snpmutator"):
        die("snpmutator not found on PATH. Install with: pip install snp-mutator")
    if not shutil.which(args.simulator):
        hint = ("conda install -c bioconda wgsim, or build from "
                "https://github.com/lh3/wgsim" if args.simulator == "wgsim"
                else "pip install InSilicoSeq")
        die(f"{args.simulator} not found on PATH. Install with: {hint}")

    recs = read_fasta(args.assembly)
    if not recs:
        die(f"no sequences parsed from {args.assembly}")
    genome_len = sum(len(s) for _, _, s in recs)
    if genome_len == 0:
        die("input assembly has zero total sequence length")
    seqid = args.seqid or recs[0][0]
    if any(c.isspace() for c in seqid):
        die(f"--seqid must not contain whitespace: {seqid!r}")

    read_len = (ISS_MODEL_READ_LENGTH[args.iss_model] if args.simulator == "iss"
                else args.read_length)
    if args.simulator == "iss" and args.read_length != 150:
        print(f"  note: --read-length is ignored with --simulator iss; "
              f"model {args.iss_model} uses {read_len} bp", file=sys.stderr)
    if read_len > genome_len:
        die(f"read length {read_len} exceeds genome length {genome_len}")

    os.makedirs(args.outdir, exist_ok=True)
    genome_dir = os.path.join(args.outdir, "genomes")
    fastq_dir = os.path.join(args.outdir, "fastq")
    vcf_dir = os.path.join(args.outdir, "vcf")
    for d in (genome_dir, fastq_dir, vcf_dir):
        os.makedirs(d, exist_ok=True)

    rng = random.Random(args.seed)
    snp_counts = [rng.randint(args.min_snps, args.max_snps) for _ in range(args.variants)]
    # Per-variant seeds derived from the master seed keep each variant
    # independently reproducible.
    seeds = [rng.randrange(1, 2**31 - 1) for _ in range(args.variants)]

    print(f"assembly    {args.assembly}")
    print(f"            {len(recs)} contig(s), {genome_len:,} bp total")
    if len(recs) > 1:
        print(f"            ! snp-mutator will concatenate these into one sequence")
    print(f"variants    {args.variants}, SNPs drawn from "
          f"{args.min_snps}..{args.max_snps} (seed {args.seed})")
    print(f"reads       {args.simulator}, {args.depth:g}x, {read_len} bp paired")
    print()

    width = max(2, len(str(args.variants)))
    manifest, truth = [], []
    ref_out = os.path.join(args.outdir, "reference.fasta")

    for idx, (n_snps, vseed) in enumerate(zip(snp_counts, seeds), start=1):
        sample = f"{args.prefix}{idx:0{width}d}"
        work = os.path.join(genome_dir, f".work_{sample}")
        os.makedirs(work, exist_ok=True)
        vcf = os.path.join(vcf_dir, f"{sample}.vcf")
        summary = os.path.join(vcf_dir, f"{sample}_summary.tsv")

        # ---- 1. introduce SNPs ------------------------------------------- #
        # One snpmutator call per variant: -s is a fixed count per invocation,
        # so a per-variant call is what lets each variant take a different
        # number of SNPs from the requested range.
        cmd = ["snpmutator", args.assembly,
               "-n", "1",
               "-s", str(n_snps),
               "-i", str(args.insertions),
               "-d", str(args.deletions),
               "-r", str(vseed),
               "-I", seqid,
               "-F", work,
               "-o", summary,
               "-v", vcf]
        if idx == 1:
            cmd += ["-R", ref_out]
        run(cmd, "snpmutator")

        produced = sorted(glob.glob(os.path.join(work, "*_mutated_*.fasta")))
        if not produced:
            die(f"snpmutator produced no FASTA in {work}")
        genome = os.path.join(genome_dir, f"{sample}.fasta")
        shutil.move(produced[0], genome)
        shutil.rmtree(work, ignore_errors=True)

        for chrom, pos, ref, alt in parse_vcf(vcf):
            truth.append({"sample": sample, "chrom": chrom, "pos": pos,
                          "ref": ref, "alt": alt})

        # ---- 2. simulate reads ------------------------------------------- #
        prefix = os.path.join(fastq_dir, sample)
        if args.simulator == "wgsim":
            # wgsim -N counts read PAIRS; each pair contributes 2 x read_len bases.
            pairs = max(1, round(args.depth * genome_len / (2 * read_len)))
            r1, r2 = prefix + "_R1.fastq", prefix + "_R2.fastq"
            run(["wgsim",
                 "-N", str(pairs),
                 "-1", str(read_len), "-2", str(read_len),
                 "-r", "0", "-R", "0",          # no wgsim mutations; SNPs come from snpmutator
                 "-e", str(args.error_rate),
                 "-d", str(args.fragment_length), "-s", str(args.fragment_sd),
                 "-S", str(vseed),
                 genome, r1, r2], "wgsim")
        else:
            # iss --n_reads counts reads across both mates.
            n_reads = max(2, round(args.depth * genome_len / read_len))
            cmd = ["iss", "generate",
                   "--genomes", genome,
                   "--n_reads", str(n_reads),
                   "--cpus", str(args.threads),
                   "--seed", str(vseed),
                   "--fragment-length", str(args.fragment_length),
                   "--fragment-length-sd", str(args.fragment_sd),
                   "--output", prefix,
                   "--quiet",
                   "--mode", args.iss_mode]
            if args.iss_mode == "kde":
                cmd += ["--model", args.iss_model]
            run(cmd, "iss")
            r1, r2 = prefix + "_R1.fastq", prefix + "_R2.fastq"
            for stray in glob.glob(prefix + ".iss.tmp*") + glob.glob(prefix + "_abundance.txt"):
                os.remove(stray)
        for f in (r1, r2):
            if not os.path.isfile(f):
                die(f"{args.simulator} did not write {f}")

        # ---- 3. verify realized depth ------------------------------------ #
        obs_bases = fastq_bases(r1) + fastq_bases(r2)
        obs_depth = obs_bases / genome_len

        # Compress only after read simulation: iss reads FASTA through
        # Biopython, which will not take a gzipped file.
        if not args.no_gzip:
            r1, r2 = compress(r1), compress(r2)
            genome = compress(genome)

        manifest.append({
            "sample": sample, "n_snps": n_snps, "seed": vseed,
            "genome_length": genome_len,
            "target_depth": f"{args.depth:g}",
            "realized_depth": f"{obs_depth:.2f}",
            "total_bases": obs_bases, "read_length": read_len,
            "simulator": args.simulator,
            "fastq_1": os.path.abspath(r1), "fastq_2": os.path.abspath(r2),
            "genome_fasta": os.path.abspath(genome),
            "vcf": os.path.abspath(vcf),
        })
        print(f"  {sample}  {n_snps:>5} SNPs   {obs_depth:6.2f}x   "
              f"{os.path.basename(r1)}, {os.path.basename(r2)}")

    # ---- write metadata -------------------------------------------------- #
    man_path = os.path.join(args.outdir, "manifest.tsv")
    man_cols = ["sample", "n_snps", "seed", "genome_length", "target_depth",
                "realized_depth", "total_bases", "read_length", "simulator",
                "fastq_1", "fastq_2", "genome_fasta", "vcf"]
    with open(man_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=man_cols, delimiter="\t")
        w.writeheader()
        w.writerows(manifest)

    truth_path = os.path.join(args.outdir, "truth_snps.tsv")
    with open(truth_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "chrom", "pos", "ref", "alt"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(truth)

    depths = [float(m["realized_depth"]) for m in manifest]
    print(f"\nSNPs per variant   min {min(snp_counts)}, max {max(snp_counts)}, "
          f"total {sum(snp_counts)} introduced, {len(truth)} rows in truth set")
    print(f"realized depth     {min(depths):.2f}x - {max(depths):.2f}x "
          f"(target {args.depth:g}x)")
    print("\nwritten")
    print(f"  {man_path}")
    print(f"  {truth_path}")
    print(f"  {ref_out}")
    print(f"  {vcf_dir}/  ({args.variants} VCF + summary)")
    print(f"  {fastq_dir}/  ({2 * args.variants} FASTQ)")
    print(f"  {genome_dir}/  ({args.variants} mutated assembly)")


if __name__ == "__main__":
    main()