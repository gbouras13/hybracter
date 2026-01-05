"""
Script to combine hybracter's flye and plassembler GFAs into one output

Some may find this useful so I have included it, but I (George) do not have the development time to continue development/improvements/fix bugs, so caveat emptor when using this

Contributed by bananabenana https://github.com/gbouras13/hybracter/issues/123 with some minor modifications - thanks so much Ben

Requires minimap2

"""

#!/usr/bin/env python3
import argparse
import pathlib
import subprocess
import tempfile
import sys
from collections import defaultdict
import shutil

def parse_gfa_raw(gfa_path):
    """Parse GFA files"""
    segments = {}
    lines = []

    with open(gfa_path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                parts = line.rstrip().split("\t")
                segments[parts[1]] = parts[2]
            lines.append(line.rstrip())

    return segments, lines

def parse_fasta(fasta_path):
    """Parse FASTA and use only the first token of the header as the sequence ID."""
    seqs = {}

    with open(fasta_path) as fh:
        sid, seq = None, []
        for line in fh:
            if line.startswith(">"):
                if sid:
                    seqs[sid] = "".join(seq)
                sid = line[1:].strip().split()[0]  # take only the first token - stupid fasta names
                seq = []
            else:
                seq.append(line.strip())
        if sid:
            seqs[sid] = "".join(seq)
    return seqs


def run_minimap(query_fa, ref_fa, paf_out, threads):
    """Run minimap"""
    cmd = ["minimap2", "-x", "asm10", "-t", str(threads), ref_fa, query_fa]
    with open(paf_out, "w") as out:
        subprocess.check_call(cmd, stdout=out)


def parse_paf(paf_path, segment_origin, min_mapq, min_cov):
    """
    Parses the paf file from minimap2. Return a dict: {GFA_segment_id: final_fasta_id} for segments to keep.
    Prioritize plassembler hits over assembly hits if multiple segments map to same final contig.
    """
    hits_by_ref = defaultdict(list)

    with open(paf_path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            rname = parts[5]
            matches = int(parts[9])
            mapq = int(parts[11])

            if mapq < min_mapq:
                continue
            cov = (qend - qstart) / qlen
            if cov < min_cov:
                continue

            hits_by_ref[rname].append({
                "query": qname,
                "matches": matches,
                "origin": segment_origin[qname]  # assembly | plassembler
            })

    kept = {}
    for ref, hits in hits_by_ref.items():
        # Separate plassembler vs assembly
        pls_hits = [h for h in hits if h["origin"] == "plassembler"]
        asm_hits = [h for h in hits if h["origin"] == "assembly"]

        if pls_hits:
            winner = max(pls_hits, key=lambda h: h["matches"])
        else:
            winner = max(asm_hits, key=lambda h: h["matches"])
        kept[winner["query"]] = ref  # map GFA segment -> final.fasta contig

    return kept


def rewrite_gfa_filtered(lines, final_seqs, kept):
    """
    Replace sequences for kept segments and remove segments/links/paths not in kept.
    kept: dict {GFA_segment_id -> final_fasta_contig_id}
    """
    out = []

    # Start with segments that mapped to final.fasta
    kept_set = set(kept.keys())

    for line in lines:
        if line.startswith("S\t"):
            parts = line.rstrip().split("\t")
            sid = parts[1]
            if sid not in kept_set:
                continue  # drop segments not kept
            final_id = kept[sid]
            seq = final_seqs.get(final_id)
            if seq is None:
                print(f"[WARN] No sequence found in final.fasta for {sid} -> {final_id}", file=sys.stderr)
                kept_set.remove(sid)  # remove from kept_set to prevent dangling links
                continue
            parts[2] = seq
            out.append("\t".join(parts))

        elif line.startswith("L\t"):
            parts = line.rstrip().split("\t")
            # only keep links if both segments exist in kept_set
            if parts[1] not in kept_set or parts[3] not in kept_set:
                continue
            out.append("\t".join(parts))

        elif line.startswith("P\t"):
            parts = line.rstrip().split("\t")
            new_path = []
            for p in parts[2].split(","):
                sid, orient = p[:-1], p[-1]
                if sid in kept_set:
                    new_path.append(sid + orient)
            if not new_path:
                continue
            parts[2] = ",".join(new_path)
            out.append("\t".join(parts))

        else:
            out.append(line)

    return out


def validate_gfa(lines):
    """Make sure there are no dangling links in the newly created gfa file - needs to match the input fasta basically"""
    segs = {l.split("\t")[1] for l in lines if l.startswith("S\t")}
    for l in lines:
        if l.startswith("L\t"):
            p = l.split("\t")
            if p[1] not in segs or p[3] not in segs:
                raise RuntimeError(f"Dangling link: {l}")


def process_genome(genome, indir, outdir, args):
    """Function which processes the actual pipeline"""
    asm_gfa = indir / "processing" / "assemblies" / genome / "assembly_graph.gfa"
    pls_gfa = indir / "processing" / "plassembler" / genome / "plassembler_plasmids.gfa"

    # Check assembly GFA
    if not asm_gfa.exists():
        print(f"[WARN] Missing assembly GFA for {genome}.", file=sys.stderr)
        return

    # Plassembler GFA optional
    if not pls_gfa.exists():
        pls_segs, pls_lines = {}, []
    else:
        pls_segs, pls_lines = parse_gfa_raw(pls_gfa)

    # Locate final.fasta (complete or incomplete)
    final_fa = None
    for status in ("complete", "incomplete"):
        p = indir / "FINAL_OUTPUT" / status / f"{genome}_final.fasta"
        if p.exists():
            final_fa = p
            break

    if not final_fa:
        print(f"[WARN] Missing final.fasta for {genome}.", file=sys.stderr)
        return

    # Load final.fasta sequences
    final_seqs = parse_fasta(final_fa)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = pathlib.Path(tmpdir)

        # Parse assembly GFA only (plassembler already handled above)
        asm_segs, asm_lines = parse_gfa_raw(asm_gfa)

        # Combine segments for minimap
        all_segs = {**asm_segs, **pls_segs}
        seg_fa = tmp / "segments.fasta"
        with open(seg_fa, "w") as fh:
            for sid, seq in all_segs.items():
                fh.write(f">{sid}\n{seq}\n")


        # Map segments to final.fasta
        paf = tmp / "map.paf"
        run_minimap(query_fa=seg_fa, ref_fa=final_fa, paf_out=paf, threads=args.threads)

        # Annotate segment origin
        segment_origin = {sid: "assembly" for sid in asm_segs}
        segment_origin.update({sid: "plassembler" for sid in pls_segs})

        # Parse PAF to find hits
        hits_by_ref = defaultdict(list)
        with open(paf) as fh:
            for line in fh:
                parts = line.rstrip().split("\t")
                qname = parts[0]
                qlen = int(parts[1])
                qstart = int(parts[2])
                qend = int(parts[3])
                rname = parts[5]
                matches = int(parts[9])
                mapq = int(parts[11])
                if mapq < args.min_mapq:
                    continue
                cov = (qend - qstart) / qlen
                if cov < args.min_cov:
                    continue
                hits_by_ref[rname].append({
                    "query": qname,
                    "matches": matches,
                    "origin": segment_origin[qname]
                })

        # Select best hit per final contig, prioritizing plassembler
        kept = {}
        for ref, hits in hits_by_ref.items():
            pls_hits = [h for h in hits if h["origin"] == "plassembler"]
            asm_hits = [h for h in hits if h["origin"] == "assembly"]
            if pls_hits:
                winner = max(pls_hits, key=lambda h: h["matches"])
            else:
                winner = max(asm_hits, key=lambda h: h["matches"])
            kept[winner["query"]] = ref

        # Combine all GFA lines
        combined_lines = asm_lines + pls_lines

        # Rewrite GFA: replace sequences with final.fasta, keep only segments in kept
        final_lines = []
        kept_set = set(kept.keys())
        for line in combined_lines:
            if line.startswith("S\t"):
                parts = line.rstrip().split("\t")
                sid = parts[1]
                if sid not in kept_set:
                    continue
                final_id = kept[sid]
                seq = final_seqs.get(final_id)
                if seq is None:
                    print(f"[WARN] No sequence in final.fasta for {sid} -> {final_id}", file=sys.stderr)
                    kept_set.remove(sid)
                    continue
                parts[2] = seq
                final_lines.append("\t".join(parts))
            elif line.startswith("L\t"):
                parts = line.rstrip().split("\t")
                if parts[1] not in kept_set or parts[3] not in kept_set:
                    continue
                final_lines.append("\t".join(parts))
            elif line.startswith("P\t"):
                parts = line.rstrip().split("\t")
                new_path = []
                for p in parts[2].split(","):
                    sid, orient = p[:-1], p[-1]
                    if sid in kept_set:
                        new_path.append(sid + orient)
                if not new_path:
                    continue
                parts[2] = ",".join(new_path)
                final_lines.append("\t".join(parts))
            else:
                final_lines.append(line)

        # Validate links
        segs_in_final = {l.split("\t")[1] for l in final_lines if l.startswith("S\t")}
        filtered_lines = []
        for l in final_lines:
            if l.startswith("L\t"):
                parts = l.split("\t")
                if parts[1] not in segs_in_final or parts[3] not in segs_in_final:
                    continue
            filtered_lines.append(l)

        # Write final GFA
        out_gfa = outdir / f"{genome}_final.gfa"
        with open(out_gfa, "w") as out:
            for l in filtered_lines:
                out.write(l + "\n")

    print(f"[OK] {genome}")


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--indir", required=True, type=pathlib.Path)
    ap.add_argument("-o","--outdir", required=True, type=pathlib.Path)
    ap.add_argument("-t","--threads", type=int, default=4)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-cov", type=float, default=0.5)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    return args


def main():
    args = parse_args()

    # checkf or minimap2

    mm2 = shutil.which("minimap2")
    if mm2 is None:
        sys.exit("ERROR: minimap2 is not installed or not found in $PATH. Please install it to use concat-gfas")

    try:
        result = subprocess.run(
            ["minimap2", "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        version = result.stdout.strip()
        print(f"minimap2 version {version} is installed")
    except Exception as e:
        sys.exit(f"ERROR: minimap2 found at {mm2} but failed to get version: {e}")


    # Get list of genomes and their relevant files
    genomes = sorted(
        p.name for p in (args.indir / "processing" / "assemblies").iterdir()
        if p.is_dir()
    )

    # Run each genome in the genome list through the process_genome function
    for genome in genomes:
        process_genome(genome, args.indir, args.outdir, args)

if __name__ == "__main__":
    main()
