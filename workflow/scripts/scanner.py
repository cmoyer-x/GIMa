#!/usr/bin/env python3
"""
mabs_island_scanner.py  —  Mycobacterium abscessus Genomic Island Scanner

Scans bacterial genomes for genomic islands using five independent
compositional and structural evidence signals, without requiring
comparative genomics (PPanGGOLiN) or external databases.

Evidence signals:
  1. GC content deviation   — foreign GC composition (z-score vs genome mean)
  2. CAI deviation          — codon usage divergence from host (age estimate)
  3. Mobility genes         — integrases, transposases, recombinases
  4. tRNA proximity         — integration at tRNA attB sites
  5. Direct repeats         — integration scars at island boundaries

Each candidate island is scored 0–1 across all five signals and assigned
a confidence tier: high (≥4 signals), moderate (3 signals), low (1-2 signals).

Validated on M. abscessus subsp. abscessus clinical isolates. Designed to
find both variable islands (also detectable by PPanGGOLiN) and fixed islands
present in all strains of a population.

Usage:
    # Single genome
    python mabs_island_scanner.py \\
        --fasta GD538.fasta \\
        --gff   GD538_prokka.gff \\
        --out   GD538_islands.tsv

    # Batch mode (all genomes in directories)
    python mabs_island_scanner.py \\
        --fasta_dir ATCC19977.FASTA \\
        --gff_dir   all_gffs \\
        --out_dir   island_predictions \\
        --threads   4

    # With existing RGP validation set
    python mabs_island_scanner.py \\
        --fasta_dir ATCC19977.FASTA \\
        --gff_dir   all_gffs \\
        --out_dir   island_predictions \\
        --validate  defense_rgp_intersection.tsv

Authors: Casey Moyer, University of Pittsburgh
"""

import os, re, csv, math, json, argparse, sys
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed


# ── Constants ─────────────────────────────────────────────────────────────────
VERSION = "0.1.0"

# Sliding window parameters
WINDOW_BP    = 8000    # window size in bp
STEP_BP      = 2000    # step size (75% overlap)
MIN_GENES    = 4       # minimum CDS in a candidate window
MERGE_GAP    = 10000   # merge adjacent candidate windows within this gap
MIN_ISLAND   = 4000    # minimum island length after merging
MAX_ISLAND   = 400000  # maximum island length (avoid calling whole chromosomes)

# Scoring thresholds
GC_Z_THRESH   = 2.0    # |z-score| >= this for GC signal (used only when std_gc reliable)
GC_ABS_THRESH = 0.01   # absolute GC deviation >= 1% from genome mean triggers GC signal
                       # More reliable than z-score for complete/near-complete genomes
CAI_THRESH    = 0.92   # CAI ratio <= this for CAI signal
TRNA_WINDOW   = 15000  # bp search window for tRNA proximity (15kb covers most attB sites)
DR_FLANK      = 150    # bp flank for direct repeat search
DR_MIN        = 8      # minimum direct repeat length
DR_MAX        = 25     # maximum direct repeat length

# Evidence weights for composite score (tuned on M. abscessus dataset)
WEIGHTS = {
    "gc":       0.25,
    "cai":      0.20,
    "mobility": 0.25,
    "trna":     0.15,
    "dr":       0.15,
}


# ── Codon table ───────────────────────────────────────────────────────────────
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
SYN_GROUPS = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    if aa != '*': SYN_GROUPS[aa].append(codon)


# ── Mobility gene patterns ────────────────────────────────────────────────────
MOBILITY_PATTERNS = {
    "integrase":   [r"integrase", r"tyrosine.recombinase", r"serine.recombinase",
                    r"site.specific.recombinase", r"phage.int"],
    "transposase": [r"transposase", r"\btnp[abc]?\b", r"insertion.sequence",
                    r"\bis\d+", r"transpos"],
    "recombinase": [r"recombinase", r"resolvase", r"invertase", r"xer[cd]"],
    "rdf":         [r"recombination.directional", r"\brdf\b", r"\bxis\b", r"excisionase"],
    "relaxase":    [r"relaxase", r"mob[abc]\b", r"mobilization.protein"],
    "conjugation": [r"conjugal", r"type.iv.secretion", r"virb\d", r"tra[bcde]\b"],
}
MOBILITY_COMPILED = {k: [re.compile(p, re.IGNORECASE) for p in v]
                     for k, v in MOBILITY_PATTERNS.items()}


def classify_mobility(product, gene=""):
    text = f"{product} {gene}".lower()
    for mtype, patterns in MOBILITY_COMPILED.items():
        for pat in patterns:
            if pat.search(text):
                return mtype
    return None


def is_ribosomal(product):
    prod = product.lower()
    if "methyltransferase" in prod or "pseudouridine" in prod or "biogenesis" in prod:
        return False
    return any(k in prod for k in
               ["ribosomal protein", "ribosomal subunit protein",
                "30s ribosomal", "50s ribosomal",
                "small subunit protein", "large subunit protein"])


# ── FASTA loading ─────────────────────────────────────────────────────────────
def load_fasta(path):
    seqs, cid, cseq = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cid: seqs[cid] = "".join(cseq).upper()
                cid = line[1:].split()[0]; cseq = []
            else:
                cseq.append(line)
    if cid: seqs[cid] = "".join(cseq).upper()
    return seqs


# ── GFF parsing ───────────────────────────────────────────────────────────────
def load_gff(path):
    """Parse Prokka GFF3. Returns (cds_list, trna_list)."""
    cds_list, trna_list = [], []
    with open(path) as f:
        for line in f:
            if "##FASTA" in line: break
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            feature = parts[2]
            if feature not in ("CDS","tRNA","tmRNA","rRNA"): continue
            try:
                start = int(parts[3]) - 1
                end   = int(parts[4])
            except: continue
            contig = parts[0].replace("gnl|XXX|","")
            attrs  = {}
            for a in parts[8].split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k.strip()] = v.strip()
            rec = {
                "contig":    contig,
                "start":     start,
                "end":       end,
                "strand":    parts[6],
                "product":   attrs.get("product","").replace("%2C",","),
                "gene":      attrs.get("gene",""),
                "locus_tag": attrs.get("locus_tag",""),
            }
            if feature == "CDS":
                cds_list.append(rec)
            else:
                trna_list.append(rec)
    return cds_list, trna_list


# ── Sequence helpers ──────────────────────────────────────────────────────────
def gc_content(seq):
    if not seq: return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join(comp.get(b,'N') for b in reversed(seq))


def get_seq(seqs, contig, start, end, strand="+"):
    seq = seqs.get(contig)
    if not seq:
        base = re.sub(r'_\d+$', '', contig)
        for k in seqs:
            if re.sub(r'_\d+$', '', k) == base:
                seq = seqs[k]; break
    if not seq and len(seqs) == 1:
        seq = list(seqs.values())[0]
    if not seq: return None
    s = seq[start:end]
    return reverse_complement(s) if strand == "-" else s


# ── CAI calculation ───────────────────────────────────────────────────────────
def count_codons(cds_seq):
    counts = Counter()
    for i in range(0, len(cds_seq)-2, 3):
        codon = cds_seq[i:i+3]
        if len(codon)==3 and CODON_TABLE.get(codon,"") not in ("*",""):
            counts[codon] += 1
    return counts


def build_w_table(ref_codon_counts):
    total = Counter()
    for c in ref_codon_counts: total.update(c)
    rscu = {}
    for aa, codons in SYN_GROUPS.items():
        aa_total = sum(total[c] for c in codons)
        n = len(codons)
        for c in codons:
            rscu[c] = (total[c]/aa_total)*n if aa_total > 0 else 1.0
    w = {}
    for aa, codons in SYN_GROUPS.items():
        mx = max(rscu.get(c,0) for c in codons) or 1.0
        for c in codons:
            w[c] = max(rscu.get(c,0)/mx, 0.001)
    return w


def calc_cai(cds_seq, w_table):
    codons = [cds_seq[i:i+3] for i in range(0,len(cds_seq)-2,3)
              if len(cds_seq[i:i+3])==3
              and CODON_TABLE.get(cds_seq[i:i+3],"") not in ("*","")]
    if not codons: return None
    return math.exp(sum(math.log(w_table.get(c,0.001)) for c in codons)/len(codons))


# ── Direct repeat detection ───────────────────────────────────────────────────
def is_low_complexity(seq):
    if len(set(seq)) <= 2: return True
    for nt in "ACGT":
        if seq.count(nt)/len(seq) > 0.70: return True
    return False


def find_direct_repeats(left_seq, right_seq, min_k=DR_MIN, max_k=DR_MAX):
    if not left_seq or not right_seq: return None
    for k in range(max_k, min_k-1, -1):
        left_kmers = {}
        for i in range(len(left_seq)-k+1):
            kmer = left_seq[i:i+k]
            if kmer not in left_kmers:
                left_kmers[kmer] = i
        for j in range(len(right_seq)-k+1):
            kmer = right_seq[j:j+k]
            if kmer in left_kmers and not is_low_complexity(kmer):
                return {"repeat_seq": kmer, "repeat_len": k,
                        "left_pos": left_kmers[kmer], "right_pos": j}
    return None


# ── RGP seed loading ─────────────────────────────────────────────────────────
def load_rgp_seeds(rgp_file):
    """
    Load RGP coordinates from PPanGGOLiN output as seed regions.
    Returns dict: strain -> list of {contig, start, stop, rgp_id}
    """
    seeds = defaultdict(list)
    with open(rgp_file) as f:
        raw_header = f.readline().strip().split("\t")
    col = {h: i for i, h in enumerate(raw_header)}
    with open(rgp_file) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 8: continue
            try:
                strain = row[col["genome"]]
                contig = row[col["contig"]].replace("gnl|XXX|","")
                start  = int(row[col["start"]])
                stop   = int(row[7])
                rgp_id = row[col["region"]]
                seeds[strain].append({
                    "contig": contig,
                    "start":  start,
                    "stop":   stop,
                    "rgp_id": rgp_id,
                })
            except: continue
    total = sum(len(v) for v in seeds.values())
    print(f"  Loaded {total} RGP seeds across {len(seeds)} strains")
    return seeds


def analyze_rgp_seed(strain, seed, seqs, gene_scores, trna_list,
                     genome_mean_cai, genome_mean_gc, genome_std_gc,
                     w_table, extend_bp=20000):
    """
    Score a single RGP seed region by extending its boundaries and
    computing all five evidence signals across the full region.
    Returns an island dict or None.
    """
    contig    = seed["contig"]
    rgp_id    = seed["rgp_id"]

    # Get genome sequence
    seq = None
    for k, s in seqs.items():
        if re.sub(r'_\d+$','',k) == re.sub(r'_\d+$','',contig) or k == contig:
            seq = s; break
    if not seq and len(seqs) == 1:
        seq = list(seqs.values())[0]
    if not seq: return None

    # Start with RGP boundaries then extend
    isl_start = max(0, seed["start"] - extend_bp)
    isl_end   = min(len(seq), seed["stop"] + extend_bp)

    # Trim to gene boundaries — find outermost genes with signal
    contig_genes = [g for g in gene_scores
                    if re.sub(r'_\d+$','',g["contig"])==re.sub(r'_\d+$','',contig)
                    and g["start"] >= isl_start and g["end"] <= isl_end]

    if not contig_genes:
        # Fall back to raw RGP coordinates
        isl_start = seed["start"]
        isl_end   = seed["stop"]
        contig_genes = [g for g in gene_scores
                        if re.sub(r'_\d+$','',g["contig"])==re.sub(r'_\d+$','',contig)
                        and g["start"] >= isl_start and g["end"] <= isl_end]

    # Trim boundaries to first/last gene
    if contig_genes:
        isl_start = min(g["start"] for g in contig_genes)
        isl_end   = max(g["end"]   for g in contig_genes)

    isl_len = isl_end - isl_start
    if isl_len < MIN_ISLAND or isl_len > MAX_ISLAND:
        # Use exact RGP coords if extension gives bad bounds
        isl_start = seed["start"]
        isl_end   = seed["stop"]
        isl_len   = isl_end - isl_start
        contig_genes = [g for g in gene_scores
                        if re.sub(r'_\d+$','',g["contig"])==re.sub(r'_\d+$','',contig)
                        and g["start"] >= isl_start and g["end"] <= isl_end]

    # ── Evidence signals ──────────────────────────────────────────────────────

    # GC content signal
    isl_seq = seq[isl_start:isl_end]
    isl_gc  = gc_content(isl_seq)
    gc_z     = (isl_gc - genome_mean_gc) / genome_std_gc if genome_std_gc > 0 else 0
    gc_dev   = abs(isl_gc - genome_mean_gc)
    gc_score = min(1.0, gc_dev / 0.03)  # 3% deviation = full score

    # CAI signal
    gene_cais = [g["cai"] for g in contig_genes if g["cai"] is not None]
    if gene_cais:
        mean_cai   = sum(gene_cais) / len(gene_cais)
        cai_ratio  = mean_cai / genome_mean_cai
        cai_score  = max(0, min(1.0, (CAI_THRESH - cai_ratio) / 0.15))
    else:
        mean_cai  = genome_mean_cai
        cai_ratio = 1.0
        cai_score = 0.0

    # Mobility signal
    mob_genes  = [g for g in contig_genes if g["is_mob"]]
    mob_types  = list(set(g["mob_type"] for g in mob_genes if g["mob_type"]))
    mob_score  = min(1.0, len(mob_genes) / 2.0)

    # tRNA proximity signal
    trna_flanked  = False
    trna_product  = ""
    trna_distance = None
    for trna in trna_list:
        if re.sub(r'_\d+$','',trna["contig"]) != re.sub(r'_\d+$','',contig):
            continue
        trna_mid = (trna["start"] + trna["end"]) // 2
        dist_l   = abs(trna_mid - isl_start)
        dist_r   = abs(trna_mid - isl_end)
        min_dist = min(dist_l, dist_r)
        if min_dist <= TRNA_WINDOW:
            if trna_distance is None or min_dist < trna_distance:
                trna_flanked  = True
                trna_product  = trna["product"]
                trna_distance = min_dist
    trna_score = 1.0 if trna_flanked else 0.0

    # Direct repeat signal
    left_seq  = seq[max(0,isl_start-DR_FLANK):isl_start+DR_FLANK]
    right_seq = seq[max(0,isl_end-DR_FLANK):isl_end+DR_FLANK]
    dr_result = find_direct_repeats(left_seq, right_seq)
    dr_score  = 1.0 if dr_result else 0.0

    # Composite score
    island_score = (WEIGHTS["gc"]       * gc_score +
                    WEIGHTS["cai"]      * cai_score +
                    WEIGHTS["mobility"] * mob_score +
                    WEIGHTS["trna"]     * trna_score +
                    WEIGHTS["dr"]       * dr_score)

    # Evidence lines
    # GC: use absolute deviation (more reliable than z-score for complete genomes)
    gc_deviation = abs(isl_gc - genome_mean_gc)
    evidence = []
    if gc_deviation >= GC_ABS_THRESH:     evidence.append("gc_foreign")
    if cai_ratio  <= CAI_THRESH:          evidence.append("cai_deviation")
    if mob_genes:                         evidence.append("mobility_gene")
    if trna_flanked:                      evidence.append("trna_flanked")
    if dr_result:                         evidence.append("direct_repeat")

    n_evidence = len(evidence)
    if n_evidence >= 4:     confidence = "high"
    elif n_evidence >= 3:   confidence = "moderate"
    elif n_evidence >= 2:   confidence = "low"
    else:                   confidence = "very_low"

    if cai_ratio < 0.85:      age = "very_recent"
    elif cai_ratio < 0.92:    age = "recent"
    elif cai_ratio < 0.97:    age = "moderate"
    else:                     age = "old"

    n_genes  = len(contig_genes)
    n_mob    = len(mob_genes)
    n_hypo   = sum(1 for g in contig_genes if "hypothetical" in g["product"].lower())
    pct_hypo = round(100*n_hypo/n_genes, 1) if n_genes > 0 else 0

    return {
        "strain":        strain,
        "contig":        contig,
        "start":         isl_start,
        "end":           isl_end,
        "length":        isl_len,
        "rgp_seed":      rgp_id,
        "n_genes":       n_genes,
        "n_mob_genes":   n_mob,
        "mob_types":     "; ".join(mob_types),
        "pct_hypo":      pct_hypo,
        "isl_gc":        round(isl_gc*100, 2),
        "gc_z":          round(gc_z, 3),
        "mean_cai":      round(mean_cai, 4),
        "cai_ratio":     round(cai_ratio, 4),
        "age_estimate":  age,
        "trna_flanked":  "Yes" if trna_flanked else "No",
        "trna_product":  trna_product,
        "trna_distance": trna_distance if trna_distance else "",
        "has_dr":        "Yes" if dr_result else "No",
        "dr_seq":        dr_result["repeat_seq"] if dr_result else "",
        "dr_len":        dr_result["repeat_len"] if dr_result else "",
        "n_evidence":    n_evidence,
        "evidence":      "; ".join(evidence),
        "island_score":  round(island_score, 4),
        "confidence":    confidence,
    }


# ── Per-genome analysis ───────────────────────────────────────────────────────
def analyze_genome(strain, seqs, cds_list, trna_list):
    """
    Full genomic island scan for one genome.
    Returns list of island candidate dicts.
    """

    # ── Step 1: Build CAI reference table ────────────────────────────────────
    ref_counts = []
    for cds in cds_list:
        if not is_ribosomal(cds["product"]): continue
        seq = get_seq(seqs, cds["contig"], cds["start"], cds["end"], cds["strand"])
        if seq and len(seq) >= 60:
            counts = count_codons(seq)
            if sum(counts.values()) >= 20:
                ref_counts.append(counts)

    if len(ref_counts) < 5:
        return []  # Not enough reference genes

    w_table = build_w_table(ref_counts)

    # ── Step 2: Per-gene scoring ──────────────────────────────────────────────
    gene_scores = []
    all_cai     = []

    for cds in cds_list:
        seq = get_seq(seqs, cds["contig"], cds["start"], cds["end"], cds["strand"])
        cai = None
        if seq and len(seq) >= 60:
            cai = calc_cai(seq, w_table)
            if cai: all_cai.append(cai)

        mob_type = classify_mobility(cds["product"], cds["gene"])

        gene_scores.append({
            **cds,
            "cai":      cai,
            "mob_type": mob_type,
            "is_mob":   mob_type is not None,
        })

    if not all_cai: return []
    genome_mean_cai = sum(all_cai) / len(all_cai)

    # ── Step 3: Per-contig GC statistics ─────────────────────────────────────
    # Compute robust GC statistics using median + MAD
    # This prevents islands from inflating the reference std_gc
    all_window_gcs_flat = []
    contig_gc_stats = {}
    for contig_name, seq in seqs.items():
        if len(seq) < 1000: continue
        window_gcs = []
        for i in range(0, len(seq)-WINDOW_BP, STEP_BP):
            window_gcs.append(gc_content(seq[i:i+WINDOW_BP]))
        if len(window_gcs) < 10: continue
        all_window_gcs_flat.extend(window_gcs)
        # Use median for robustness
        sorted_gcs = sorted(window_gcs)
        median_gc  = sorted_gcs[len(sorted_gcs)//2]
        # MAD (median absolute deviation) — robust std estimator
        mad = sorted([abs(g - median_gc) for g in sorted_gcs])[len(sorted_gcs)//2]
        # Scale MAD to approximate std (factor 1.4826 for normal distribution)
        robust_std = max(mad * 1.4826, 0.002)
        contig_gc_stats[contig_name] = (median_gc, robust_std)

    # Genome-wide robust stats
    if all_window_gcs_flat:
        sorted_all = sorted(all_window_gcs_flat)
        genome_mean_gc = sorted_all[len(sorted_all)//2]  # median
        mad_all = sorted([abs(g-genome_mean_gc) for g in sorted_all])[len(sorted_all)//2]
        genome_std_gc  = max(mad_all * 1.4826, 0.002)
    else:
        genome_mean_gc = 0.64
        genome_std_gc  = 0.015

    # ── Step 4: Sliding window scan ───────────────────────────────────────────
    # Group genes by contig
    genes_by_contig = defaultdict(list)
    for g in gene_scores:
        genes_by_contig[g["contig"]].append(g)

    candidate_windows = []

    for contig_name, seq in seqs.items():
        genes = genes_by_contig.get(contig_name, [])
        # Also try Prokka-suffixed name
        if not genes:
            for k in genes_by_contig:
                if re.sub(r'_\d+$','',k) == re.sub(r'_\d+$','',contig_name):
                    genes = genes_by_contig[k]; break
        if not genes: continue

        mean_gc, std_gc = contig_gc_stats.get(contig_name,
                          (genome_mean_gc, genome_std_gc))

        for win_start in range(0, max(1, len(seq)-WINDOW_BP), STEP_BP):
            win_end   = min(win_start + WINDOW_BP, len(seq))
            win_seq   = seq[win_start:win_end]
            win_genes = [g for g in genes
                         if g["start"] >= win_start and g["end"] <= win_end]

            if len(win_genes) < MIN_GENES:
                continue

            # GC signal
            win_gc   = gc_content(win_seq)
            gc_z     = (win_gc - mean_gc) / std_gc if std_gc > 0 else 0
            gc_score = min(1.0, abs(win_gc - mean_gc) / 0.03)  # 3% abs deviation = full score

            # CAI signal
            win_cais = [g["cai"] for g in win_genes if g["cai"] is not None]
            if win_cais:
                win_mean_cai = sum(win_cais)/len(win_cais)
                cai_ratio    = win_mean_cai / genome_mean_cai
                cai_score    = max(0, min(1.0, (CAI_THRESH - cai_ratio) / 0.15))
            else:
                cai_ratio = 1.0; cai_score = 0.0

            # Mobility signal
            mob_genes  = [g for g in win_genes if g["is_mob"]]
            mob_score  = min(1.0, len(mob_genes) / 2.0)

            # Composite window score
            composite = (WEIGHTS["gc"]  * gc_score +
                         WEIGHTS["cai"] * cai_score +
                         WEIGHTS["mobility"] * mob_score)

            # Only keep windows with at least weak signal
            if composite < 0.10 and not mob_genes:
                continue

            candidate_windows.append({
                "contig":     contig_name,
                "start":      win_start,
                "end":        win_end,
                "win_gc":     win_gc,
                "gc_z":       gc_z,
                "gc_score":   gc_score,
                "cai_ratio":  cai_ratio,
                "cai_score":  cai_score,
                "n_mob":      len(mob_genes),
                "mob_score":  mob_score,
                "mob_types":  list(set(g["mob_type"] for g in mob_genes if g["mob_type"])),
                "n_genes":    len(win_genes),
                "composite":  composite,
            })

    if not candidate_windows:
        return []

    # ── Step 5: Merge adjacent candidate windows ──────────────────────────────
    # Sort by contig then start
    candidate_windows.sort(key=lambda w: (w["contig"], w["start"]))

    merged = []
    current = [candidate_windows[0]]

    for win in candidate_windows[1:]:
        last = current[-1]
        if (win["contig"] == last["contig"] and
                win["start"] - last["end"] <= MERGE_GAP):
            current.append(win)
        else:
            merged.append(current)
            current = [win]
    merged.append(current)

    # ── Step 6: Build island candidates from merged windows ───────────────────
    islands = []

    for window_group in merged:
        contig    = window_group[0]["contig"]
        isl_start = window_group[0]["start"]
        isl_end   = window_group[-1]["end"]

        # Extend boundaries to capture full flanking genes
        # Find all genes within MERGE_GAP of the window boundaries
        contig_genes = [g for g in gene_scores
                        if re.sub(r'_\d+$','',g["contig"])==re.sub(r'_\d+$','',contig)]
        flanking = [g for g in contig_genes
                    if g["end"] >= isl_start - MERGE_GAP//2
                    and g["start"] <= isl_end + MERGE_GAP//2]
        if flanking:
            isl_start = min(isl_start, min(g["start"] for g in flanking))
            isl_end   = max(isl_end,   max(g["end"]   for g in flanking))

        isl_len = isl_end - isl_start
        if isl_len < MIN_ISLAND or isl_len > MAX_ISLAND:
            continue

        # Aggregate scores across merged windows
        n_wins   = len(window_group)
        mean_gc_score  = sum(w["gc_score"]  for w in window_group) / n_wins
        mean_cai_score = sum(w["cai_score"] for w in window_group) / n_wins
        max_mob_score  = max(w["mob_score"] for w in window_group)
        max_gc_z       = max(abs(w["gc_z"]) for w in window_group)
        mob_types      = list(set(t for w in window_group for t in w["mob_types"]))
        total_mob      = sum(w["n_mob"] for w in window_group)

        mean_cai_ratio = sum(w["cai_ratio"] for w in window_group) / n_wins

        # Get island sequence for GC and direct repeat computation
        seq = None
        for k, s in seqs.items():
            if re.sub(r'_\d+$','',k) == re.sub(r'_\d+$','',contig) or k == contig:
                seq = s; break
        if not seq and len(seqs)==1:
            seq = list(seqs.values())[0]

        isl_gc  = gc_content(seq[isl_start:isl_end]) if seq else 0.64
        mean_gc = sum(w["win_gc"] for w in window_group) / n_wins
        gc_z    = (isl_gc - mean_gc) / genome_std_gc if genome_std_gc > 0 else 0

        # tRNA proximity signal
        trna_flanked  = False
        trna_product  = ""
        trna_distance = None

        for trna in trna_list:
            if re.sub(r'_\d+$','',trna["contig"]) != re.sub(r'_\d+$','',contig):
                continue
            trna_mid = (trna["start"] + trna["end"]) // 2
            dist_l   = abs(trna_mid - isl_start)
            dist_r   = abs(trna_mid - isl_end)
            min_dist = min(dist_l, dist_r)
            if min_dist <= TRNA_WINDOW:
                if trna_distance is None or min_dist < trna_distance:
                    trna_flanked  = True
                    trna_product  = trna["product"]
                    trna_distance = min_dist

        trna_score = 1.0 if trna_flanked else 0.0

        # Direct repeat signal
        dr_result = None
        if seq:
            left_seq  = seq[max(0,isl_start-DR_FLANK):isl_start+DR_FLANK]
            right_seq = seq[max(0,isl_end-DR_FLANK):isl_end+DR_FLANK]
            dr_result = find_direct_repeats(left_seq, right_seq)

        dr_score  = 1.0 if dr_result else 0.0

        # Composite island score
        island_score = (WEIGHTS["gc"]       * min(1.0, abs(gc_z)/5.0) +
                        WEIGHTS["cai"]      * mean_cai_score +
                        WEIGHTS["mobility"] * max_mob_score +
                        WEIGHTS["trna"]     * trna_score +
                        WEIGHTS["dr"]       * dr_score)

        # Evidence count (lines of independent evidence)
        # GC: use absolute deviation (more reliable than z-score for complete genomes)
        gc_deviation = abs(isl_gc - genome_mean_gc)
        evidence = []
        if gc_deviation >= GC_ABS_THRESH:     evidence.append("gc_foreign")
        if mean_cai_ratio <= CAI_THRESH:      evidence.append("cai_deviation")
        if total_mob > 0:                     evidence.append("mobility_gene")
        if trna_flanked:                      evidence.append("trna_flanked")
        if dr_result:                         evidence.append("direct_repeat")

        n_evidence = len(evidence)

        # Confidence tier
        if n_evidence >= 4:     confidence = "high"
        elif n_evidence >= 3:   confidence = "moderate"
        elif n_evidence >= 2:   confidence = "low"
        else:                   confidence = "very_low"

        # Age estimate from CAI
        if mean_cai_ratio < 0.85:      age = "very_recent"
        elif mean_cai_ratio < 0.92:    age = "recent"
        elif mean_cai_ratio < 0.97:    age = "moderate"
        else:                          age = "old"

        # Gene count in island
        island_genes = [g for g in gene_scores
                        if re.sub(r'_\d+$','',g["contig"])==re.sub(r'_\d+$','',contig)
                        and g["start"] >= isl_start and g["end"] <= isl_end]
        n_genes    = len(island_genes)
        n_mob      = sum(1 for g in island_genes if g["is_mob"])
        n_hypo     = sum(1 for g in island_genes if "hypothetical" in g["product"].lower())
        pct_hypo   = round(100*n_hypo/n_genes, 1) if n_genes > 0 else 0

        islands.append({
            "strain":        strain,
            "contig":        contig,
            "start":         isl_start,
            "end":           isl_end,
            "length":        isl_len,
            "n_genes":       n_genes,
            "n_mob_genes":   n_mob,
            "mob_types":     "; ".join(mob_types),
            "pct_hypo":      pct_hypo,
            "isl_gc":        round(isl_gc*100, 2),
            "gc_z":          round(gc_z, 3),
            "mean_cai":      round(sum(w["cai_ratio"] for w in window_group)/n_wins, 4),
            "cai_ratio":     round(mean_cai_ratio, 4),
            "age_estimate":  age,
            "trna_flanked":  "Yes" if trna_flanked else "No",
            "trna_product":  trna_product,
            "trna_distance": trna_distance if trna_distance else "",
            "has_dr":        "Yes" if dr_result else "No",
            "dr_seq":        dr_result["repeat_seq"] if dr_result else "",
            "dr_len":        dr_result["repeat_len"] if dr_result else "",
            "n_evidence":    n_evidence,
            "evidence":      "; ".join(evidence),
            "island_score":  round(island_score, 4),
            "confidence":    confidence,
        })

    return islands


# ── Validation against known defense RGPs ────────────────────────────────────
def validate_against_rgps(all_islands, validation_file):
    """
    Compare predicted islands to known defense RGPs.
    Reports sensitivity (recall) and specificity.
    """
    known = []
    with open(validation_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("in_rgp") == "Yes":
                try:
                    known.append({
                        "strain":  row["strain"],
                        "start":   int(row["rgp_start"]),
                        "end":     int(row["rgp_stop"]),
                        "system":  row["subtype"],
                    })
                except: continue

    # Group predictions by strain
    pred_by_strain = defaultdict(list)
    for isl in all_islands:
        pred_by_strain[isl["strain"]].append(isl)

    tp, fn = 0, 0
    for k in known:
        strain_preds = pred_by_strain.get(k["strain"], [])
        # Check if any predicted island overlaps this known RGP
        hit = any(
            p["start"] <= k["end"] and p["end"] >= k["start"]
            for p in strain_preds
        )
        if hit: tp += 1
        else:   fn += 1

    sensitivity = tp / (tp + fn) if (tp+fn) > 0 else 0
    print(f"\n  Validation against known defense RGPs:")
    print(f"    Known defense RGPs:    {len(known)}")
    print(f"    Correctly predicted:   {tp}  (sensitivity {sensitivity:.1%})")
    print(f"    Missed:                {fn}")
    print(f"    Total islands called:  {len(all_islands)}")

    return {"tp": tp, "fn": fn, "sensitivity": sensitivity}


# ── Output ────────────────────────────────────────────────────────────────────
FIELDNAMES = [
    "strain","contig","start","end","length","rgp_seed","n_genes","n_mob_genes",
    "mob_types","pct_hypo","isl_gc","gc_z","mean_cai","cai_ratio",
    "age_estimate","trna_flanked","trna_product","trna_distance",
    "has_dr","dr_seq","dr_len","n_evidence","evidence","island_score","confidence",
]


def write_tsv(islands, outpath):
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(islands)


def print_summary(strain, islands):
    conf_counts = Counter(i["confidence"] for i in islands)
    high   = conf_counts.get("high",0)
    mod    = conf_counts.get("moderate",0)
    low    = conf_counts.get("low",0)+conf_counts.get("very_low",0)
    print(f"  {strain:30}: {len(islands):4} islands  "
          f"[high={high} mod={mod} low={low}]")


# ── Single-genome entry point ─────────────────────────────────────────────────
def run_single(fasta_path, gff_path, out_path, validate=None,
               rgp_seeds=None, extend_bp=20000):
    strain = os.path.basename(fasta_path).rsplit(".",1)[0]
    print(f"Scanning {strain}...")

    seqs            = load_fasta(fasta_path)
    cds_list, trnas = load_gff(gff_path)

    # Guided mode: score known RGP seeds + de novo scan
    if rgp_seeds and strain in rgp_seeds:
        islands = run_guided(strain, seqs, cds_list, trnas,
                             rgp_seeds[strain], extend_bp)
        # Also run de novo to catch islands not in RGP file
        de_novo = analyze_genome(strain, seqs, cds_list, trnas)
        # Add de novo islands that don't overlap with guided islands
        for isl in de_novo:
            overlaps = any(
                isl["start"] <= g["end"] and isl["end"] >= g["start"]
                and isl["contig"] == g["contig"]
                for g in islands
            )
            if not overlaps:
                isl["rgp_seed"] = ""
                islands.append(isl)
    else:
        islands = analyze_genome(strain, seqs, cds_list, trnas)
        for isl in islands:
            isl["rgp_seed"] = ""

    write_tsv(islands, out_path)
    print_summary(strain, islands)

    if validate and islands:
        validate_against_rgps(islands, validate)

    return islands


def run_guided(strain, seqs, cds_list, trnas, seeds, extend_bp):
    """Score all RGP seeds for this strain."""
    # Build gene scores first (needed by analyze_rgp_seed)
    # Reuse the gene scoring logic from analyze_genome

    # Build CAI reference
    ref_counts = []
    for cds in cds_list:
        if not is_ribosomal(cds["product"]): continue
        seq = get_seq(seqs, cds["contig"], cds["start"], cds["end"], cds["strand"])
        if seq and len(seq) >= 60:
            counts = count_codons(seq)
            if sum(counts.values()) >= 20:
                ref_counts.append(counts)

    if len(ref_counts) < 5: return []
    w_table = build_w_table(ref_counts)

    # Score all genes
    gene_scores = []
    all_cai = []
    for cds in cds_list:
        seq = get_seq(seqs, cds["contig"], cds["start"], cds["end"], cds["strand"])
        cai = None
        if seq and len(seq) >= 60:
            cai = calc_cai(seq, w_table)
            if cai: all_cai.append(cai)
        mob_type = classify_mobility(cds["product"], cds["gene"])
        gene_scores.append({**cds, "cai": cai, "mob_type": mob_type,
                            "is_mob": mob_type is not None})

    if not all_cai: return []
    genome_mean_cai = sum(all_cai) / len(all_cai)

    # Robust GC stats
    all_gcs = []
    for seq in seqs.values():
        if len(seq) < 1000: continue
        for i in range(0, len(seq)-WINDOW_BP, STEP_BP):
            all_gcs.append(gc_content(seq[i:i+WINDOW_BP]))

    if all_gcs:
        sorted_gcs    = sorted(all_gcs)
        genome_mean_gc = sorted_gcs[len(sorted_gcs)//2]
        mad = sorted([abs(g-genome_mean_gc) for g in sorted_gcs])[len(sorted_gcs)//2]
        genome_std_gc  = max(mad * 1.4826, 0.002)
    else:
        genome_mean_gc = 0.64; genome_std_gc = 0.015

    # Score each seed
    islands = []
    for seed in seeds:
        isl = analyze_rgp_seed(strain, seed, seqs, gene_scores, trnas,
                               genome_mean_cai, genome_mean_gc, genome_std_gc,
                               w_table, extend_bp)
        if isl: islands.append(isl)

    return islands


# ── Batch entry point ─────────────────────────────────────────────────────────
def run_batch(fasta_dir, gff_dir, out_dir, threads=1, validate=None, rgp_seeds=None, extend_bp=20000):
    os.makedirs(out_dir, exist_ok=True)

    # Build strain list
    tasks = []
    for fname in sorted(os.listdir(fasta_dir)):
        if not fname.endswith((".fasta",".fa",".fna")): continue
        strain    = fname.rsplit(".",1)[0]
        fasta_path = os.path.join(fasta_dir, fname)
        # GFF files may have _WGS/_hybrid stripped — try multiple name variants
        strain_base = re.sub(r'(_WGS|_hybrid|_UNCUT)$', '', strain, flags=re.IGNORECASE)
        gff_path = os.path.join(gff_dir, f"{strain}_prokka.gff")
        if not os.path.exists(gff_path):
            gff_path = os.path.join(gff_dir, f"{strain_base}_prokka.gff")
        if not os.path.exists(gff_path):
            gff_path = os.path.join(gff_dir, f"{strain}.gff")
        if not os.path.exists(gff_path):
            gff_path = os.path.join(gff_dir, f"{strain_base}.gff")
        if not os.path.exists(gff_path):
            continue
        out_path = os.path.join(out_dir, f"{strain}_islands.tsv")
        tasks.append((strain, fasta_path, gff_path, out_path))

    print(f"Scanning {len(tasks)} genomes (threads={threads})...")

    all_islands = []

    if threads > 1:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            futures = {
                ex.submit(run_single, t[1], t[2], t[3], None, rgp_seeds, extend_bp): t[0]
                for t in tasks
            }
            for future in as_completed(futures):
                strain = futures[future]
                try:
                    islands = future.result()
                    all_islands.extend(islands)
                    print_summary(strain, islands)
                except Exception as e:
                    print(f"  ERROR {strain}: {e}")
    else:
        for strain, fasta_path, gff_path, out_path in tasks:
            try:
                islands = run_single(fasta_path, gff_path, out_path,
                                     None, rgp_seeds, extend_bp)
                all_islands.extend(islands)
            except Exception as e:
                print(f"  ERROR {strain}: {e}")

    # Write combined output
    combined_path = os.path.join(out_dir, "all_islands_combined.tsv")
    write_tsv(all_islands, combined_path)
    print(f"\nCombined output: {combined_path}")

    # Summary stats
    conf = Counter(i["confidence"] for i in all_islands)
    print(f"\n{'='*55}")
    print(f"  Total islands predicted:  {len(all_islands)}")
    print(f"  High confidence:          {conf.get('high',0)}")
    print(f"  Moderate confidence:      {conf.get('moderate',0)}")
    print(f"  Low / very low:           {conf.get('low',0)+conf.get('very_low',0)}")
    print(f"  Strains processed:        {len(tasks)}")
    print(f"{'='*55}")

    if validate and all_islands:
        validate_against_rgps(all_islands, validate)

    return all_islands


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--fasta",     help="Single genome FASTA")
    mode.add_argument("--fasta_dir", help="Directory of FASTA files (batch mode)")

    p.add_argument("--gff",      help="Single genome GFF (required with --fasta)")
    p.add_argument("--gff_dir",  help="Directory of GFF files (required with --fasta_dir)")
    p.add_argument("--out",      help="Output TSV (single mode)",
                   default="islands.tsv")
    p.add_argument("--out_dir",  help="Output directory (batch mode)",
                   default="island_predictions")
    p.add_argument("--validate", help="Defense RGP intersection TSV for sensitivity test",
                   default=None)
    p.add_argument("--rgp_guided", help="RGP TSV from PPanGGOLiN to use as seed regions",
                   default=None)
    p.add_argument("--rgp_extend", type=int, default=20000,
                   help="bp to extend beyond RGP boundaries when guided (default 20000)")
    p.add_argument("--threads",  type=int, default=1,
                   help="Parallel threads for batch mode (default 1)")

    # Tunable parameters
    p.add_argument("--window",   type=int, default=WINDOW_BP,
                   help=f"Sliding window size bp (default {WINDOW_BP})")
    p.add_argument("--step",     type=int, default=STEP_BP,
                   help=f"Window step size bp (default {STEP_BP})")
    p.add_argument("--min_ev",   type=int, default=2,
                   help="Minimum evidence lines to report island (default 2)")

    return p.parse_args()


def main():
    args = parse_args()

    # Override globals from args
    global WINDOW_BP, STEP_BP
    WINDOW_BP = args.window
    STEP_BP   = args.step

    # Load RGP seeds if guided mode requested
    rgp_seeds = None
    if args.rgp_guided:
        print(f"  Loading RGP seeds from {args.rgp_guided}...")
        rgp_seeds = load_rgp_seeds(args.rgp_guided)

    if args.fasta:
        if not args.gff:
            print("ERROR: --gff required with --fasta", file=sys.stderr)
            sys.exit(1)
        islands = run_single(args.fasta, args.gff, args.out, args.validate,
                             rgp_seeds, args.rgp_extend)
        # Filter by min evidence
        filtered = [i for i in islands if i["n_evidence"] >= args.min_ev]
        if len(filtered) < len(islands):
            write_tsv(filtered, args.out)
            print(f"  Filtered to {len(filtered)} islands (n_evidence >= {args.min_ev})")
    else:
        if not args.gff_dir:
            print("ERROR: --gff_dir required with --fasta_dir", file=sys.stderr)
            sys.exit(1)
        run_batch(args.fasta_dir, args.gff_dir, args.out_dir,
                  args.threads, args.validate, rgp_seeds, args.rgp_extend)


if __name__ == "__main__":
    main()
