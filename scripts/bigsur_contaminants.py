"""Contaminant gene families to drop before BigSur feature selection.

The first cohort run showed feature selection being driven by carry-over rather
than tumour biology: haemoglobin topped the mcFano ranking in 25/39 samples,
crystallin in 21/39 and immunoglobulin in 20/39 (see results/bigsur/cohort_summary.csv).
Crystallins are lens structural proteins -- in enucleated-eye specimens they are
almost certainly lens carry-over, not expression (CRYGC~CRYGD correlate at
mcPCC 0.998, which is a physical-contamination signature).

Because mcFano is computed per gene against a global depth fit, a handful of
enormously overdispersed contaminant genes both crowd out the HVG list and skew
the fit for everything else, so we drop them up front rather than post-filtering.

The patterns are deliberately narrow. Genes that merely start with the same
letters but are real biology are NOT matched: HBEGF, HBP1, IGF1, IGHMBP2 (a
helicase), IGLON* (neuronal adhesion), CRYL1, CRYZ and CRYM (not lens-structural).

Caveat worth knowing: CRYAB is also HSPB5, a genuine stress chaperone expressed
in retina and glia. It is included here because in this cohort it tracks with the
other lens crystallins, but a sample where CRYAB matters biologically would need
`--exclude-contaminants` narrowed to the groups that apply.
"""

import re

CONTAMINANT_PATTERNS = {
    # HBA1/2, HBB, HBD, HBE1, HBG1/2, HBM, HBQ1, HBZ -- erythrocyte carry-over.
    "haemoglobin": re.compile(r"^HB[ABDEGMQZ]\d?$"),
    # Lens structural crystallins only (alpha, beta, gamma).
    "crystallin": re.compile(r"^CRY(A[AB]\d?|B[AB]\d|BB\d|G[A-DNS])$"),
    # Ig constant + variable segments and the J chain -- plasma-cell / serum carry-over.
    "immunoglobulin": re.compile(
        r"^(IGH[VDJ]\d.*|IGH(G\d|A\d|M|D|E)|IGK[VC].*|IGL[VCJ]\d.*|IGLL\d|JCHAIN)$"),
}

GROUPS = tuple(CONTAMINANT_PATTERNS)


def contaminant_genes(var_names, groups=GROUPS):
    """Return {group: [genes]} for the gene names matching each requested group."""
    hits = {}
    for g in groups:
        pat = CONTAMINANT_PATTERNS[g]
        hits[g] = [n for n in var_names if pat.match(n)]
    return hits


def drop_contaminants(adata, groups=GROUPS, verbose=True):
    """Drop contaminant families from `adata` in place-ish; returns the filtered adata."""
    hits = contaminant_genes(adata.var_names, groups)
    drop = {g for genes in hits.values() for g in genes}
    if verbose:
        for group, genes in hits.items():
            shown = ",".join(sorted(genes)[:8]) + ("..." if len(genes) > 8 else "")
            print(f"  drop {group:15s} n={len(genes):3d}  {shown}", flush=True)
    keep = [n for n in adata.var_names if n not in drop]
    return adata[:, keep].copy()
