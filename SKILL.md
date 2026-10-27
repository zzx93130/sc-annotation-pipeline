---
name: sc-annotation-pipeline
description: Execute the lab's standard cell-annotation pipeline for human or mouse scRNA-seq data, from raw multi-batch counts to tiered, ontology-mapped cell-type labels. Use this whenever the user wants to annotate cell types in single-cell data, process new scRNA-seq samples, run QC-to-annotation end to end, integrate batches before annotation, or says things like "annotate my dataset", "run the annotation pipeline", "process these h5ad files", "assign cell types to my samples", or hands over raw CellRanger/.h5ad output and asks what cell types are in it. Also use it when the user asks to run only part of the pipeline (just QC, just integration, just the annotation loop) — the stages are modular and this skill says how each one must be done.
metadata:
  authored_via: "claude-science"
  published_by: "agent"
  published_at: "2026-10-25T02:06:02.706Z"
  authored_by: "unknown"
  last_modified_by: "agent"
  last_modified_at: "2026-10-25T02:06:02.603Z"
---

# Single-cell cell annotation pipeline

This skill executes the lab's pipeline design (the full rationale lives in
`references/design.md` — read it before the first run in a session, and consult
the relevant section whenever a decision below needs its "why"). For a
multi-study atlas, read `references/batch_effects.md` before deciding any
batch-correction key: dataset-confined clusters have at least three distinct
causes needing three different answers, and one global correction destroys a real
population in one of them. The design's
core principles, which resolve most judgment calls that come up mid-run:

1. Annotation quality is bounded by preprocessing quality — QC, ambient
   correction, and doublet removal are annotation stages, each with a gate.
2. Never let one method decide: every label needs two independent lines of
   evidence, and disagreements are surfaced, not silently resolved.
3. Raw counts are immutable (`layers["counts"]`); every transform gets its own
   layer/obsm slot.
4. Annotate hierarchically, and **run the loop once per level**: label the
   coarse compartments on the whole object, then subset each compartment and
   run the same filter/vote/triage/subcluster loop inside it. Each level gets
   its own Leiden resolution, its own marker panel graded to that level, and
   its own Cell Ontology term and evidence tier. A marker is only as specific
   as the set of populations it is asked to discriminate — `PTPRC` separates
   immune from epithelium and says nothing about a macrophage subtype, so
   level-1 markers must not be scored in a level-3 comparison.
5. Uncertainty is an output: `unresolved` is a legitimate label; forced
   assignment is not.
6. Everything is logged to `adata.uns["pipeline_log"]`: parameters, versions,
   per-stage cell counts, and every removal.

## Before starting

Collect from the user (ask only for what's missing; proceed with stated
assumptions if they're away): the input files (.h5ad or 10x .h5 per sample), a
sample sheet (sample_id, batch, condition, species, tissue), the tissue (this
selects the CellTypist model, reference atlas, and marker dictionary), and any
lab marker YAML. Species is a config switch — never mix human and mouse in one
run. **`sample_id` is the prioritized `batch_key`** throughout.

For **human lung / IPF**, `references/lung_ipf_marker_panels.md` is a reproduced
run record: the curated per-compartment marker panels actually scored, the final
reviewed subtype labels (77 across myeloid / epithelial / T / B / NK-ILC /
fibroblast / endothelial), the genes that discriminated each subtype in an
8-study cohort, and the IPF differential-abundance signature. Use it as a
scaffold — re-grade each gene inside its compartment; markers keep their source.

Companion skills to use where noted below, when they are available in the
session: `single-cell-rna-qc` (Stage 2), `scvi-tools` (Stages 5 and 7a — its
scripts cover training, label transfer, and clustering), and
`sc-cluster-annotation-loop` (Stage 6 filter and the whole of Stage 8; its
kernel helpers are validated — prefer them over reimplementing).

## Environment gotchas

**UMAP segfaults at `import umap` (arm64 / macOS).** `umap-learn` imports its
optional `parametric_umap` submodule at package load, which imports TensorFlow;
on Apple-silicon conda builds TF's native `preload_check()` aborts (SIGSEGV,
exit 139) at import. Because it's a native abort, not a Python `ImportError`,
umap's own try/except cannot catch it and the whole process dies BEFORE any
embedding runs — so the crash looks scale-related (it fires at 15k cells and at
840k alike) but is purely an import-time fault. This is NOT a
numba/pynndescent/llvmlite problem; don't chase the optimizer, subsample, or
switch `init`.
Permanent fix: `manage_packages(mode="uninstall", environment=<env>,
packages=["tensorflow"])` — nothing in this pipeline uses TF; the crashing core
package goes and `import tensorflow` becomes a catchable error, so umap loads
and full UMAP fits run (~12 s for 15k cells).
In-session fix (no uninstall rights): `import sys;
sys.modules["tensorflow"]=None` BEFORE `import umap` turns the abort into a
catchable ImportError.
Diagnostic tell: run the embedding in a SUBPROCESS (bash `python script.py`),
not the kernel — a bare exit 139 with no Python traceback is the signature; a
`faulthandler.enable()` traceback ending in `tensorflow/python/platform/
self_check.py preload_check` confirms it.

## Stage execution

Run the stages in order. Stages 1–3 are per sample; 4+ on the pooled object.
After each stage, record the gate result in the pipeline log. When a gate
fails, stop and show the user the evidence rather than pushing through.

**Stage 0 — Input contract.** Verify raw integer counts, Ensembl IDs + symbols,
species (check MT-* vs mt-* prefixes against the sample sheet), complete sample
sheet, no duplicate barcodes.

**Stage 1 — Ambient RNA.** CellBender on the raw matrix (SoupX if only the
filtered matrix exists). Keep corrected and uncorrected counts. Gate: log
contamination per sample; flag > ~20%.

**Stage 2 — QC.** MAD-based, permissive (5 MADs counts/genes, 3 MADs + hard cap
for mito%; species-correct gene patterns). Gate: before/after plots + retention
table; halt if a sample loses > 40% of cells.

**Stage 3 — Doublets (3-method consensus).** Run scDblFinder, Scrublet, and
DoubletFinder per sample; store each call plus `n_doublet_calls`. ≥ 2/3 →
`doublet_consensus`, removed before integration; single-method calls stay
flagged for cluster-level review.

**Stage 4 — Normalization & HVGs.** Median log1p into `layers["lognorm"]`.
HVGs seurat_v3, 2000–4000, `batch_key=sample_id`; exclude mito/ribo/hb from the
HVG set.

**Stage 4b — Regress & scale (mandatory for every scanpy PCA embedding).** Every
scanpy clustering embedding runs `sc.pp.regress_out(["total_counts", "pct_counts_mt"])`
on the log-normalised HVG matrix, then `sc.pp.scale(max_value=10)`, then PCA (arpack).
This applies at level 1, at every level-2 per-compartment re-embedding, in every
subcluster and QC-loop re-embedding, and in joint atlas+reference spaces (compute the same
covariates for the reference cells). Use the sidecar `regress_scale_pca(adata_hvg)`. It
regresses in 500-gene blocks: `regress_out` fits each gene independently, so the PCs are
identical to one full call but memory stays bounded.
- Compute the covariates with `sc.pp.calculate_qc_metrics` on raw counts of **all** genes,
  before HVG subsetting.
- **Never regress cell-cycle or proliferation scores**: cycling is signal the pipeline
  needs downstream.
- Only an unregressed *control* embedding may skip this step (`keys=[]`), and it must be
  labelled as the control.
- Check that the regression worked with `qc_pc_correlation(adata)`. In practice, the max
  |r| of the top 10 PCs with log library size fell from 0.45 to 0.17 in a 35k-cell
  fibroblast test.
- scVI latents skip this stage (the model handles library size itself). Marker scoring
  and DE always use the unscaled `lognorm` layer.

**Pitfall seen in practice (LUNG_IPF atlas, 742k cells, 8 studies).** Regressing `total_counts` at the
level-1 (compartment) embedding moved 1,306 cells (0.18%, against 155 without regression) into the
wrong compartment. 1,196 of them came from one deep-sequenced study and were very shallow (median about
730 UMIs against 3,685 for that study). They then received wrong-lineage labels. After routing, run
`compartment_drift(obs, "<unregressed or manual compartment>", "compartment", batch_col="<study>")`.
If the changed cells are shallow and concentrated in one batch, add a depth floor (about 1,000 UMIs)
for that batch before routing. Regression also dissolved small substates in that run (DC1 merged into
migratory DC; IFN-response, proliferating and HSP-high T cells; MT-high alveolar macrophages). When a
substate disappears, check whether it was QC-driven before you drop it. If Leiden resolution is chosen
by a data-driven rule, re-run the rule on the regressed embedding; do not reuse the old resolution.

**Stage 5 — Integration.** scVI (batch_key=sample_id) and/or Harmony on the
Stage 4b PCA — running both is legitimate: curate on Harmony's embedding,
train scVI/scANVI for label transfer. The scVI step includes its optimization
part: hyperparameter search (n_latent, n_layers, n_hidden, dropout, ZINB vs
NB; scvi.autotune or a small grid), early stopping on validation ELBO, and
model selection on scIB-style metrics — not loss alone. Gate: batch mixing +
bio-conservation; confirm always-distinct types (B vs T) stay separated.

**Stage 6 — Clustering (iterative, per level).** Leiden at 0.2 / 0.5 / 1.0 on
the chosen latent. The **coarse resolution is not a diagnostic — it is the
level-1 annotation resolution**: pick the resolution whose clusters correspond
one-to-one with compartments (typically ~0.2), and keep the finer resolutions
for use inside each compartment subset at level 2. Loop: cluster-size barplot →
small clusters with no markers and no cycling genes in their top expressed genes
are removed → redo Stage 4 HVGs → 4b/5 → 6, until clean (cap ~5 iterations; log
every removal). Then apply the four-criterion cluster filter from
`sc-cluster-annotation-loop` (size floor, donor breadth, per-cell dual-lineage
test calibrated to this dataset's median, prior designations) with the
ambient-vs-doublet complexity diagnosis. Re-apply the same filter inside each
compartment subset — a fragment that passes the size floor globally can fail it
within its compartment.

**Stage 7a — Reference engines.** CellTypist with the **same-tissue model**
(generic Immune_All only as a flagged fallback) + scANVI/scArches transfer from
the tissue-matched atlas (HLCA, Tabula Sapiens / Tabula Muris Senis; popV
optionally as a third engine). Per-cell calls, computed once.

**Stage 7b — Markers, graded by level.** Score the user-defined,
tissue-matched marker dictionary (public DBs are scaffold only, never override
tissue-matched markers); Wilcoxon rank_genes_groups per cluster on lognorm;
check DE hits against artifact signatures. Marker z-scores are recomputed every
annotation round, never carried forward.

Every marker carries a **level**, assigned from how broadly the sources claim
it — not by hand:

| level | claimed for | example | scored at |
|---|---|---|---|
| L1 super-compartment | many types, one super-compartment | `PTPRC`, `CD53`, `CLDN5` | coarse resolution only |
| L2 compartment | ≥4 types, one compartment | `EPCAM`, `ELF3`, `CLEC7A` | coarse resolution only |
| L3 cell type | 2–3 types, one compartment | `CD3E`, `MARCO` | within-compartment |
| L4 subtype | exactly one type | `SFTPC`, `CLEC9A`, `FOXJ1` | within-compartment |
| L0 uninformative | spans >1 super-compartment | `LAMP3`, `HOPX`, `COL1A1` | excluded globally; may be reinstated inside a subset where it no longer straddles |

Two rules follow. **The level decides what a gene may be scored *against*, not
whether to score it.** At level 1 a compartment's score set is the **union of
every marker of every type it contains** — restricting to L1/L2 genes leaves
endothelial and LEC with identical sets, likewise fibroblast and SMC/pericyte,
and mesothelium with none. At level 2 the L1/L2 genes must be dropped, because
every population inside the compartment carries them and they add a constant to
each competing set. And **re-derive specificity inside each subset**: a gene's
promiscuity is relative to the populations present, so `LAMP3` is globally
ambiguous (AT2 lamellar bodies and migratory DC) yet discriminating within the
myeloid subset.

**Keep the panels small.** More markers is not more evidence: AUCell averages
over the set, so every non-discriminating gene pulls the score toward the
compartment mean and shrinks the gap between competing sets. Use **3 genes per
compartment at level 1** and **at most 5 per cell type at level 2**. Measured on
the lung atlas: 27 level-1 genes reproduced a 402-gene union panel exactly
(14/16 clusters, 98.63% of cells, questioning the same two clusters), and in the
myeloid subset a cap of 5 gave median call margin 0.810 against 0.756 using every
available gene and 0.678 at a cap of 8 — the larger cap readmits `LAMP3` to AT2
and the AT2/DC boundary degrades. Pick the cap on measured margin, not by feel.
A bigger panel can also be confidently wrong where a small one is correctly
uncertain: one 1,141-cell myeloid cluster came out Mast at margin 0.039 on 27
genes (flagged) and Lymphoid at margin 0.903 on 66 genes (a confident error).

**A state must not strip its parent of a shared marker.** `SFTPC` is claimed for
AT2 by one source and for AT2-transitional by another, and both are right — it is
an AT2-lineage gene high in both. Count marker breadth over **parent groups**,
keep the shared gene in the parent's set, and score the state on markers unique
to it. Where a state has no unique marker, mark it `unscoreable`: giving it the
parent's pan-lineage genes makes sibling states mathematically indistinguishable
(three monocyte-derived macrophage states ended up with byte-identical sets that
way). "Cannot be separated on marker evidence" is a result; an identical score
set is a bug.

**Grade support by uniqueness, not by corroboration count.** A type can carry a
subtype claim when ≥3 markers *in its scored set* are claimed for it alone —
count the unique markers, not the set size. Counting sources instead rewards
promiscuity and excludes the decisive markers: on this panel it dropped 255
subtype-unique markers across 35 cell types, leaving CD8 T cell with `CD8A` alone
(no `CD8B`, `GZMH`, `GZMK`), basal with `KRT5`/`TP63` alone (no `KRT14`,
`KRT15`, `DSG3`), mast with `MS4A2` alone (no `CPA3`, `TPSAB1`) and AT2 with
`MFSD2A`/`SLC39A8` alone (no `SFTPD`, `MUC1`). Record both columns so the choice
stays auditable.

Use a **rank-based scorer, not a control-gene one**. AUCell (`decoupler`
`dc.mt.aucell`) or UCell score each cell from its own gene ranking and are
therefore independent of dataset composition; Seurat's `AddModuleScore` and
`scanpy`'s `score_genes` normalise against a control set drawn from the
supplied matrix, so on a marker-only matrix the background is meaningless.
Score in the **full gene space**, never a marker-only subset.

**Stage 8 — Vote, triage, loop — once per hierarchy level.**

*Level 1, whole object.* Cluster coarse (Leiden ~0.2), then follow
`sc-cluster-annotation-loop` to label the compartments: epithelial, endothelial,
lymphatic endothelial, fibroblast, smooth muscle / pericyte, myeloid, lymphoid,
mast, mesothelial. Score the union panel per compartment. This level is where
the evidence is strongest and it should reach near-complete agreement — on a
840k-cell lung atlas, AUCell matched the reference consensus on 98.6% of cells,
with the residue being one cluster where the panel resolved LEC and the
reference said only "endothelium", and one 1,141-cell cluster at margin 0.022
that was left unresolved rather than force-called. Report the margin with every
call and flag rather than assign below the floor. Freeze the compartment label
and tier before going deeper.

*Level 2, per compartment.* **Subset each compartment from raw counts and run
the standard scanpy pipeline inside it** — CP10K → log1p → HVGs selected *within
the subset* → `regress_out(total_counts, pct_counts_mt)` → scale → PCA → batch correction → neighbours → Leiden (~0.5–1.2) →
UMAP. Then re-derive marker specificity within the subset and run the same
filter → vote → triage → subcluster loop for cell types and subtypes.

Re-clustering from counts rather than restricting the parent latent is
load-bearing, not a preference. Whole-atlas HVG selection is dominated by the
genes that separate *compartments*, so the genes separating populations *within*
one compartment are frequently absent from it: on the lung atlas, CD4 T cell and
B cell had **zero** of their subtype markers in the atlas HVG set, meaning the
integrated latent had no axis along which those two could separate at any Leiden
resolution. Subset HVG selection recovered them. It also gives the
under-clustered compartments enough clusters to support their panels for the
first time — endothelium went from 2 usable clusters to 14, fibroblast from 1 to
14.

The trade-off, which must be stated in the output rather than hidden: sub-labels
are no longer positioned in the parent embedding, so a level-2 UMAP is only
comparable within its own compartment. Carry both the level-1 and level-2
cluster ids per cell so the join back is explicit. (Restricting the parent latent
is the alternative — it keeps one shared embedding, and is the right choice when
the level-2 populations are known to separate on atlas-level HVGs. Log which was
used.) Subsetting is what buys the resolution either way: median call margin
across myeloid clusters rose 0.527 → 0.773 when the panel was re-derived inside
the myeloid subset, and 19 of 60 calls changed.

At both levels: evidence table per cluster (every voter distinguishable, legacy
labels reported separately) → parent-level rollup before calling anything
ambiguous → tiers A–D with predicate-generated criteria → subcluster tier C only
(same latent, resolution ~0.5, re-filter the fragments; drop-and-record, never
merge back) → repeat until tier C is empty. Report tier-C shrinkage AND tier-D
growth each round. Freeze label + tier + CL term at every level, and never
promote a cell to a subtype whose compartment call is not itself frozen.

**Stage 8b — Calling rules.** Three rules govern whether a scored cluster earns
its label. All three were forced by a failure they fix.

*Support floor, non-negotiable.* A call is only admissible if its own set is
actually expressed: the best-detected gene in the called set must clear ~20% of
the cluster's cells. A margin computed over sets that are all absent is
meaningless — a standardised score yields a verdict from two negative values
without ever testing expression.

*Then EITHER strong detection OR a clear margin.* Above the floor, resolve when
the called set's best gene is detected in ≥40% of cells **or** the margin to the
runner-up is ≥0.75. Both paths are needed and neither alone works: a
detection-only cut strips dropout-prone single-gene sets (`CD4` detects in only
7–34% of the clusters it scores, `FOXP3` in 9%, and CD4 RNA is genuinely
under-captured), while a large margin over a set that is not expressed is
confidently wrong. Report which path admitted each call.

*Required defining genes veto the call.* Some sets are named for a gene that
must be present — a channel marker shared with other lineages is not enough to
call an ionocyte. Give those sets an explicit `required_genes` list; if none of
the required genes clears the floor, the cluster goes unresolved with the missing
genes named, and falls through to the best set that *is* eligible.

Two invariants on required genes, both learned the hard way:

- **A required gene must be a MEMBER of the set it gates.** Detection is measured
  only over the scored set, so a required gene named outside it reads as absent
  and vetoes unconditionally — on first use this blanked six CD8 and three
  monocyte clusters, with an empty detection string as the only clue. Assert it.
- **Every gate raises the chance an UNGATED set becomes a sink.** Clusters
  excluded from their own set fall through to whatever is still ungated, and a
  set containing a pan-lineage gene will absorb them: gating two capillary sets
  sent four clusters to `Tip cell` (which contained pan-endothelial `KDR`) where
  `ESM1` read 1–13%. After any gating change, check whether an ungated set gained
  clusters it did not previously hold. Before gating a sink set, confirm the
  population exists at all — max `ESM1` across all 19 endothelial clusters was
  13.3%, below the floor, so gating it removed an unsatisfiable set rather than
  deleting cells.

**Stage 8c — Cell state as a second axis.** Proliferation, interferon response
and stress/heat-shock are **states a cell of any identity can be in**, not
identities. Hold their gene sets out of identity scoring and report them in a
separate column; a proliferating T cell is still a T cell.

Two things make this work:

- **Relative, not absolute, thresholds.** An absolute detection floor sits inside
  one lineage's baseline and never fires elsewhere. Require the state's median
  detection to clear both a floor (~30%) *and* a fold over **that lineage's own
  median** for that state (~2.5×). Stress needs a lower fold (~2.0) because it is
  a graded program rather than on/off — epithelium runs a 40.7% stress baseline
  where T/NK runs 19.7%, so one absolute cut would condemn every epithelial
  cluster or spare every real case.
- **Exclusive, not independent.** Test states against each other and assign at
  most one, or a cluster comes back carrying three at once.

**Canonicalise lineage-flavoured set names.** A panel often carries variants of
one program per compartment (`Proliferating myeloid`, `EC IFN_response`). Map them
onto a single program label before the exclusive assignment, or the same biology
is reported under different names per lineage and the variant can outrank the
canonical set. Keep the raw set that fired in its own column.

Restrict stress to **inducible** chaperones (`HSPA1A`, `HSPA1B`, `DNAJB1`,
`HSPH1`). Including constitutive ones (`HSP90AA1`, `HSPB1`) puts the baseline
near half of all cells and the fold rule can never fire. Watch the proliferation
set for ubiquitous chromatin genes — including one inflated the median enough to
mark several clusters proliferating spuriously.

The fold rule can pick a state with *lower* absolute detection when that state's
baseline is lower, and that is correct: one epithelial cluster read stress 47.4%
(1.16× its lineage median, fails) and IFN 32.6% (4.59×, passes), and IFN is the
right call.

**Stage 9 — Validate the annotation, and give the user something to review.**
Recomputed (never subset) marker dot plots and a cluster UMAP; specificity and abundance plots kept
distinct; agreement rates reported with their matching rule; composition sanity
per sample; no cluster defined by mito/ribo/stress/ambient genes.

**The review deliverable is four plots per lineage, not one, and they ship on
every run — not only when the user asks.** A finished annotation is not
reviewable until they exist. Each dot-plot row is labelled with the cluster id,
its assigned label, its status, its size and its margin (unresolved rows in red);
the UMAP direct-labels each cluster in surrounding whitespace:

| view | what it is | answers |
|---|---|---|
| **UMAP coloured by cluster** | the clustering embedding, one colour per cluster | is the partition coherent — does each called cluster hold its own territory, and do any two labels overlap |
| **panel markers** | the marker panel, grouped by the cell type each gene was written for | does the called set actually fire, and does it fire *here* rather than next door |
| **cluster-vs-rest DE** | each cluster's top 10 genes by Wilcoxon score (cluster-vs-rest), adjusted p-value shown alongside | what actually distinguishes this cluster from its lineage neighbours |
| **status genes** | HSP/stress, IFN-response and proliferation sets, grouped by status | is the status call real, and is a thin-margin identity call explained by a status program dominating the cluster |

Ship all four; they fail differently and a reader cannot substitute one for
another. If the embedding does not resolve a called population (recovered by
targeted subclustering, or partition and layout from different spaces), say so
on the UMAP. The DE/status plots are selected without the panel, but that does
not make the plotted genes disjoint from it — a cluster's own marker can top
its DE (Mast tryptases, LEC CCL21), which is corroboration; never caption "none
of these are panel genes" without computing the overlap.

**The DE view must be DE, not abundance.** A top-10 *highest-expressed*
ranking is dominated by housekeeping transcripts and carries almost no
information about identity — measured across six lineages on the lung atlas, 48
genes appearing in some cluster's abundance top-10 were shared by half or more of
that lineage's clusters, and 78.5% of those entries were housekeeping (`MALAT1`,
`B2M`, `FTH1`, `TMSB4X`, ribosomal and mitochondrial transcripts). The same
top-10 taken from the *filtered DE* list has **zero** genes shared by half the
clusters, a median sharing of 1 cluster, and 0.5% housekeeping. Abundance answers
a different question and can be shown on request; it is not the review view.

A **thin DE row is a result, never a rendering gap** — do not pad it. Four of 116
clusters had fewer than 10 passing DE genes, and in each case that is the
finding: one 33,392-cell alveolar-macrophage cluster has 5, meaning almost
nothing separates it from the other alveolar clusters at this resolution.

The DE view catches what the panel view structurally cannot. One 467-cell
B-lineage cluster's top DE genes were secretoglobins (`SCGB3A1` 78%, `SCGB1A1`
73%, `BPIFB1` 26%, `SLPI` 33%) — airway-secretory ambient soup, invisible to the
panel view because no B-lineage set contains those genes, and invisible to the
min-scored cross-lineage test because soup is *secreted* transcripts while the
level-1 sets are membrane/TF genes. Note the reading it also invites: those cells
are real B cells (`MS4A1` 82%, `CD79A` 79%) carrying ambient burden, not an
epithelial cluster to remove. Flag ambient burden as a column beside this view so
a secreted program is legible as contamination rather than identity.

Conventions, all four views: rows sorted so clusters sharing a call sit together;
**margin in every row label** and **unresolved rows in red**, so a weak call is
visible without cross-referencing a table; dot size = % of cells detecting;
colour = z of mean expression **across the rows shown**, stated as valid within
that panel only; vertical rules between gene blocks. Compute statistics on the
clusters shown — never subset a previously computed matrix.

Whichever view is shown, **print the per-gene sharing count beside it** — how
many clusters in the lineage carry each gene in their own top 10. It is what
makes "distinctive" checkable rather than asserted, and it is the number that
separates the DE view (median 1) from an abundance view (median 2–4, with 48
genes shared by half or more).

**Stage 9b — DE markers, with a specificity term.** Cluster-vs-rest Wilcoxon,
then filter on significance **and** effect size **and** specificity. Significance
alone discriminates nothing at these cell counts: on the *smallest* lineage
11.3% of all gene × cluster tests clear an adjusted p of 0.05, rising to 31.4% on
the largest — thousands of genes per cluster. Require adjusted p < 0.05, detection
in ≥25% of the cluster, log fold-change ≥ 0.5, and a **detection difference**
(pct_in − pct_out) ≥ 15 points; rank by that difference so the ordering and the
plot agree. This moves housekeeping content of the top genes from ~90% down to
around 1%.

Two things this view is good for and one it is not. It identifies clusters the
panel cannot — a small B cluster with no panel call carried an unmistakable
germinal-centre signature (`RGS13` 79/2, `LRMP` 78/15, `HMCES` 57/7). It also
surfaces a **state** as a marker set whenever the state defines the cluster,
which is why identity is called from the panel and the state reported separately.
It is *not* a substitute for the panel view: a cluster whose DE list is topped by
its identity program will show no state gene in the top 15 even when four state
genes pass every filter — so corroborate a state against the **passing** gene
set, not a top-N slice. Measured: 14 of 14 state calls had ≥1 state gene passing;
only 10 of 14 had one in the DE top 15.

**Provenance is not optional here.** Re-clustering renumbers clusters, so a DE
table computed on a previous partition will join onto the new ids and produce
lineage-inappropriate genes beside correct detection values. Matching ids are
*not* evidence the table is current — that is exactly what makes the failure
quiet. Assert that the DE cluster ids are a subset of the current label table,
and measure the overlap between the two independent gene views (panel-scored
detection vs cluster-vs-rest) per cluster; a median of ~7 of 10 is healthy.
Zero overlap is expected for a large homogeneous cluster whose cluster-vs-rest
comparison finds only ribosomal genes while the filtered DE view returns the
canonical program — verify those individually rather than assuming either
outcome.

**Stage 10 — Outputs.** `annotated.h5ad` (all layers, latents, resolutions,
labels + tiers + CL terms, flags, full log), the exclusion table and
dropped-subcluster table, the HTML report, and the versioned config + marker
dictionary. Deliver the report and h5ad to the user; the tier column is
load-bearing for any downstream differential analysis.

Also required, because they are what a user actually reviews:

- **The four dot plots per lineage** (Stage 9), and a **per-cluster review
  table** with one row per cluster: label, state, size, margin, which call path
  admitted it, per-gene detection for the called set, required-gene detection,
  the unresolved reason where applicable, the DE top genes with their in/out
  detection, and the top-10 expressed genes. One row per cluster is what makes
  a claim checkable without re-running anything.
- **Retain every cell**, including those removed as contamination, flagged with
  the lineage they were removed from and the reason. Cells leaving an analysis
  must be accounted for somewhere a reader can check.
- **A decisions record** naming every rule that was changed mid-run, what forced
  the change, and any claim that was withdrawn. A pipeline whose rules moved is
  fine; one whose rules moved silently is not.

## Judgment calls

- **Never build a marker panel from the clusters it will judge.** A panel taken
  from a previous run's cluster-vs-rest DE cannot falsify the labels it came
  from. It happened in this lab: 98 of 111 genes in a panel described as
  "authoritative and tissue-matched" were its own clusters' DE hits, 109 of 111
  were rank-1 for the type they marked, and it omitted two thirds of textbook
  lung markers. Panels come from published atlases, ontology-curated resources,
  or the bench — and each gene keeps its source. A previous run's panel is a
  comparison column, never evidence.
- **Never write marker genes from memory.** Every gene traces to a retrieved
  record with its source; a recalled list is not a source even when the genes
  are individually plausible.
- **A corroboration count is not specificity.** Requiring ≥2 independent sources
  rewards promiscuity: pan-lineage genes are listed by everyone for everything
  and score as well-supported, while a decisive single-source marker (`SFTPC`
  for AT2) gets dropped. Grade by level (Stage 7b) and keep decisive
  single-source markers, flagged.
- User asks for a single stage: run it under this skill's rules, still log to
  `pipeline_log`, and tell them which upstream assumptions you're making.
- Small dataset (< ~20k cells): Harmony-only integration is fine; skip scVI
  hyperparameter search, keep the gate.
- No GPU for CellBender/scVI: say so, use SoupX and either CPU scVI (small
  data) or Harmony; record the deviation in the log.
- Missing sample sheet fields: infer what is safely inferable (species from
  gene prefixes), ask or halt for the rest (condition, batch cannot be
  guessed).
- Anything that deviates from `references/design.md` gets a one-line entry in
  the pipeline log with the reason.
- **Test at the resolution you annotate at.** Every removal and quality rule must
  run on the partition the labels come from. This blind spot bit twice on one
  atlas: contamination was tested at Leiden 0.2 while subtyping used the finer
  resolution chosen afterwards, so clusters that only separate at fine resolution
  were never tested; and the ribosomal/low-complexity rule ran only inside the
  contamination rounds and never on the final partition at all.
- **Validate a rule change against saved per-cluster tables before spending
  compute.** Three full lineage re-runs were discarded for defects that a
  seconds-long check against the existing detection table would have caught. The
  prediction also tells you when a change does nothing: an edit adding two
  transcription factors changed no call, because the support floor keys on a
  set's *best*-detected gene and the shared marker cleared it alone.
- **Compare cluster-by-cluster against the previous run, never headline counts.**
  A run showed both a higher labelled fraction and an extra subtype while two NK
  clusters had been silently relabelled by an over-broad new set. In another, net
  subtype count was identical while a rare dendritic population had been absorbed
  into a much larger one. The headline hides exactly the regressions worth
  catching.
- **A subsample fixes which knob to turn, not what the metric will read.**
  Subsampling shrinks every cluster and so inflates any per-cluster purity or
  dominance measure — a four-arm test on a subsample picked the right
  configuration and mispredicted its effect size by 5×, in the wrong direction.
  Use it to choose; re-measure on the full data.
- **Do not tune a threshold to capture specific clusters.** This is the pattern to
  refuse most consistently: a floor lowered to admit one population, a fit
  validated on the same cases that motivated it, a rule "verified" on its own
  training set. When no rule isolates what the user wants removed, implement it
  as a **named exclusion carrying its evidence**, flagged as their judgement, with
  the rules you tested and rejected recorded — and a drift check that withholds
  the exclusion if the cluster's size moves more than ~10%, since a cluster id
  does not survive a re-partition.

## Helpers (`kernel.py`, auto-loaded)

| function | does |
|---|---|
| `regress_scale_pca(adata_hvg, keys=None, n_pcs=50, block_genes=500)` | default PCA path: blocked `regress_out(total_counts, pct_counts_mt)` → `scale(max_value=10)` → arpack PCA into `obsm["X_pca"]`. It matches the LUNG_IPF `ipf_annot.embed.scaled_pca` exactly (max PC difference 0.0). |
| `qc_pc_correlation(adata, rep="X_pca")` | max \|r\| of the top 10 PCs with log total_counts and pct_counts_mt (before/after check) |
| `compartment_drift(obs, before_col, after_col, batch_col=)` | cells whose compartment changes after regression: their depth, batch and top transitions |

