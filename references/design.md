# Cell Annotation Pipeline — Design Document

**Scope:** Human and mouse scRNA-seq · raw counts · multiple samples/batches
**Strategy:** Hybrid annotation — automated reference-based label transfer, cross-checked against marker-based manual curation
**Stack:** Python / scverse ecosystem (scanpy, scvi-tools, CellTypist), following [sc-best-practices.org](https://www.sc-best-practices.org)

---

## 1. Design principles

1. **Annotation quality is bounded by preprocessing quality.** Most "annotation problems" are actually ambient RNA, doublets, or batch effects that survived upstream. The pipeline therefore treats QC and integration as first-class annotation stages, each with an explicit quality gate.
2. **Never let one method decide.** Every label must be supported by at least two independent lines of evidence (e.g., reference transfer *and* canonical markers). Disagreements are surfaced, not silently resolved.
3. **Keep raw counts immutable.** Raw counts are preserved in `adata.layers["counts"]` forever; every transformation lives in a separate layer or `obsm` slot so any stage can be re-run.
4. **Annotate hierarchically.** Coarse compartments first (immune / epithelial / stromal / endothelial), then subtypes within each compartment. Errors at coarse level are catastrophic; errors at fine level are expected and manageable.
5. **Uncertainty is an output, not a failure.** Cells the pipeline cannot confidently label are flagged `unresolved`, never force-assigned.
6. **Everything is provenance-tracked.** Each stage writes its parameters, software versions, and per-stage cell counts into `adata.uns["pipeline_log"]`.

---

## 2. Pipeline overview

```
Stage 0  Input contract & validation
Stage 1  Ambient RNA correction            (per sample)
Stage 2  QC filtering (MAD-based)          (per sample)
Stage 3  Doublet detection                 (per sample; 3-method consensus)
Stage 4  Normalization & feature selection (pooled)
Stage 4b Regress-out & scaling             (regress_out + scale; for the PCA path)
Stage 5  Integration / batch correction    (scVI or Harmony)
Stage 6  Clustering & embedding            (Leiden, multi-resolution; loops
         back through Stages 4–6 until no junk small clusters remain)
Stage 7a Automated annotation              (CellTypist + scANVI label transfer)
Stage 7b Marker-based evidence             (scores + DE)
Stage 8  Vote → triage into tiers A–D → subcluster tier C → re-filter →
         freeze label + tier (loop until tier C is empty)
Stage 9  Validation & QC of the annotation itself
Stage 10 Outputs & reporting
```

Stages 1–3 run **per sample** (batch effects corrupt cross-sample QC); Stages 4+ run on the **pooled** object.

---

## 3. Stage-by-stage design

### Stage 0 — Input contract

The pipeline accepts one `.h5ad` per sample (or a 10x `.h5` converted on ingest) and a sample sheet (CSV) declaring: `sample_id`, `batch` (chemistry/run/site), `condition`, `species`, tissue, and expected cell-count range.

Validation gate before anything runs:

- `X` contains raw integer counts (checked, not assumed).
- Gene identifiers are Ensembl IDs with symbols in `var`; species inferred and confirmed against the sample sheet (human `MT-*` vs mouse `mt-*` prefixes are a cheap sanity check).
- No duplicated barcodes across the intended union; sample sheet fields complete.

**Species handling:** human and mouse are processed as *separate pipeline runs* sharing one codebase. Species is a config switch that selects gene patterns (`MT-`/`mt-`, `RPS|RPL`/`Rps|Rpl`, `HB[AB]`/`Hb[ab]`), the reference atlas, and the marker dictionary. Do not mix species in one integration unless the scientific question demands it (then ortholog-map first — out of scope for the default path).

### Stage 1 — Ambient RNA correction

Free-floating mRNA from lysed cells contaminates every droplet and is the leading cause of spurious marker expression (e.g., hemoglobin or lysozyme "expression" in T cells), which directly corrupts annotation.

- **Tool:** CellBender (`remove-background`) on the *raw* (unfiltered) matrix per sample; SoupX is the lighter-weight alternative when the raw matrix is unavailable.
- Store both corrected and uncorrected counts; downstream stages use corrected counts, but the report compares marker specificity before/after.
- **Gate:** contamination fraction per sample logged; samples > ~20% flagged for review rather than silently passed.

### Stage 2 — QC filtering

MAD-based (median absolute deviation) adaptive thresholds per sample, per scverse best practices — deliberately **permissive**, because aggressive QC removes rare cell types before annotation ever sees them.

- Metrics: total counts, genes detected, % mitochondrial, % ribosomal, % hemoglobin.
- Defaults: 5 MADs for counts/genes, 3 MADs + hard cap (~8–10% human PBMC-like tissue; tissue-dependent — neurons/cardiomyocytes tolerate more) for mito.
- **Gate:** before/after violin plots and per-sample retention table. A sample losing > 40% of cells halts the pipeline for human review.

### Stage 3 — Doublet detection

Doublets masquerade as "novel intermediate cell types" and are the classic annotation trap.

- **Tools (consensus approach):** run **scDblFinder** (via `scanpy` external and R bridge), **Scrublet**, and **DoubletFinder** (R) per sample, each with the expected doublet rate scaled to loaded cell count (~0.8% per 1,000 cells for 10x). Store each method's score and binary call in `obs` (`dblfinder_call`, `scrublet_call`, `doubletfinder_call`) plus a consensus count `n_doublet_calls` (0–3).
- **Consensus rule:** a cell called a doublet by **≥ 2 of the 3 methods** is marked `doublet_consensus = True` — these are high-confidence doublets. Single-method calls are kept but flagged `doublet_candidate` for cluster-level review rather than removed outright, since each method has distinct false-positive profiles (Scrublet struggles with homotypic doublets; simulation-based methods are sensitive to the assumed doublet rate).
- **Policy:** score every cell; *flag* at this stage, *remove* consensus doublets before integration. After clustering (Stage 6), any cluster enriched for `doublet_candidate` cells (or with median scores in the top decile) plus dual-lineage markers is reviewed as a probable doublet cluster.

### Stage 4 — Normalization & feature selection

- Pool samples into one object; keep `layers["counts"]`.
- Median-count log1p normalization (`sc.pp.normalize_total` + `log1p`) into `layers["lognorm"]` — used for marker scoring, DE, and visualization. (scVI consumes raw counts directly, so normalization choice does not affect integration.)
- HVGs: `flavor="seurat_v3"`, `n_top_genes=2000–4000`, computed **with `batch_key`** so batch-specific genes don't dominate — with **`sample_id` prioritized as the `batch_key`** (fall back to a coarser run/chemistry column only when per-sample groups are too small for stable HVG estimation). Exclude mito/ribo/hemoglobin and cell-cycle genes from the HVG set by default (configurable — keep cell-cycle genes if proliferation states are a target of annotation).

### Stage 4b — Regress-out & scaling (PCA path)

For **every** scanpy PCA-based embedding (the Harmony path in Stage 5, level-1 routing, per-compartment reclustering in Stage 8, subcluster/QC-loop re-embeddings, and joint atlas+reference spaces), always apply the classic scanpy pair after HVG subsetting (lab default; implemented by the `regress_scale_pca` sidecar):

- `sc.pp.regress_out` on unwanted covariates — `total_counts` and `pct_counts_mt`. **Do not regress out proliferation/cell-cycle genes** — cycling status is biological signal the annotation needs (the Stage 6 small-cluster check explicitly keeps cycling clusters, and proliferating states can be annotation targets).
- `sc.pp.scale(max_value=10)` so PCA isn't dominated by highly expressed genes; store the result in `layers["scaled"]` and run PCA on it.

Two guardrails: keep `lognorm` (unscaled) as the layer for marker scoring, DE, and dot plots — scaled values distort fold-change interpretation; and note that the scVI path **skips this stage by design**, since scVI models raw counts directly and handles technical covariates via `setup_anndata` (`continuous_covariate_keys` for mito%/counts if needed). Regression choices are recorded in the config, since regressing out a covariate that correlates with a real cell state (e.g., cell cycle in progenitor-rich tissue) can erase biology.

### Stage 5 — Integration

Two co-supported methods, selected by config:

- **scVI** on raw counts of HVGs with `batch_key` = `sample_id` (prioritized, consistent with Stage 4's HVG computation). The scVI step includes an explicit **optimization part** before the model is accepted:
  - *Hyperparameter search* over the settings that matter most for integration quality — `n_latent` (10–50), `n_layers` (1–2), `n_hidden`, dropout rate, and gene likelihood (ZINB vs NB) — using `scvi.autotune` (Ray Tune based) when compute allows, or a small manual grid otherwise.
  - *Training optimization*: early stopping on validation ELBO (`early_stopping=True` with a held-out validation split), monitoring the train/validation loss curves for divergence or under-training, and adjusting learning rate / `max_epochs` when the ELBO has not plateaued.
  - *Selection criterion*: candidate models are compared on the Stage 5 gate's scIB-style metrics (batch mixing + bio-conservation), not on the loss alone — the best-ELBO model is not always the best integration. The chosen configuration and its search results are recorded in `uns["pipeline_log"]` so the optimization is reproducible. Robust across batch severity, scales well, benchmarks at/near the top (Luecken et al. benchmark), and yields a latent space that Stage 7's scANVI label transfer can be warm-started from — this is why it remains the default when reference-based annotation is central.
- **Harmony** on the PCA from Stage 4b (regress_out + scale → PCA → `sc.external.pp.harmony_integrate`). Fast, no GPU/training required, and in practice often produces cleaner, more interpretable UMAPs with crisper cluster separation — a real advantage for the manual-curation and figure-generation steps. Preferred for small-to-medium datasets and whenever visualization quality drives the review process.

It is legitimate (and common) to run **both**: cluster and curate on the Harmony embedding for its visualization quality, while still training scVI/scANVI for the label-transfer engine in Stage 7a — the two latent spaces coexist in `obsm` (`X_harmony`, `X_scVI`) and Stage 8's evidence table notes which embedding each result came from.

- **Other fallbacks:** scANVI/sysVI if batches span technologies (e.g., 10x v2 vs v3 vs Smart-seq). Reach for these only when the Stage 5 gate shows scVI leaving technology-driven separation in the latent space.
- **Gate:** integration is *evaluated, not trusted* — scIB-style metrics (kBET/iLISI for mixing, cLISI/silhouette on preliminary labels for bio-conservation) plus UMAP colored by batch. Over-integration (distinct cell types merged) is checked explicitly by confirming that a handful of always-distinct populations (e.g., B vs T cells) remain separated.

### Stage 6 — Clustering & embedding (iterative loop with Stages 4–5)

- kNN graph on the chosen integrated latent space (scVI or Harmony); UMAP for visualization only (never for quantitative decisions).
- **Leiden at multiple resolutions** (0.2 / 0.5 / 1.0), all stored. Coarse resolutions anchor compartment-level annotation; fine resolutions support subtype work. A clustering tree (cluster membership across resolutions) documents the hierarchy.
- No attempt to find "the correct resolution" — resolution choice is deferred to the annotation hierarchy (Stage 8).

**Small-cluster cleanup loop.** Clustering is not run once — it iterates with Stages 4–6 until no junk micro-clusters remain:

1. **Check cluster sizes** and render them as a barplot (cells per cluster, per resolution); this plot goes into the report for every iteration.
2. **For small clusters** (below a configurable threshold, e.g. < 0.1–0.5% of cells or < ~50 cells), inspect their highly expressed genes (`rank_genes_groups` on the cluster vs rest). If the top genes contain **no cell-type markers and no cycling/proliferation-related genes** (i.e., the cluster is defined only by ambient, stress, mito/ribo, or otherwise uninformative genes), remove those cells. Small clusters that *do* show real markers or a cycling signature are kept — they are candidate rare populations or proliferating states, not junk.
3. **Re-run from the HVG step**: Stage 4 (HVG re-selection on the remaining cells) → Stage 4b/5 (re-scale/re-integrate) → Stage 6 (re-cluster). Repeat until no removable small cluster remains.

**Cluster-level filter (gate into annotation).** Independent of the junk-gene loop above, every cluster entering Stage 7 passes a four-criterion filter, per the lab's `sc-cluster-annotation-loop` procedure. Any one criterion is sufficient to exclude, but the criteria are evaluated independently and the one that fired is recorded per cluster — they mean different things and a reader will want to know which applied:

1. **Size floor** (default < 300 cells) — not because small is uninteresting, but because such a cluster cannot support donor-level statistics downstream. The floor is stated in the config, not pretended to be principled.
2. **Donor breadth** — fewer than ~5 donors per condition means no contrast is testable, whatever the cell count. Size and testability are *not* interchangeable: a 179-cell cluster can pass the donor test while a 226-cell one fails it — donor breadth is the binding constraint, cell count a proxy.
3. **Per-cell dual-lineage mixing** — measured per cell, never by cluster-mean marker scores (which are diluted by off genes and inflated by single promiscuous ones). Count, per cell, how many canonical markers of each lineage are detected; a cell is dual with ≥2 of each. Flag a cluster only when **both** thresholds fire: prevalence ratio (2nd lineage / 1st) ≥ ~0.80 **and** ≥ 60% of cells dual. The absolute threshold depends on the marker panel — ambient RNA gives every cluster a baseline dual fraction — so calibrate it against the median observed across this dataset's own clusters under this panel, and report that median; the ratio is the more portable criterion. (For scale: genuine doublet clusters have measured ~90% dual at ratio ~0.95 while clean clusters sat at a third of that with ratios < 0.35.)
4. **Prior technical designation** from an earlier round (ambient mixture, doublet, unassigned) — retrieve the *recorded reason*, not just the flag. A per-cell quality gate cannot license overturning an exclusion made on other grounds (e.g., confinement to one batch is a whole-cluster property that a per-cell gate says nothing about).

**Ambient vs doublet diagnosis.** Dual-lineage clusters are diagnosed, not just excluded — the mechanism tells the user something about the experiment. Ambient contamination shows **low library complexity** (few genes at low depth) with a doublet score near the atlas median; true doublets show normal complexity with elevated doublet scores (cross-check against the Stage 3 consensus flags). Compare each suspect against a *clean cluster of the same cell type*, not the atlas average — depth varies by cell type. The diagnosis goes into the exclusion record.

Convergence guardrails: cap iterations (e.g., 5) to prevent endless churn; log every removed cluster (size, top genes, iteration) in `uns["pipeline_log"]` so the cleanup is fully auditable; and if successive iterations keep producing new small junk clusters, treat it as a symptom of upstream problems (ambient contamination, under-filtered debris) and revisit Stages 1–2 rather than looping forever.

### Stage 7a — Automated annotation (two independent engines)

1. **CellTypist** — fast logistic-regression classifiers with curated immune/tissue models. Human: **always prefer the model trained on the same tissue** as the samples being annotated (e.g., a lung model for lung samples, a gut model for gut samples — browse the CellTypist model catalog by tissue); fall back to the generic `Immune_All_Low`/`High` models only when no same-tissue model exists, and flag those labels as lower-confidence; mouse: mouse models where available, else orthology-mapped human models (flagged as lower-confidence). Run with `majority_voting=True` over the Leiden clustering. Outputs per-cell labels + probabilities.
2. **scANVI label transfer from a reference atlas** — Human Cell Atlas / CELLxGENE census reference for the tissue (e.g., HLCA for lung, Tabula Sapiens generally; Tabula Muris Senis for mouse). Warm-start scANVI from the Stage 5 scVI model, treat query cells as unlabeled, transfer labels with per-cell soft probabilities. Where a pre-trained scArches reference exists (e.g., HLCA), map into it instead of retraining.
   *Optional third engine:* popV (consensus over many predictors, [Nature Genetics 2024](https://www.nature.com/articles/s41588-024-01993-3)) when maximum robustness matters more than runtime — its per-cell consensus score slots directly into the Stage 8 confidence logic.

Both engines annotate **cells**, but decisions in Stage 8 are made at **cluster** level; per-cell labels are aggregated to cluster majority + purity. Reference calls and doublet scores are computed **once, upstream, covering every cell** — that is what makes the Stage 8 loop cheap: each round's vote is a groupby over existing per-cell calls, not a classifier re-run. The one exception is the marker source: marker z-scores are cross-cluster statistics and must be **recomputed every round**, since carrying a previous round's scores forward scales each gene against clusters that no longer exist.

### Stage 7b — Marker-based evidence

Independent of any reference:

- **Curated marker dictionary** per species/tissue (versioned YAML in the repo — this is the pipeline's most valuable long-term asset). Sources, in order of trust: **user-defined marker sets** curated by the lab, prioritizing markers validated in the *same tissue* (from the lab's own previous datasets or published same-tissue studies), then lab knowledge. Public databases (PanglaoDB, CellMarker) are used only as a starting scaffold and are not treated as reliable on their own — tissue-matched markers always take precedence when they conflict.
- `sc.tl.score_genes` per marker set → per-cluster score matrix / dot plot.
- **Data-driven markers:** `sc.tl.rank_genes_groups` (Wilcoxon) per cluster on `lognorm`; top markers checked against the dictionary and against known artifact signatures (stress/dissociation genes `FOS/JUN/HSP*`, ambient signatures from Stage 1).

### Stage 8 — Vote, triage & iterative subclustering (tiers A–D)

The heart of the hybrid design, structured as the loop from the lab's `sc-cluster-annotation-loop` procedure: **cluster → filter → vote → triage → subcluster tier C → re-filter → freeze**. The organizing idea: a cluster label is a claim, and claims vary in strength — so every frozen label carries an **evidence tier** recording how well the evidence converged, instead of all labels being treated as equal.

**The vote.** For each cluster, assemble an **evidence table** with every source kept distinguishable: CellTypist majority label + confidence; scANVI label + median transfer probability; popV consensus (when run); marker-score hits (recomputed this round); top DE genes; doublet score; sample/donor composition. Any prior/legacy label is kept as one voter but **reported separately** — agreement with it is partly built in, so where it is the minority view it must not carry a cluster it no longer describes.

**A split vote usually means granularity, not disagreement.** When reference engines disagree at leaf level, map their labels to a common parent and re-test before doing anything else: if the vote concentrates at parent level (e.g., three fibroblast subtypes summing to a decisive fibroblast call), assign at parent level — subclustering there would chase reference granularity, not biology. Do **not** reach for "the reference vocabulary lacks this cell type" without querying the model vocabularies directly; that explanation is usually available and usually wrong.

**Triage tiers:**

| Tier | Criterion | Supports |
|---|---|---|
| **A** | References concentrate at leaf level **and** internal label purity is high | Subtype-level claims |
| **B** | References concentrate at **parent** level, parents agree | Parent-level claims; subtype only with marker support |
| **C** | **Low internal purity** — the cluster splits its own cells across labels | Nothing yet: subcluster, re-vote |
| **D** | Purity adequate but references disagree with each other | Compartment-level claims only |

Tier C is about **internal coherence, not classifier disagreement** — a cluster can have references agreeing at parent level and still be tier C because its own cells carry different labels. Subclustering can act on internal heterogeneity; it cannot act on two references disagreeing. Per cluster, emit a string naming the test that placed it in its tier, **generated from the predicate, not hand-written** — prose about the code eventually contradicts the code. The former `unresolved` category maps onto this scheme: a clean cluster no reference or marker set fits is tier D at compartment level, or `unresolved` outright when even the compartment is unsupported (candidate rare/novel population).

**Subclustering (tier C only), three constraints:**

- **Reuse the same latent** — rebuild the neighbour graph on the existing scVI/Harmony embedding restricted to those cells; recomputing an embedding per subcluster makes sub-labels incomparable with the rest of the atlas.
- **Go coarse** (resolution ~0.5) — the goal is to separate the populations the cluster is mixing, not re-atomize it.
- **Re-apply the Stage 6 cluster filter afterwards** — splitting *manufactures* small clusters. Fragments failing the filter are **dropped and recorded** (size, donor counts, criterion failed), never merged back into the parent — the parent was tier C precisely because it was heterogeneous, so returning them undoes the split.

**Stopping rule.** Loop back to the filter and re-vote until **tier C is empty** — a real termination condition, unlike "sweep resolutions until it looks right." The loop is not free: expect **tier D to grow as tier C shrinks**, because splitting converts internal heterogeneity into reference disagreement over smaller fragments (a real run went C: 4 → 1 while D: 4 → 11). Report both movements every round — whether another round is worth it is a judgement about that trade, and the analyst cannot make it if only the tier-C reduction is shown.

**Freeze & hierarchy.**

- Labels are recorded at **three levels** (`compartment` → `cell_type` → `cell_subtype`), each mapped to a **Cell Ontology (CL) term**, and each with its tier.
- Manual curation happens in one place: a review notebook that renders each cluster's evidence table and lets the analyst accept/override with a required free-text justification, written back into `uns["pipeline_log"]`.
- Subtype annotation recurses: subset a compartment, re-cluster within it (same latent), repeat 7a/7b with subtype-level references and markers — the filter and tiers apply at every level.

### Stage 9 — Validation of the annotation itself

Annotation gets its own QC, checked before sign-off:

- **Marker dot plot** of the final labels vs the curated dictionary — the single most persuasive figure; every label must "light up" its markers. Two rules: whenever clusters have been removed or split, **recompute the plot — never subset** a previously computed matrix, because z-scores and colour scales are cross-cluster statistics and subsetting leaves every dot scaled against clusters no longer shown; and keep the **marker-panel plot (specificity)** distinct from any **top-expressed-genes plot (abundance)** — the latter is dominated by housekeeping genes, and removing the worst offenders does not convert an abundance ranking into a specificity ranking.
- **Marker agreement rate**, reported with its matching rule stated. A rate is defined by its rule (exact label vs prefix match); comparing rates across figures requires the same rule, and if a rate is quoted alongside an earlier one, keep both columns so the comparison stays auditable — a rate can appear to fall purely because the rule tightened.
- **Composition sanity:** cell-type proportions per sample vs literature expectations; a cell type present in only one batch is a red flag for residual batch effect.
- **Label agreement metrics:** cluster purity of transferred labels, CellTypist-vs-scANVI Cohen's kappa (per compartment); a drop localizes where the annotation is soft.
- **Negative controls:** no cluster should be defined chiefly by mito/ribo/stress genes or by an ambient signature.
- **Held-out check (optional, recommended once per tissue):** annotate one sample with the pipeline blinded, compare to an expert's manual annotation.

### Stage 10 — Outputs

- `annotated.h5ad`: counts + lognorm layers, scVI/scANVI (and Harmony) latents, all Leiden resolutions, three-level labels + CL terms, **the evidence tier per label**, per-cell confidence, `unresolved`/doublet flags, full `pipeline_log`. The tier is the load-bearing column: a donor-level differential analysis run on tier-D labels supports a compartment-level conclusion and not a subtype one, and that constraint is invisible once the label is just a string.
- **Exclusion table** with per-criterion columns (size floor / donor breadth / dual-lineage / prior designation, plus the ambient-vs-doublet diagnosis) for every filtered cluster, and a **dropped-subcluster table** (name, size, donor counts, criterion failed) — cells leaving the analysis are accounted for somewhere a reader can check.
- Annotation report (HTML/PDF): per-stage gates, evidence tables, marker dot plot, UMAPs (by label / batch / confidence), composition tables, curation justifications.
- Versioned marker dictionary and config file — the reproducibility contract: same config + same inputs ⇒ same labels.

---

## 4. Key design decisions & rationale

**Why hybrid rather than reference-only?** References go stale, miss tissue-specific or disease states, and fail silently on populations absent from the atlas. Markers alone are slow and analyst-dependent. Cross-checking makes each method's failure mode visible to the other.

**Why cluster-level decisions from cell-level predictions?** Per-cell labels are noisy; aggregating to clusters denoises, while retaining per-cell probabilities preserves the ability to catch heterogeneous clusters (low purity ⇒ subcluster).

**Why scVI/scANVI as the integration+transfer backbone?** One latent space serves integration, clustering, *and* label transfer, so annotation and batch correction are consistent by construction; scArches extends it to pre-trained community references without retraining.

**Why per-sample QC before pooling?** MAD thresholds computed on pooled data conflate technical batch differences with quality differences and systematically over-filter the sample with the lowest depth.

**Why "unresolved" as a first-class label?** Forced assignment is the most common source of downstream false conclusions. An honest `unresolved` bucket (typically 1–5% of cells) is a feature.

**Why evidence tiers instead of flat labels?** A label is a claim, and claims vary in strength. Sorting clusters by how well the evidence converges — and recording the tier next to the label — lets a downstream differential analysis ask what a label is good for instead of assuming they are all equal. It also gives the annotation loop a real stopping rule (tier C empty) in place of the unfalsifiable "sweep resolutions until it looks right."

---

## 5. Tooling summary

| Stage | Primary tool | Alternative |
|---|---|---|
| Ambient RNA | CellBender | SoupX |
| QC | scanpy + MAD filtering | — |
| Doublets | scDblFinder + Scrublet + DoubletFinder (≥2/3 consensus) | — |
| Normalization | log1p (median) | scran pooling |
| Regress-out & scaling (PCA path) | sc.pp.regress_out + sc.pp.scale | — |
| Integration | scVI or Harmony (config-selected, can run both) | sysVI (cross-tech) |
| Clustering | Leiden (multi-res) | — |
| Auto-annotation | CellTypist + scANVI/scArches | popV, Azimuth (human) |
| Markers/DE | sc.tl.score_genes on user-defined, tissue-matched marker sets; rank_genes_groups | decoupler; PanglaoDB/CellMarker as scaffold only |
| Cluster filter & tier triage | `sc-cluster-annotation-loop` kernel (`dual_lineage_fraction`, `complexity_vs_doublet`, `rollup_vote`, `triage_tiers`, `refilter_subclusters`) | — |
| Ontology | Cell Ontology (CL) via obonet | — |
| Reporting | Jupyter → HTML | quarto |

---

## 6. Failure modes this design guards against

1. Ambient RNA read as real marker expression → Stage 1 + Stage 7b artifact-signature check.
2. Doublets annotated as novel intermediates → Stage 3 flags + Stage 8 dual-lineage rule.
3. Over-integration merging real cell types → Stage 5 bio-conservation gate.
4. Reference bias (atlas lacks your population) → marker evidence + `unresolved` path.
5. Resolution dogmatism → multi-resolution clustering + hierarchical labels.
6. Silent parameter drift across reruns → versioned config, marker dictionary, and pipeline log.
7. Species mix-ups (human model on mouse data) → Stage 0 species inference + config switch.
8. Split reference votes read as ambiguity (or as "the reference lacks this cell type") → parent-level rollup and a direct vocabulary query before any subclustering.
9. Subclustering manufacturing sub-floor fragments that get tiered as legitimate → the Stage 6 filter re-applies after every split, and dropped fragments are recorded, not merged back.
10. Cluster-mean lineage scores mistaking ambient baseline for doublets (or missing real mixing) → the per-cell two-threshold dual-lineage test, calibrated against the dataset's own median.
11. Tier definitions drifting from the code that computes them → tier criteria are emitted from the predicates themselves, and any iterative loop is run for two full rounds before being trusted — the untested path is the one the design rests on.

---

## 7. Suggested phasing for implementation

1. **Phase 1:** Stages 0–6 (preprocessing through clustering) — deterministic, easiest to gate automatically.
2. **Phase 2:** Stage 7a/7b engines + evidence-table generation.
3. **Phase 3:** Curation notebook, consensus rules, CL mapping.
4. **Phase 4:** Validation suite + report; then a blinded held-out benchmark on one well-understood tissue (PBMC is the canonical choice) before production use.

---

*References: [Single-cell best practices — Annotation](https://www.sc-best-practices.org/cellular_structure/annotation.html) · [popV consensus annotation, Nat Genet 2024](https://www.nature.com/articles/s41588-024-01993-3) · [scvi-tools docs](https://docs.scvi-tools.org/) · [OSCA — Cell type annotation](https://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html)*
