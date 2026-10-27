# Dataset and batch effects in a multi-study atlas

A cluster that is almost entirely one dataset is a warning, not a verdict. This
document is the diagnostic that tells you **which** kind of problem you have,
because on one 9-study lung atlas there were three distinct causes and each
needed a different answer. Applying one global correction would have fixed two
and destroyed a real population in the third.

## First: is it enrichment, or just abundance?

Score every cluster's dominant dataset share against **that dataset's share of
the same lineage subset**, not against its share of the atlas. A dataset
supplying 57% of a compartment will produce clusters at 80–88% by chance alone —
about 1.5x its own base rate, which is nothing. The same 85% share is a 14.7x
enrichment for a dataset supplying 6.7%.

This one ratio separated the real cases from the artefacts of abundance in every
lineage tested. Compute it before anything else.

## Cause 1 — intron-inclusive quantification

**Signature.** One dataset puts a much larger fraction of its counts into
long-intron genes, and a correspondingly *smaller* fraction into intron-free
genes. Measured within a single compartment (so composition cannot explain it):
26.6% of counts in top-quartile-intron genes against ~10% for every other study,
and simultaneously the lowest share in intron-free protein-coding genes.

**Cause.** The alignment was run with intron inclusion (`cellranger
--include-introns`) while the others were not, so intronic reads inflate genes
with large introns and dilute single-exon ones.

**Verify with the two-sided control.** Intron-free genes cannot gain counts from
intron inclusion, so their share must be *depressed* in the affected dataset.
Restrict that control to intron-free **protein-coding** genes excluding
mitochondrial and ribosomal ones — the unrestricted intron-free set is dominated
by a handful of extremely abundant single-exon transcripts whose share is driven
by biology, and the control fails to corroborate.

**Fix: dataset as an explicit covariate, at BOTH steps.** Use the study id as the
batch key for HVG selection *and* for batch correction. Selecting variable genes
per sample leaves the intron-inclusive regime free to nominate its inflated
genes, so correcting only the embedding is not enough.

**Result where it applies.** Myeloid 5.79% → 0.00% of cells in a dominated
cluster with all 11 subtypes retained and labelled fraction rising; T/NK 19.87%
→ 0.01% with all four at-risk CD8/CD4 labels kept and none resting on one
dataset. Both needed a resolution change alongside (below).

**Cost, which must be reported.** The correction removes the axis it corrects
along. On T/NK an ILC population merged into NK and did **not** re-separate at
any resolution (35% purity at 1.2, 1.6 and 2.0 alike) — its cells remain
distinct within the merged cluster (`KIT` 20.9% vs 1.3%, `IL7R` 69.4% vs 17.7%,
cytotoxic markers *lower*), so this is a real loss, not a hairline artefact being
corrected. Check the merged cluster's internals before accepting a merge.

## Cause 2 — tissue handling and sorting protocol

**Signature.** The dominant dataset's intron load is ordinary, and the clusters
carry a real but *induced* program — inflammatory cytokines on an otherwise
normal lineage core.

**Diagnosis.** Score dissociation-stress genes (`FOS`, `FOSB`, `JUN`, `JUNB`,
`EGR1`, `ATF3`, `HSPA1A`, `HSPA1B`, `DNAJB1`, `HSPB1`, `SOCS3`, `ZFP36`, `IER2`,
`KLF6`) per dataset. On the lung atlas the outlier dataset ranked first at 86.0%
mean detection and its cytokine level tracked stress at Spearman 0.79 across
studies. Its published protocol explained it: explant tissue held *on wet ice
overnight*, then minced, digested and FACS-sorted into stromal fractions — both
established inducers of immediate-early and cytokine programs. The comparison
dataset was whole-lung dissociate with no sorting step.

**Exclude the disease confound, and check its direction.** A study contributing a
different disease could produce the same clusters. Within that dataset the
inflammatory fraction was *higher in controls* (48.5%) than in cases (23.7%) —
the opposite of disease biology, which is what settles it.

**Do NOT apply the covariate here.** Correcting on the study id would merge a
genuine inflammatory-fibroblast phenotype into the adventitial/alveolar bulk. The
phenotype is real; what is not comparable is its **abundance** across datasets.
Instead: keep the sample-level embedding, name the clusters for what they are
(`Inflammatory / stress-activated fibroblast`), and mark them unusable for
cross-dataset abundance comparison.

## Cause 3 — dominance that survives correction

A population can stay dataset-confined after the covariate is applied — on the B
lineage a 2,344-cell cluster at 92% one study simply fragmented into three
clusters at 87/90/96% at finer resolution, the same cells each time. Record it as
unexplained rather than reaching for either of the causes above.

## When the covariate makes things WORSE

Test it, do not assume it. Compare **at matched resolution**: on epithelium the
covariate embedding carried more dominance at every resolution tested
(0.17→0.89, 1.18→4.08, 0.16→13.67). The reason was Cause-0 abundance — that
compartment was 57% one study, whose clusters were barely enriched — so the
correction spent itself on the wrong variable while licensing a finer clustering
that fragmented the abundant dataset into more small, easily-dominated clusters.

## Resolution interacts with all of this

Finer clustering mechanically produces more small clusters and small clusters are
more easily dominated, so a dominance figure is meaningless without its
resolution. Two consequences:

- Compare configurations **at matched resolution**, or the comparison is
  confounded. A covariate run that chose a finer resolution than its baseline
  looked worse (5.79% → 7.01%) while being better at matched resolution (5.79% →
  4.09%).
- Add a **dominance term to the resolution selector**, ranked *after* subtype
  recovery and *before* cluster count. With dataset corrected, every resolution
  from 0.8 up recovered the same subtypes, so cluster count alone climbed past
  the resolution that recovered them with zero dominated clusters and the best
  margins.

## Gene-length-aware DE: measure before correcting

CQN and length-regression exist to remove a coupling between transcript length
and counts that arises in **full-length bulk** protocols. UMI 3' counting
registers one molecule per transcript regardless of length, so the coupling may
not exist — and on this atlas it was *negative* in all six lineages (Spearman
−0.047 to −0.120), the opposite sign to bulk length bias, explaining little
between-gene variance.

So: **measure the effect on your own data first and report the measurement
beside any correction.** Where a length-aware version is wanted, treat it as a
sensitivity analysis, not a replacement for the primary DE.

Three implementation traps, all of which produced wrong numbers before being
fixed:

- **CQN needs count-scale input.** Feeding it summed log-normalised values is
  silently wrong; pseudobulk from raw counts and assert integrality.
- **A plain LOWESS on fold-changes fits the extremes.** Tiny immunoglobulin
  segments with pathological fold-changes at the sparse length extremes drove an
  absurd "correction". Use a robust fit restricted to genes where the fold-change
  is meaningful, and flag the excluded ones.
- **Centre the trend within each cluster.** An uncentred fit absorbs the
  cluster's overall effect-size level into its own intercept — a quantity
  comparable in magnitude to the genuine length slope, and real biology. Not
  centring caused 31% gene turnover; centring cut it to 4.6%.

If the covariate of interest is intron length specifically (the Cause-1 case),
note that intron length is only weakly coupled to exon length (log-log r ≈ 0.57),
so regressing exon length cannot remove it. And a within-dataset replication test
is the honest check on whether the association is artefact: on this atlas, pooled
DE genes that replicated across *both* quantification regimes had **longer**
introns than those replicating in neither, so intron length carries biology as
well as regime and regressing it out removes real signal.

## What to hand the user

A per-cluster table with the dominant dataset, its share, that dataset's share of
the lineage, and the enrichment ratio; the per-dataset diagnostic scores
(intron load, stress score) that assigned the cause; and a per-lineage decision
naming which configuration was adopted and why. The decision differs by lineage —
that is the finding, not a failure to converge.
