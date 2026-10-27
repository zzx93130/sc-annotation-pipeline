# Lung IPF atlas — validated marker panels & annotation record

Run record from a full human-lung IPF-vs-Control annotation (8 public scRNA-seq
studies, benchmark_v3_no_GSE161685). 840,419 cells at level 1 → **745,281 cells
retained across 5 compartments and 77 reviewed subtypes** after per-compartment
cleaning. Every panel below is the *curated starting panel* actually scored; the
"validated discriminators" are the genes that separated each subtype in this
cohort (read from the frozen per-cluster evidence, cross-study, ≥8 datasets).
Markers keep their source — do not treat this as a from-memory list; it is a
reproduced result to reuse as a scaffold, re-graded inside each compartment.

Cohort: IPF 404,220 / Control 341,061 cells · 8 datasets · 239 samples · 186 donors.

---

## Level-1 compartments (major lineage)
Score the union panel per compartment at Leiden ~0.3–0.5; assign the major lineage,
then subset and re-derive HVGs/markers inside each. Compartments recovered here:
Myeloid (343,345), Epithelial (218,482), Lymphoid (90,310 → T/B/NK-ILC), Endothelial
(60,017), Fibroblast/stromal (33,127). Lymphatic EC (LEC) is its own L1 compartment,
separate from blood EC.

---

## Myeloid (14 subtypes; L2 panel)
- **Alveolar macrophage** — FABP4, MARCO, MRC1
- **Alveolar macrophage (metallothionein-high)** — MT1E, MT1X
- **Interstitial macrophage** — FOLR2, LYVE1
- **Classical monocyte** — FCN1, S100A8, S100A12, CD14, S100A9
- **Non-classical monocyte** — FCGR3A, LILRB2, MS4A7, FCN1
- **cDC2** — CD1C, FCER1A, CLEC10A, CD1E
- **mDC** — LAMP3, CCR7, FSCN1
- **pDC** — IRF7, TCF4, PLD4, SERPINF1
- **Proliferating myeloid** — TYMS, MKI67, STMN1
- **Profibrotic macrophage (SPP1+)** — SPP1, MATK, SDC2, LILRB4, CHI3L1, MMP9
- **Lysosomal / lipid-handling macrophage** — LGMN, RNASE1
- **cDC1 (CLEC9A+XCR1+)** — CLEC9A, XCR1, IRF8

Final reviewed labels: AM (alveolar macrophage), CCL2+ IM, Classical monocyte, HSP-high macrophage, IFN-response macrophage, IM (interstitial macrophage), JARID2+ macrophage, LSAMP+ AM, Metallothionein-high macrophage, Non-classical monocyte, SPP1+MMP9+ macrophage, cDC1, cDC2, mDC.
IPF-relevant: **SPP1+MMP9+ (profibrotic) macrophage** (SPP1/MMP9/CHI3L1/LILRB4) is
the single most abundance-shifted subtype (strongly expanded in IPF); alveolar
macrophage (FABP4/MARCO) and classical monocyte are depleted.

## Epithelial (15 subtypes; L2 panel + required-gene gates)
- **AT0 / pre-terminal bronchiole secretory** — SFTPB, SCGB3A2
- **AT1** — AGER, HOPX, RTKN2, PDPN, CLIC5
- **AT2** — SFTPC, SFTPD, SFTPA1, NAPSA
- **AT2 (activated/transitional)** — CDKN3, BIRC5, KRT8
- **Aberrant basaloid** — COL1A1, KRT17, CDH2
- **Basal (KRT5+/TP63+)** — KRT5, KRT14, TP63
- **Airway secretory (club/goblet, MUC5B+)** — MUC5B, TFF3, MUC5AC
- **Multiciliated (ciliated)** — FOXJ1, TPPP3
- **PNEC** — ASCL1, CALCA, CHGA
- **ionocyte** — CFTR, FOXI1, ASCL3
- **IPF aberrant secretory (MMP7+)** — MMP7, CTSE, CEACAM6, PDZK1IP1, SAA1

Required-gene gates (a subtype call is vetoed unless these are detected):
{
"ionocyte": [
"FOXI1",
"ASCL3"
],
"PNEC": [
"ASCL1",
"CHGA"
],
"IPF aberrant secretory (MMP7+)": [
"MMP7"
],
"AT1": [
"PDPN",
"RTKN2",
"CLIC5"
]
}

Validated discriminators (this cohort):
  - IFN-response AT2: IFN-response AT2 (ISG15/IFI6/MX1) — state
  - AT2: AT2 (SFTPA1/NAPSA/SFTPD)
  - Transitional AT1: AT1/AT2 transitional (HOPX+SFTPB 98) — gate held
  - AT2: AT2 (PGC/LAMP3/ABCA3)
  - AT2: AT2 (SFTPA/FABP5)
  - CCL20+ AT2: AT2 inflammatory (CXCL2/CCL20/CSF3/ICAM1) — NF-kB substate
  - AT1: AT1 (AGER/CAV1/CAV2/SPOCK2) — clean
  - AT0: AT0/pre-terminal secretory (SCGB3A2/MGP)
  - Airway secretory club: Airway secretory (BPIFB1/MUC5B/TSPAN8)
  - Proliferating epi: AT2 activated/transitional, cycling (STMN1/TYMS/PCLAF)
  - Ciliated: Multiciliated (RSPH1/C20orf85)
  - Ciliated: Multiciliated (SNTN/CCDC78)
  - SOX5+ Ciliated: Multiciliated (AGBL4/CFAP54 axonemal)
  - Aberrant basaloid: Aberrant basaloid (COL1A1/KRT17/LAMB3/FN1/MMP7) — IPF hallmark
  - CCNO+ Ciliated: Multiciliated deuterosomal (CDC20B/CCNO/FOXJ1)
  - Airway secretory club: Airway secretory (MDK/TSPAN8/AQP5/S100P)
  - Basal: Basal KRT5+/TP63+ (KRT17/S100A2/KRT5)
  - SOX5+ AT0: AT0/pre-terminal (PDE4D/SOX5 low-complexity)
  - PNEC: PNEC / ionocyte (mixed rare; TFF3/HEPACAM2/ASCL3)

IPF hallmark: **Aberrant basaloid** (COL1A1/KRT17/CDH2/LAMB3/FN1/MMP7) — near-absent
in control, the largest fold-change in the whole atlas. AT2 and AT1 are depleted
(alveolar epithelial loss); basal and PNEC expand (bronchiolization). An AT1/AT2
transitional population sits at the boundary and is held by the AT1 required gate.

## Lymphoid — T (14 subtypes; curated Pan-T + lineage + state panel)
- **CD4 T** — CD4, CD40LG, TRAT1
- **CD8 T** — CD8A, CD8B
- **gd T** — TRDC, TRGC1, TRGC2, TRDV2, TRGV9
- **MAIT** — SLC4A10, TRAV1-2, KLRB1
- **Treg** — FOXP3, IL2RA, IKZF2, CTLA4
- **Tfh** — CXCL13, PDCD1, ICOS, TOX2, CD200
- **Naive T** — CCR7, SELL, LEF1, TCF7, FOXP1, SATB1, PDE3B
- **Proliferating T** — MKI67, TOP2A, TYMS, STMN1
States: Stress / heat-shock (HSPA1A, HSPA1B, DNAJB1, HSPH1)

Final reviewed labels: Activated T, CD4 Th, Cytotoxic CD8 T, GZMK+ CD8 T, IFN-response T, MAIT, Naive T, Proliferating T, RORA+ CD4 Th17, TRBV9+ CD8 T, Tfh, Treg, Trm, high-Meta T.
Notes: γδ markers (TRDC/TRGV9) present but sub-floor in the Cytotoxic CD8 cluster →
resolved to CD8 (not γδ). "high-Meta T" = metallothionein-high substate (MT1X/MT1E/MT2A).

## Lymphoid — B (9 subtypes)
- **Naive B cell** — IGHD, TCL1A
- **GC B cell** — RGS13, LMO2
- **Memory B cell** — CD27, TNFRSF13B
- **Plasmablast** — XBP1, PRDM1
- **Plasma cell** — MZB1, JCHAIN

Final reviewed labels: Atypical memory B, GC B, IFN-response B, IGLL1+ Plasma B, Memory B, Naive B, Plasma B, Proliferating B, high-meta B.
Notes: plasmablast vs plasma is gated by proliferation (MKI67/TOP2A/STMN1) — only the
cycling cluster is plasmablast. Plasma cells carry ambient myeloid transcript (LYZ 20–33%)
— NOT a doublet when MZB1-high/MS4A1-negative. B+T/NK and B+myeloid doublets show
exclusive foreign programs (CD3E/NKG7 or LYZ/S100A8) and were removed.

## Lymphoid — NK/ILC (7 subtypes)
- **CD56bright NK** — NCAM1, SELL, XCL1, XCL2, GZMK, KLRC1, IL7R
- **CD56dim NK** — FCGR3A, GNLY, GZMB, PRF1, FGFBP2, CX3CR1, SPON2
- **Adaptive NK** — KLRC2, GZMH, IL32, CD52
- **ILC** — IL7R, KIT, KLRB1, AREG, TNFRSF18
- **ILC2** — GATA3, PTGDR2, IL1RL1, HPGDS
- **ILC3** — RORC, IL23R, AHR
States: Proliferating, IFN-response, high-MT, high-HSP

Final reviewed labels: CD56bright NK, CD56dim NK, HLA-DR+ NK, HSP-high NK, IFN-response NK, IFNG+ NK, ILC.
Validated discriminators (this cohort):
  - IFNG+ NK: IFNG 40 + chemokine CCL3 75/CCL4 96; activated cytokine NK (note c1 IFNG 46 comparable)
  - HSP-high NK: HSP-high: HSPA1A/HSPA1B/DNAJB1 stress program
  - ILC: ILC: KIT 13/GATA3 23/IL7R 55/AREG 82, cytotoxic LOW (FCGR3A 12); 1341/1422 r1-ILC cells; splits cleanly at Leiden1.0 c11 (KIT 38)
  - CD56dim NK: resting CD56dim: FCGR3A 69/GNLY 93/GZMB 94
  - CD56bright NK: CD56bright: FCGR3A LOW 21/GZMK 51/XCL1 61/KLRC1 55/SELL
  - CD56dim NK: terminal cytotoxic CD56dim: FCGR3A 87/SPON2 87/FGFBP2 88/GZMB 98
  - HLA-DR+ NK: HLA-DR+ NK: NKG7 99/KLRD1 89/GNLY 96 + CD74 91/HLA-DRB1 67; CD3E 61 ambient only (CD3D 12/CD3G 9) NOT T
  - IFN-response NK: IFN-response: ISG15/MX1/IFI6/IFIT3
CRITICAL: AUCell over-calls "ILC" from AREG/IL7R, which NK also express — the real
ILC discriminator is **KIT/RORC** (~0 in NK). Only one genuine ILC population here.
HLA-DR+ NK has CD3E ambient (61%) but NO full CD3 complex (CD3D 12/CD3G 9) — it is NK,
not T contamination.

## Fibroblast / stromal (10 subtypes)
- **Alveolar fibroblast** — NPNT, TCF21, GPC3, SCN7A, LIMCH1
- **Adventitial fibroblast** — MFAP5, PI16, CD34, SFRP2, SERPINF1
- **Myofibroblast (CTHRC1+)** — CTHRC1, POSTN, COL1A1, COL3A1, COL5A2, TNC
- **Inflammatory fibroblast** — CCL2, CXCL12, CXCL14, PDGFRA
- **Immune-recruiting fibroblast** — CXCL13, CCL19, CCL21
- **HAS1-high fibroblast** — HAS1, TWIST1, PLIN2
- **Proliferating fibroblast** — STMN1, TYMS, MKI67, TOP2A
- **Smooth muscle / pericyte** — ACTA2, TAGLN, MYH11, PDGFRB, NOTCH3
- **Mesothelial** — MSLN, UPK3B, WT1, CALB2
States: IFN-response, high-HSP, high-MT

Final reviewed labels: Alveolar fibroblast, CD55+ Adventitial fibroblast, DPT+ Adventitial fibroblast, HAS1+ Subpleural fibroblast, Mesothelial, Myofibroblast, Proliferating fibroblast, SFRP1+ Adventitial fibroblast, SMC-like, Stressed fibroblast.
Validated discriminators (this cohort):
  - Alveolar fibroblast: NPNT 41/TCF21 73/GPC3 77/LIMCH1 74/ITGA8 55 — alveolar
  - SMC-like: TAGLN 80/ACTA2 56/ASPN 72 (MYH11 30) — SMC-like (contractile, weak MYH11)
  - Mesothelial: KRT18 91/CALB2 79/UPK3B 76/MSLN 71 — mesothelial
  - Myofibroblast: CTHRC1 82/POSTN 77/FAP 57/COL1A1 98/ACTA2 61 — IPF pathological myofibroblast
  - Proliferating fibroblast: STMN1 91/TYMS 76/TOP2A 49/MKI67 46 — proliferating
  - DPT+ Adventitial fibroblast: DPT 78/LEPR 67/ASPN 63/LTBP2 71/TCF21 70 — DPT+ adventitial
  - CD55+ Adventitial fibroblast: CD55 80/MFAP5 97/PI16 69/CD34 74/SCARA5 74 — CD55+ adventitial
  - SFRP1+ Adventitial fibroblast: SFRP1 50/C3 96/MFAP5 68/WISP2 67/PLA2G2A 75 — SFRP1+ adventitial
  - Stressed fibroblast: HSPA1A 89/DNAJB1 89/FOS 90/JUN 93/ATF3 80 — stress/AP-1 state
  - HAS1+ Subpleural fibroblast: HAS1 62 (z+1.6)/IL6 72/PTX3 73/CXCL8 62/CXCL2 81 — inflammatory-activated program; 'subpleural' is the reviewer's anatomical assignment (no single marker in-panel discriminates it; SFRP2 62% is at cross-cluster average, NOT enriched)
IPF hallmarks: **Myofibroblast (CTHRC1+)** (CTHRC1/POSTN/FAP/COL1A1) and **HAS1+
subpleural** (HAS1/IL6/PTX3/CXCL8) — both strongly expanded. Three adventitial
subtypes split on CD55 / SFRP1(C3) / DPT+LEPR. "Stressed fibroblast" = HSP/AP-1
immediate-early state (HSPA1A/DNAJB1/FOS/JUN), likely dissociation-induced, not a
distinct identity. SMC-like has weak MYH11 (contractile, not textbook smooth muscle).

## Endothelial (8 subtypes)
- **EC capillary — gCap (general)** — APLNR, EDNRB, RGCC, IL7R, FCN3
- **EC capillary — aerocyte (aCap)** — CA4, HPGD, TMEM100, SOSTDC1, EDA
- **EC arterial** — GJA5, DKK2, IGFBP3, BMX, SERPINE2
- **EC venous (pulmonary)** — ACKR1, SELP, CPE, PLVAP, HDAC9
- **EC venous (systemic/peribronchial)** — COL15A1, VWA1, PLVAP, SPRY1
- **Tip cell** — ESM1, APLN, ANGPT2, DLL4, KDR
- **EC lymphatic (LEC)** — CCL21, TFF3, LYVE1, PROX1, MMRN1
- **EC IFN-response** — ISG20, CXCL10, CXCL11, IFIT3, ISG15
States: Proliferating, high-HSP, high-MT

Final reviewed labels: ACKR1+ venous EC, Arterial EC, IFN-response EC, Proliferating EC, Pulmonary venous EC, Systemic venous EC, aCap, gCap.
Validated discriminators (this cohort):
  - Systemic venous EC: COL15A1 52/CD34 70/PLVAP 73/SPRY1 72/VWA1 57 — systemic/peribronchial venous
  - Systemic venous EC: COL15A1 47/ADAMTS9 71/ACKR1 83/POSTN 68 — systemic venous (ADAMTS9+; ER 5.13 study-skewed)
  - gCap: FCN3 94/RGCC 86/IL7R 67/TMEM100 72 — general capillary (gCap)
  - Pulmonary venous EC: ACKR1 97/CPE 63/CLU 94/C7 67/SELP 50 — pulmonary venous
  - ACKR1+ venous EC: ACKR1 97/PLVAP 77/SPRY1 79/POSTN 79/CCL14 69 — venous (ACKR1-high)
  - Pulmonary venous EC: CPE 85/EFEMP1 66/PTGIS 62/CLU 99 — pulmonary venous (PTGIS+)
  - Arterial EC: GJA5 53/DKK2 58/HEY1 61/EFNB2 59 — arterial
  - aCap: CA4 67/HPGD 81/EDNRB 70/IL1RL1 65 — aerocyte (aCap)
  - Arterial EC: GJA5 50/DKK2 56/BMX 43/MECOM 54/FN1 63 — arterial (BMX+/MECOM+ subtype)
  - Proliferating EC: STMN1 92/TYMS 70/TOP2A 47/MKI67 45 — proliferating
  - IFN-response EC: GBP1 77/IFIT3 67/CXCL10 58/CXCL11 54/ISG20 62 — IFN-response state
Full lung EC taxonomy: gCap (FCN3/RGCC/IL7R) vs aCap/aerocyte (CA4/HPGD/EDNRB),
arterial (GJA5/DKK2/HEY1), and three venous programs (ACKR1-high, pulmonary CPE/CLU,
systemic/peribronchial COL15A1+). No LEC in the blood-EC compartment (PROX1 <1%) —
lymphatic EC are a separate L1 compartment. ACKR1+/systemic venous expand in IPF.

---

## IPF differential-abundance signature (per-donor proportion, Mann-Whitney, BH-FDR)
46/77 subtypes significant (FDR<0.05). Canonical fibrotic signature, all compartments:
**Expanded** — Aberrant basaloid, CTHRC1+ myofibroblast, HAS1+ subpleural fibroblast,
SPP1+MMP9+ macrophage, Basal/PNEC epithelial, Tfh/Treg/GC B/Memory B/IFN-response T
(tertiary lymphoid structures), ACKR1+/systemic venous EC.
**Depleted** — AT2, AT1, alveolar macrophage, classical monocyte, CCL20+ AT2.

## Method conventions settled in this run (lung IPF)
- **Integration:** dataset-aware HVG (per-DatasetID seurat_v3 consensus, restricted to
  datasets ≥1000 cells when a compartment is small; raise the floor if the per-batch
  loess is ill-conditioned) + SampleID Harmony (theta=2). Standard theta gave the best
  cross-study mixing (median enrichment ratio ~1.0–1.3) while preserving subtypes.
- **Iterative cleaning:** remove flagged doublet/contamination clusters, then re-run the
  full build unchanged; repeat until clean. Study-confined-but-real islands (e.g.
  GSE136831 memory B) are kept + flagged, not removed — removal only regenerates them.
- **AUCell is scaffold, not verdict:** it mis-calls where a channel marker is shared
  (ILC from AREG, LEC from ACKR1, plasma-as-memory). Always confirm calls on exact
  per-cluster detection% + cluster-vs-rest DE, read programmatically (never eyeball a
  wide detection table — column read-shifts are the recurring error).
- **var_names:** the level-1 counts h5ad stores Ensembl IDs; re-index to
  `gene_symbol_unique` before scoring markers.
