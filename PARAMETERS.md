# EFP model v2 — parameter provenance

All parameters in `sigmoid_model_v2.py`, with sources and links. Native units: E2 pmol/L,
FSH/LH IU/L, InhB pg/mL, P4 nmol/L, GnRH a.u. (A0=1), time min.

Confidence key: **solid** = directly measured/derived · **anchored** = set at operating level by
half-max convention (order-right, not independently measured) · **soft** = estimated · **placeholder** =
not yet sourced.

---

## 1. Clearances — `k = ln2 / t½`

| Hormone | t½ (min) | k (/min) | Source | Conf. |
|---|---|---|---|---|
| E2 | 25 | 0.0277 | IV estradiol terminal t½ 27.45±5.65 min, n=8 — White 1998, *Pharmacotherapy* [[DOI]](https://doi.org/10.1002/j.1875-9114.1998.tb03157.x) | solid |
| FSH | 240 | 0.0029 | canonical ~4 h; glycoform 343–757 min — Bousfield 2022, *JCEM* [[link]](https://academic.oup.com/jcem/article/107/10/e4058/6653068) | solid |
| LH | **90** (effective) | 0.0077 | slow glycoform terminal 80–196 min — Bousfield 2022 [[link]](https://academic.oup.com/jcem/article/107/10/e4058/6653068); also Wide 2009 [[PMC]](https://pmc.ncbi.nlm.nih.gov/articles/PMC2681272/). *Interpulse tail; fast phase ~20 min over-sharpened pulses* | moderate |
| P4 | 5 | 0.1386 | MCR 2100–2800 L/day → rapid — Little 1966, *J Clin Invest* [[link]](https://www.jci.org/articles/view/105405) | moderate |
| InhB | 120 | 0.0058 | inhibin B t½ 1.5 h & 3 h, **n=2 granulosa-tumor patients** — F&S 2000 [[DOI]](https://doi.org/10.1016/S0015-0282(00)00453-2) *(citation from prior research — verify)* | **shaky** |
| GnRH | 3 | 0.23 | GnRH t½ ~2–4 min (standard) — see pulsatility review [[PMC]](https://pmc.ncbi.nlm.nih.gov/articles/PMC4307809/) | solid |

---

## 2. Half-saturation constants `K` (grouped by ON/OFF in EFP)

| K | Value | Driver | Basis | Source | Conf. | EFP |
|---|---|---|---|---|---|---|
| K_FSH_E2 | 6.5 | FSH | EFP FSH median (half-max) | FSH ranges — Medscape [[link]](https://emedicine.medscape.com/article/2089048-overview) | anchored | **ON** |
| K_FSHinh | 6.5 | FSH | FSH median (clone of E2 loop) | same | anchored | **ON** |
| K_inh | 80 | InhB | EFP InhB median | PMC8081350 [[link]](https://pmc.ncbi.nlm.nih.gov/articles/PMC8081350/) *(assay-dependent; lit 10–150)* | soft | **ON** |
| K_G | 1.0 | GnRH | GnRH pulse **peak** (=A0); MM half-max at peak | convention/anchor (this project) | anchored | **ON** |
| K_amp | 315 | E2 | Rasgon amplitude-law midpoint (8.2−0.013·E2, half at ~315) | Rasgon 2003 (repo PDF) | soft | PARTIAL |
| K_sens | 734 | E2 | surge threshold 200 pg/mL × 1000/272.4 | Kauffman 2022, *Front Neurosci* [[link]](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.953252/full) | derived | OFF |
| K_freq | 500 | E2 | mid-follicular E2 (frequency-rise midpoint) | pulse data — Endotext [[link]](https://www.ncbi.nlm.nih.gov/books/NBK279070/) | soft | OFF |
| K_LH_P4 | 7 | LH | LH median (luteal param; V_P=0 in EFP) | PMC8042396 [[link]](https://pmc.ncbi.nlm.nih.gov/articles/PMC8042396/) | anchored | OFF |

---

## 3. Production amplitudes `V` — solved from steady state `production = k·X`

| V | Value | Solved from | Conf. |
|---|---|---|---|
| V_E | 7.764 | `k_E2·E2 / Hp(FSH;6.5,4)` = 0.0277·140 / 0.5 | solid |
| V_inh | 0.924 | `k_inh·InhB / Hp(FSH;6.5,4)` = 0.0058·80 / 0.5 | soft (inherits k_inh) |
| V_L | 2.17 | `k_LH·LH / ⟨MM⟩` = 0.0077·7.1 / 0.0252 | solid |
| V_F | 0.574 | `(k_FSH·FSH − b_FSH) / (⟨MM⟩·Hm(InhB;80,4))` = (0.0188−0.0116)/(0.0252·0.5) | solid |
| V_P (EFP) | 0 | anchored — no granulosa LH receptors in EFP (P4 receptor-absence exception) | design |
| V_P (MLP) | ~10.77 | `(k_P4·P4_mlp − b_P4)/Hp(LH_mlp;7,4)` — **luteal, out of EFP** | pending recal |

`⟨MM⟩ = 0.0252` = cycle-averaged GnRH coupling `⟨MM(GnRH;K_G=1)⟩` = `ln(1+A0/K_G)/(k_GnRH·T)`,
set by the generator (period T≈120 min at EFP), independent of the V's.

---

## 4. Basals `b`, GnRH generator, E2 sensitivity

### Basals — additive, GnRH-independent floors
| b | Value | = | Basis | Conf. |
|---|---|---|---|---|
| b_FSH | 0.0116 | `k_FSH · FSH_nadir(4)` | constitutive FSH (GnRH blockade: LH dies, FSH persists); number from late-foll nadir | moderate |
| b_P4 | 0.029 | `k_P4 · P4_efp(0.21)` | adrenal/constitutive EFP P4 (no LH-R yet) | solid |

### GnRH generator
| Param | Value | Basis | Source | Conf. |
|---|---|---|---|---|
| A0 | 1.0 | reference kick (convention, absorbed into K_G) | — | fixed |
| k_GnRH | 0.23 /min | t½ ~3 min | (§1) | solid |
| f_min | 0.0083 /min | 1 pulse/120 min (early-foll floor) | Endotext [[link]](https://www.ncbi.nlm.nih.gov/books/NBK279070/); pulsatility [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/2364566/) | anchored |
| f_max | 0.0167 /min | 1 pulse/60 min (late-foll ceiling) | same | anchored |
| K_freq | 500 | E2 half-max of freq rise | (§2) | soft |
| n_freq | 4 | default cooperativity | — | free |
| K_amp | 315 | E2 half-suppression of amplitude | Rasgon 2003 | soft |
| n_amp | 2 | gentle (≈ linear Rasgon law) | — | free |
| E2_ref | 140 | EFP E2 (amplitude normalized to A0 here) | (§5) | — |

### E2 pituitary priming (gain = 1 + A·Hp(E2); ≈off in EFP)
| Param | Value | Basis | Source | Conf. |
|---|---|---|---|---|
| K_sens | 734 | surge threshold (§2) | Kauffman 2022 [[link]](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.953252/full) | derived |
| n_sens | 8 | steep (switch-like surge) | — | free |
| A_LH | 4.0 | priming fold; LH>FSH but no number (Veldhuis JCI 1986) | **placeholder** | placeholder |
| A_FSH | 1.0 | smaller than A_LH | **placeholder** | placeholder |

---

## 5. Operating levels (anchoring targets & initial conditions)

| Hormone | EFP level | Source | Conf. |
|---|---|---|---|
| E2 | 140 pmol/L | early-foll (day-3); PMC8042396 foll range 114–332 (median 198 = mid-foll) [[link]](https://pmc.ncbi.nlm.nih.gov/articles/PMC8042396/) | solid |
| FSH | 6.5 IU/L | early-foll median — Medscape [[link]](https://emedicine.medscape.com/article/2089048-overview) | solid |
| LH | 7.1 IU/L | foll median 7.14 — PMC8042396 [[link]](https://pmc.ncbi.nlm.nih.gov/articles/PMC8042396/) | solid |
| P4 | 0.21 nmol/L | foll median 0.212 — PMC8042396 | solid |
| InhB | 80 pg/mL | PMC8081350 [[link]](https://pmc.ncbi.nlm.nih.gov/articles/PMC8081350/) *(assay-dependent)* | soft |
| FSH nadir | 4 IU/L | late-foll (pins b_FSH) — Medscape/Mayo ranges | moderate |

**Exponents n** (all free sensitivity-analysis knobs): default n=4; n_sens=8; n_inh=4; n_amp=2.

---

### Full source list
- White 1998 *Pharmacotherapy* (E2 t½): https://doi.org/10.1002/j.1875-9114.1998.tb03157.x
- Bousfield 2022 *JCEM* (FSH/LH glycoform t½): https://academic.oup.com/jcem/article/107/10/e4058/6653068
- Wide 2009 *JCEM* (LH/FSH sialylation t½): https://pmc.ncbi.nlm.nih.gov/articles/PMC2681272/
- Little 1966 *J Clin Invest* (P4 MCR): https://www.jci.org/articles/view/105405
- Fertility & Sterility 2000 (inhibin B t½, n=2): https://doi.org/10.1016/S0015-0282(00)00453-2
- PMC8042396 (E2/LH/P4 cycle medians): https://pmc.ncbi.nlm.nih.gov/articles/PMC8042396/
- PMC8081350 (inhibin B levels): https://pmc.ncbi.nlm.nih.gov/articles/PMC8081350/
- Kauffman 2022 *Front Neurosci* (E2 surge threshold): https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.953252/full
- Endotext NBK279070 (GnRH/LH pulse frequencies): https://www.ncbi.nlm.nih.gov/books/NBK279070/
- Pulsatility review PMC4307809: https://pmc.ncbi.nlm.nih.gov/articles/PMC4307809/
- LH pulsatility across cycle, PubMed 2364566: https://pubmed.ncbi.nlm.nih.gov/2364566/
- Medscape FSH reference ranges: https://emedicine.medscape.com/article/2089048-overview
- Rasgon et al. 2003 (amplitude law; source PDF in repo)
