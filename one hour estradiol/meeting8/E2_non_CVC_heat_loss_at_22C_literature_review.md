# Estradiol-driven heat loss at 22 C, excluding cutaneous vascular conductance

## Question

At an ambient temperature of 22 C, is there evidence that estradiol (E2) increases genuine heat loss through a mechanism other than altered cutaneous vascular conductance (CVC), and is there a defensible coefficient for such an E2 term?

## Bottom line

**There is no published, positive, mechanism-specific non-CVC E2 heat-loss coefficient for a mouse at 22 C.** The coefficient that is strictly supported for an explicitly non-CVC E2 loss arm is therefore:

`beta_E2_nonCVC = 0`

This is not a claim that every possible non-CVC effect is biologically zero. It means that none has been shown at 22 C with a heat-flux measurement that can identify an E2 coefficient.

The closest quantitative fallback is an old rat gradient-layer calorimetry experiment. At 20 C, after both groups were completely shaved, physiological chronic E2 replacement increased measured whole-body **dry heat loss** from `7.91 +/- 0.34` to `8.87 +/- 0.08 W/kg`, a difference of:

`Delta q_dry(20 C) = +0.96 W/kg` (E2 replacement minus OVX control)

This is genuine measured heat leaving the animal, not a temperature proxy and not a reduction in thermogenesis. A follow-up thyroid-blockade experiment supported the interpretation that the extra loss can be primary and thermogenesis compensatory during severe cold. However, the mechanism remained unresolved, and the authors specifically considered altered peripheral vascular control. Therefore this value can be used only as an **empirical E2-associated whole-body dry-loss prior**, not as a proven non-CVC mechanism.

## Mechanism-by-mechanism findings

| Candidate route | Evidence near 22 C | Quantitative result | Can it supply a non-CVC E2 coefficient? |
|---|---|---:|---|
| Evaporative loss (saliva spreading, grooming, panting, respiratory and cutaneous water loss) | Direct rat measurements at 15 and 25 C found no E2 or estradiol-benzoate effect on evaporative heat loss. E2 increased evaporation only during severe heat exposure at 38 C, when saliva spreading was recruited. | Hormone main effect at 15/25 C: `P > 0.05`. Pooled evaporation was `0.61 +/- 0.50 W/kg` at 15 C and `1.24 +/- 0.16 W/kg` at 25 C in the large cage, but these are total evaporation values, not E2 effects. | **No. Set the E2-specific evaporative coefficient to zero at 22 C.** |
| Thermoregulatory behavior | In rats, low estrogen increased heat-escape behavior only during a hot-chamber test (32-40 C). E2 instead facilitated heat acquisition/tail-hiding behavior during cold exposure. In adult mice, gonadectomy did not alter thermal preference. | No E2-dependent heat-loss flux at 22 C was reported. | **No.** The available direction is neutral or heat-conserving at cool temperatures. |
| Fur/pelage insulation | Chronic E2 visibly retarded fur growth in rats. Complete shaving did not eliminate the dry-loss difference. | At 20 C after shaving: E2 `8.87 +/- 0.08 W/kg`; control `7.91 +/- 0.34 W/kg`. | **No separate fur coefficient.** Fur was not sufficient to explain the effect. |
| Body mass and surface-area-to-mass ratio | E2 reduced body mass, which can increase heat loss per kg. The original study calculated that mass could predict an 8% increase in experiment 1 and 4% in experiment 2, versus observed 14-17% increases in heat production; adjustment for mass left the hormone effect significant. | Mediated, chronic, and animal-specific; not an acute E2 flux. | **Not as a direct E2 term.** Model measured body mass or geometry separately. |
| Basal ventilation / respiratory evaporation | E2 replacement did not affect basal ventilation in normoxic wild-type mice in a direct ventilatory study. The thermoregulation literature does not report an E2-dependent respiratory heat flux at 22 C. | No usable flux difference. | **No.** |
| Whole-body dry heat loss | Direct gradient-layer calorimetry in dark-phase OVX rats detected greater dry heat loss with chronic E2 across ambient temperatures; the exact post-shaving result at 20 C is given above. | `+0.96 W/kg` at 20 C after shaving; plasma E2 `49.4 +/- 5.4 pg/mL` in the experiment. | **Yes as an empirical phenotype-level prior, but not as a proven non-CVC pathway.** |
| Whole-body thermal conductance in newer work | At 25-34 C, a later rat study found no pointwise difference in whole-body thermal conductance between OVX and OVX+E2 groups in either light or dark phase. | No group difference at each tested temperature. | This null result prevents treating the old `0.96 W/kg` as a universal coefficient. |

## What the exact 22 C mouse experiment does and does not show

Krajewski-Hall et al. exposed OVX CD-1 mice with or without chronic E2 to a randomized ambient-temperature series that included exactly 22 C. At 22 C, E2 lowered core temperature but did not change mean tail temperature, mean heat-loss index, or the within-hour HLI range. Their home-cage experiment also found no E2 effect on tail temperature, HLI, or activity. The study did not measure oxygen consumption, direct calorimetry, evaporative water loss, or respiratory heat loss.

Therefore the lower mouse core temperature at 22 C cannot be assigned to genuine heat loss rather than altered heat production. It supplies **no non-CVC heat-loss coefficient**.

## The direct-calorimetry coefficient

### Measured quantity

Laudenslager et al. used a gradient-layer calorimeter, not tail temperature, to measure whole-body dry heat loss. Rats were OVX, chronically replaced with E2 through a SILASTIC capsule or given an empty capsule, and tested during the dark phase. The main E2 effect on dry heat loss was significant across ambient temperatures (`P < 0.005`).

The cleanest exact comparison is the post-shaving test at 20 C:

- E2: `8.87 +/- 0.08 W/kg`
- OVX control: `7.91 +/- 0.34 W/kg`
- Difference: `0.96 W/kg`
- Approximate SE of the difference: `sqrt(0.08^2 + 0.34^2) = 0.35 W/kg`
- Plasma E2: `49.4 +/- 5.4 pg/mL`; control was below assay detection

This result rules out pelage as the full explanation but does not rule out altered skin blood flow elsewhere on the body.

### Directionality

In a follow-up study at -5 C, the same group suppressed thermogenic capacity using propylthiouracil or surgical thyroidectomy. Dry heat loss remained high in E2-treated rats while heat production fell below heat loss and core temperature declined. This supports the interpretation that E2-associated loss can be primary and thermogenesis compensatory. It does not establish the route at 22 C. The follow-up paper concluded that the mechanism was unknown and discussed peripheral vasculature, body mass, and altered hypothalamic control.

## Model-ready coefficients

### 1. Strict, mechanism-specific recommendation

If the model term is explicitly labeled **non-CVC E2 heat loss at 22 C**, use:

`E2_NONCVC_LOSS_GAIN = 0`

This is the only coefficient that does not overstate the literature.

### 2. Empirical whole-body dry-loss sensitivity analysis

If a phenotype-level E2 loss term is required despite the unresolved mechanism, use the calorimetry result as a rat-derived prior.

The measured 20 C coefficient is:

`q_E2,dry = 0.96 W/kg * E2_state`

where `E2_state = 1` represents chronic physiological replacement and `0` represents OVX control. Do not interpret `E2_state` as pg/mL; the experiment was not a dose-response study.

A simple dry-gradient extrapolation from 20 to 22 C, using a representative core temperature of 37.5 C, gives:

`0.96 * (37.5 - 22) / (37.5 - 20) = 0.850 W/kg`

This extrapolation is a modeling assumption, not a directly measured 22 C value.

For the project's `22.91 g` mouse:

`Q_E2,dry(22 C) = 0.850 W/kg * 0.02291 kg = 0.01948 W = 0.01675 kcal/h`

With the project's whole-body heat capacity `C_b = 0.01846324 kcal/C`, this is equivalent to:

`-(0.01675 / 0.01846324) = -0.907 C/h`

when the empirical E2 replacement state is fully on.

Thus the simplest empirical term would be:

`dT/dt ... - 0.907 * E2_state   [C/h]`

or, in the heat balance before division by heat capacity:

`Q_loss ... + 0.01675 * E2_state   [kcal/h]`

### 3. If the existing E2 Hill gate must be reused

The current project gate is:

`H(E2) = E2^0.684908 / (20.146774^0.684908 + E2^0.684908)`

At the source experiment's `49.4 pg/mL`, `H = 0.648921`. To reproduce the 22 C extrapolated flux at that concentration:

`E2_DRY_LOSS_GAIN_KCAL_H = 0.0167498 / 0.648921 = 0.0258118 kcal/h`

and:

`Q_E2,dry = 0.0258118 * H(E2)   [kcal/h]`

This is a model translation only. Neither the Hill midpoint nor the Hill exponent was measured for dry heat loss.

## Recommendation for the present model

Use `0` for the primary non-CVC E2 loss coefficient. If the model needs to test whether an unobserved E2-associated dry-loss arm could counter the existing thermogenic term, add the rat-derived coefficient only as a clearly labeled sensitivity case:

- Base model: `E2_NONCVC_LOSS_GAIN_KCAL_H = 0`
- Empirical rat-prior sensitivity: `E2_DRY_LOSS_GAIN_KCAL_H = 0.0258118` multiplying the existing Hill gate
- Direct binary interpretation at 22 C: `0.01675 kcal/h` for physiological E2 replacement versus OVX

Do not describe the sensitivity term as evaporative loss, fur loss, behavioral loss, or non-CVC loss. The evidence supports only **whole-body E2-associated dry heat loss of unresolved mechanism**.

## Primary sources

1. Laudenslager ML, Wilkinson CW, Carlisle HJ, Hammel HT. Energy balance in ovariectomized rats with and without estrogen replacement. *Am J Physiol*. 1980;238:R400-R405. [DOI](https://doi.org/10.1152/ajpregu.1980.238.5.R400)
2. Laudenslager ML, Carlisle HJ, Calvano SE. Increased heat loss in ovariectomized hypothyroid rats treated with estradiol. *Am J Physiol*. 1982;243:R70-R76. [DOI](https://doi.org/10.1152/ajpregu.1982.243.1.R70)
3. Krajewski-Hall SJ, Blackmore EM, McMinn JR, Rance NE. Estradiol alters body temperature regulation in the female mouse. *Temperature*. 2018;5:56-69. [DOI](https://doi.org/10.1080/23328940.2017.1384090)
4. Baker MA, Dawson DD, Peters CE, Walker AM. Effects of estrogen on thermoregulatory evaporation in rats exposed to heat. *Am J Physiol*. 1994;267:R673-R677. [DOI](https://doi.org/10.1152/ajpregu.1994.267.3.R673)
5. Hosono T, Chen XM, Miyatsuji A, et al. Effects of estrogen on thermoregulatory tail vasomotion and heat-escape behavior in freely moving female rats. *Am J Physiol Regul Integr Comp Physiol*. 2001;280:R1341-R1347. [DOI](https://doi.org/10.1152/ajpregu.2001.280.5.R1341)
6. Marui S, Masuda Y, Kato I, Nagashima K. Influence of exogenous and endogenous estrogen on thermoregulatory responses to mild heat and the interaction with light and dark phases. *J Physiol Sci*. 2020;70:56. [DOI](https://doi.org/10.1186/s12576-020-00782-x)
7. Kaikaew K, Steenbergen J, Themmen APN, Visser JA, Grefhorst A. Sex difference in thermal preference of adult mice does not depend on presence of the gonads. *Biol Sex Differ*. 2017;8:24. [DOI](https://doi.org/10.1186/s13293-017-0145-7)

## Scope limits

- The only exact near-room-temperature flux coefficient comes from rats, not mice.
- It reflects chronic replacement, not minute-to-minute E2 fluctuations.
- The source comparison is replacement versus OVX, not a dose-response curve.
- The 22 C conversion is a linear dry-gradient extrapolation from 20 C.
- Direct mouse evidence at exactly 22 C does not identify a heat-loss route.
- Pregnancy cannot validate this coefficient because body mass, geometry, progesterone, activity, metabolism, and insulation change simultaneously.
