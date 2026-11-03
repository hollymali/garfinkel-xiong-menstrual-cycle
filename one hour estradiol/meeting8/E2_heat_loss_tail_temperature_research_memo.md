# Estradiol, heat loss, and mouse-tail temperature

## Evidence review and modeling plan for Meeting 8

**Prepared:** September 2026  
**Question:** How strongly does estradiol (E2) alter heat loss, and how can paired mouse tail and core temperatures be used to infer vasomotor state in pregnancy, inducible *Esr1* knockout, and cycling nonpregnant conditions?

---

## Executive answer

**Takeaway for Ivy's mice at 22 °C:** strong constriction is the expected baseline, but the literature does **not** justify assuming that every freely moving female mouse is continuously or maximally constricted. Two non-E2 mouse studies support near-complete constriction at 22 °C: restrained male CD-1 mice had tail temperature near ambient below 24 °C,[27] and freely moving male C57BL/6J mice had a mean tail temperature of 23.5 ± 0.9 °C at 22 °C.[9] In contrast, the closest female E2 study found appreciable within-hour HLI range at 22 °C and explicitly did not see the low variability expected from constant maximal constriction.[1] Thus the state at 22 °C is protocol-, strain-, sex-, and context-dependent; Ivy's own time series must show whether there is enough dynamic range to analyze events.

Continuous low-dose E2 from a SILASTIC capsule did not change the **third-hour mean** tail temperature, mean HLI, or one-hour HLI range in mice at exactly 22 °C.[1] The capsule initially produced physiological levels, but exposure declined during the 9-day chamber series (serum E2 was 6.8 ± 0.46 pg/mL on day 17). The study also did not count natural dilation episodes, measure their dwell time, or quantify their area. Thus it shows no detectable E2 effect on those summary endpoints at 22 °C; it does not establish that E2 has no effect on event frequency or duration.

Use paired tail, core, and ambient temperatures to calculate HLI and analyze natural tail-warming events as a secondary, prespecified phenotype. HLI is a tail vasomotor proxy, not cutaneous vascular conductance itself: true vascular conductance requires blood-flow and pressure measurements. The literature supports asking whether pregnancy, cycling E2, or *Esr1* deletion changes event occupancy, frequency, duration, magnitude, or core-temperature coupling, but it does not supply a published mouse event-rate coefficient at 22 °C.

For Meeting 8:

- Do not add a fixed E2 multiplier to whole-body conductance.
- Add tail temperature/HLI as an observation of a latent vasomotor state.
- Add a bounded tail heat-loss term only if HLI improves prediction of later core-temperature change.

This fits the closest mouse evidence: continuous E2 produced no home-cage mean tail/HLI change at 23–25 °C and no mean or range effect at exactly 22 °C, but it affected tail/HLI at selected nearby chamber temperatures.[1] Tail removal at 22 °C also had little effect on whole-body thermal balance, and intact tails were only about 1.5 °C above ambient at baseline, showing that the mouse tail is usually a small, strongly constricted heat exchanger in the cold.[9] The null E2 findings followed sustained capsule exposure over days, not a single injection.[1]

---

## 1. What “E2 affects heat loss” can mean

Several distinct quantities are called heat loss in this literature. They should not be treated as interchangeable.

1. **Tail skin temperature, `Ttail`.** This is an accessible surface measurement. It rises with increased warm blood delivery, but it also changes passively with ambient temperature (`Ta`), core temperature (`Tcore`), sensor position, posture, contact with bedding, and local convection.

2. **Heat-loss index, HLI.** The standard normalized index is

   `HLI = (Ttail - Ta) / (Tcore - Ta)`.

   It is approximately 0 when the tail approaches ambient temperature during strong constriction and approaches 1 when tail temperature approaches core temperature during strong dilation.[10] HLI removes much of the passive effect of ambient and core temperature, although it is still not blood flow or an absolute heat flux.

3. **Tail heat flux.** A first physical approximation is

   `Qtail ≈ Atail [hc (Ttail - Ta) + εσ((Ttail + 273.15)^4 - (Ta + 273.15)^4)]`,

   with additional conductive loss when the tail contacts bedding or a surface. This requires tail area, emissivity, convection, posture/contact, and environmental calibration.

4. **Whole-body dry thermal conductance.** A common indirect-calorimetry estimate is metabolic heat production divided by the core-to-ambient gradient. It includes the entire animal, not just the tail.

5. **A thermoregulatory threshold.** E2 may change the core or ambient temperature at which vasodilation begins, without changing maximum dilation or maximum constriction. This is the pattern most clearly supported by the controlled rat study.[2]

6. **Evaporative loss.** At high ambient temperature, salivation and grooming can dominate. Tail temperature alone will not capture this route. Dacks and Rance explicitly raised increased evaporative cooling as a possible explanation for lower E2-treated core temperature at very high ambient temperature.[2]

Consequently, “no effect on whole-body conductance” does not imply “no effect on tail vasomotion,” and “lower tail temperature” does not automatically imply a meaningful change in whole-body heat loss.

---

## 2. Primary rodent evidence

### 2.1 Evidence matrix

| Study | Model and E2 exposure | Endpoint | Quantitative result | Interpretation |
|---|---|---|---|---|
| Krajewski-Hall et al., 2018[1] | OVX CD-1 mice; subcutaneous SILASTIC E2 capsule, initially 180 μg/mL; higher 360 and 720 μg/mL capsule conditions studied sequentially | Home-cage `Ttail`, `Tcore`, HLI; 20–36 °C chamber series | At exactly 22 °C, E2 did not change the third-hour mean `Ttail`, mean HLI, or one-hour HLI range. The HLI range remained appreciable and was not suppressed at low ambient temperatures. HLI was lower with E2 at 23, 25, and 26 °C. No intact-cycle tail/HLI effect in 6 mice across 16 cycles. | The mouse tail was still dynamically variable at 22 °C, but the published summaries do not show an E2 effect there. No natural events were counted. |
| Ding et al., 2013[19] | OVX 6-month-old 3×TgAD mice; E2 50 μg/kg subcutaneously every day for 5 weeks, started at OVX | Tail temperature by infrared thermometer twice weekly during the sleep phase | At week 5, OVX mice were 23.24 ± 0.08 °C vs 22.28 ± 0.07 °C in sham mice, a +0.96 °C difference; E2 prevented the OVX rise. | Supports a delayed/chronic mouse tail effect, but the disease-model background, sparse sampling, and lack of paired core/ambient measurements limit transferability. |
| Shuto et al., 2005[20] | OVX ICR mice; intramuscular estradiol valerate 1 mg/kg once weekly for 3 weeks, starting 7 days after OVX | Tail-temperature response after 10-min treadmill exercise at 25 °C | Estrogen reduced the integrated 1-6-min postexercise tail-temperature rise by 81.3 ± 11.9%. | E2 strongly suppresses an evoked dilation-like tail transient even when baseline mean effects can be difficult to detect. |
| Fleischer et al., 2021[22] | OVX C57BL/6 mice; E2 0.2 mg/kg by daily oral gavage for 63 days | Baseline thermal-imaging `Ttail`; senktide-evoked response | E2 lowered baseline tail temperature after 63 days (`p<0.05`) and reduced senktide-evoked tail warming at multiple postinjection time points. | Long exposure can reveal a mouse baseline effect and suppress an evoked vasodilatory response; absolute effect size was not tabulated in the text. |
| Dacks & Rance, 2010[2] | OVX Sprague-Dawley rats; two subcutaneous SILASTIC capsules containing E2, recorded 9-25 days after implantation | HLI over 13-34 °C; thermoneutral-zone boundaries | E2 lowered HLI overall; pointwise differences at 24.4 and 29 °C, but not 22 °C. Lower vasodilation boundary shifted from 21.9 to 26.1 °C (+4.2 °C). Basal core temperature below 32.5 °C was essentially unchanged. | Strong evidence for a threshold shift toward sustained constriction, not a constant conductance multiplier. |
| Williams et al., 2010[3] | OVX rats; continuous subcutaneous E2 capsule, monitored for 21 days at 21.1-22.5 °C | Tail and core temperature every 5 min | E2 reduced dark-phase tail temperature from days 2-21, while core temperature was not different. Adding progesterone did not reverse this effect. | Prolonged recording detects an active-phase vasomotor effect even near 22 °C. |
| Marui et al., 2020[4] | OVX Wistar rats; subcutaneous silicone tube with 50-60 mg E2 powder, constant exposure for at least 14 days; separate intact diestrus/proestrus comparison | Home-cage tail temperature; HLI and whole-body conductance during 25-34 °C heat challenges | At 25 °C, dark-phase tail temperature was 29.7 ± 0.8 °C with E2 vs 33.5 ± 0.7 °C in OVX controls (−3.8 °C). Proestrus was 1.3 °C lower than diestrus in the dark phase. During 28, 31, and 34 °C challenges, HLI and whole-body conductance did not differ between OVX and E2 groups. | E2 strongly altered baseline dark-phase tail temperature but not the response endpoint during warmer challenges. The sign and detectability depend on phase and assay. |
| Okada et al., 1997[5] | OVX Wistar rats; subcutaneous osmotic minipump delivering 0.3 or 1.0 μg E2 per animal per day for 14 days | Tail skin temperature | OVX elevated tail temperature by day 6; both continuous E2 doses significantly reduced it after treatment. | Continuous E2 restored a lower-tail-temperature phenotype consistent with greater vasoconstrictor tone. |
| Bowe et al., 2006[21] | OVX Wistar rats; E2 4 μg/kg subcutaneously every day; separate oral E2 75 μg/day condition | Dark-phase tail temperature | Subcutaneous E2 significantly lowered tail temperature by treatment day 2; it returned to baseline about 4 days after stopping. Oral E2 was significant by day 2 and returned to baseline about 2 days after stopping. | Demonstrates an onset/offset time course and receptor dependence, not an instantaneous one-dose effect. |
| Hosono et al., 2001[6] | OVX Wistar rats; subcutaneous Silastic tube containing 2 mg estradiol benzoate; serum E2 95.6 ± 5.6 pg/mL | Freely moving tail/core temperature every minute for 3 h at 13–31 °C | E2 significantly lowered mean `Ttail-Ta` at exactly 22 °C. “Large” dilation was classified as a >5 °C session range, but only the number of rats with at least one such excursion was reported. More than half of E2-treated rats still met this criterion at 22 °C. | E2 shifts the ambient-temperature window for on/off tail switching. It does not abolish natural dilation at 22 °C, and this paper gives incidence over 3 h—not events/hour. |
| Opas et al., 2006[7] | OVX WT, ERαKO, and ERβKO mice; subcutaneous estradiol cypionate, stated as an approximately 0.1 mg/kg/day dose | Basal tail skin temperature | Ovariectomy elevated tail temperature; estradiol cypionate suppressed it dose-dependently. Suppression remained in either single-receptor knockout. | A tail response to E2 is not a clean readout of ERα alone; receptor redundancy or compensation matters. |
| Zhang et al., 2021[8] | Female WT and hypothalamic *Esr1* conditional-KO mice; tamoxifen 0.1 mg/kg subcutaneously every day for 4 weeks vs corn oil | Tail/core telemetry-derived HLI | WT treatment effect: `F(1,14)=15.52`, `p=0.0015`; 24-h `p=0.005`. Conditional KO: `F(1,13)=3.109`, `p=0.1013`. | Mechanistic evidence that hypothalamic ERα is required for the chronic tamoxifen HLI effect. It is not an E2 dose-response study. |

### 2.2 The cleanest quantitative E2 result: a vasodilation-threshold shift

Dacks and Rance provide the most directly useful mechanistic number for model design.[2] They ovariectomized female rats and supplied E2 continuously through two subcutaneous SILASTIC implants. The mean serum E2 concentration was 24.8 ± 0.6 pg/mL. Tail and core temperatures were assessed over 23 ambient temperatures from 13 to 34 °C during the light phase, 9-25 days after implantation.

Using overlapping confidence intervals to define the thermoneutral region, the lower boundary for tail vasodilation was 21.9 °C in OVX controls and 26.1 °C with E2. Thus, E2 shifted the onset of tail vasodilation upward by about 4.2 °C. Importantly:

- HLI was significantly lower at 24.4 and 29 °C.
- The groups were not significantly different at 22 °C.
- Basal core temperature was nearly identical below 32.5 °C: 36.87 ± 0.03 °C in OVX controls and 36.91 ± 0.03 °C with E2.
- At the cold and hot extremes, maximal constriction/dilation appeared similar; the major effect lay in the transition region.

For a model, this favors a representation such as

`Pr(dilated) = logistic[(signal - threshold(E2, phase, genotype)) / slope]`

over

`whole-body conductance = baseline conductance × (1 - constant × E2)`.

The 4.2 °C value should **not** be inserted directly into a mouse model. It is a rat, light-phase, chronic-replacement result and is expressed on the ambient-temperature axis. It is best used as a prior expectation that E2 may shift a switching threshold by several degrees in a permissive thermal range.

### 2.3 Why rat studies near 22-25 °C can appear inconsistent

The apparent contradictions are informative:

- Dacks and Rance found no pointwise difference at 22 °C in a light-phase ambient-temperature protocol.[2]
- Williams et al. observed lower tail temperature throughout the **dark phase** with continuous E2 at 21.1-22.5 °C using five-minute sampling across 21 days.[3]
- Marui et al. found a large dark-phase difference at 25 °C but little light-phase difference, and no group difference in HLI or whole-body conductance during separate warm challenges.[4]

This suggests that phase, arousal/activity, acclimation, measurement duration, and whether the endpoint is a mean, a transition, or a response slope can determine whether an E2 effect is detected. A single daytime average is particularly likely to miss an active-phase effect.

### 2.4 Mouse-specific evidence

Krajewski-Hall et al. is the most relevant direct mouse E2 study.[1] Adult OVX CD-1 females received a subcutaneous continuous-release E2 capsule. Serum E2 was 25.3 ± 4.2 pg/mL on day 7, 17.8 ± 1.6 pg/mL on day 14, and 6.8 ± 0.46 pg/mL on day 17, showing that the implant exposure declined over time.

During home-cage monitoring, E2 lowered light-phase core temperature but did not significantly alter tail temperature or HLI. Reported means for the sequential OVX, 360 μg/mL, and 720 μg/mL capsule conditions were:

| Endpoint | OVX | E2 360 μg/mL | E2 720 μg/mL |
|---|---:|---:|---:|
| Light `Tcore` (°C) | 36.6 ± 0.04 | 36.3 ± 0.07 | 36.3 ± 0.05 |
| Light `Ttail` (°C) | 31.4 ± 0.32 | 30.5 ± 0.22 | 31.5 ± 0.38 |
| Light HLI | 0.63 ± 0.02 | 0.58 ± 0.01 | 0.65 ± 0.03 |
| Dark `Ttail` (°C) | 29.8 ± 0.46 | 29.2 ± 0.15 | 30.1 ± 0.48 |
| Dark HLI | 0.49 ± 0.03 | 0.47 ± 0.01 | 0.50 ± 0.03 |

Those home-cage comparisons were null for the tail endpoints. Yet, during an ambient-temperature series, E2-treated mice had lower tail temperature at 21, 23, 25, and 26 °C and lower HLI at 23, 25, and 26 °C. E2 also abolished the acute tail-temperature rise produced by the neurokinin-B receptor agonist senktide. In OVX mice, senktide increased tail temperature by 2.64 ± 0.96 °C at approximately 13 minutes and decreased core temperature by 1.74 ± 0.25 °C at approximately 21 minutes. This is useful validation that their tail sensor captured physiologically meaningful dilation even though the chronic home-cage mean was null.

The Krajewski-Hall experiments answer several different questions and should not be collapsed into one result:

- **Home cage, 23–25 °C:** every 10 minutes for five days; published tests used six-hour light/dark blocks. E2 altered light-phase core temperature but not mean tail temperature or HLI. Two- and fourfold higher capsule concentrations also failed to produce a monotonic tail/HLI effect.
- **Ambient series, exact 22 °C included:** every five minutes for three hours; only the third hour was analyzed. At 22 °C, E2 changed core temperature but not mean tail temperature, mean HLI, or HLI range. The nonzero HLI range and the authors' statement that low ambient temperatures did not suppress variability argue against constant maximal constriction. A caveat is that the 180 µg/mL capsule's output declined during testing: serum E2 was 17.8 ± 1.6 pg/mL on day 14 but only 6.8 ± 0.46 pg/mL on day 17, and the exact day on which each randomized ambient temperature was tested was not reported.
- **Intact cycles, home cage 23–25 °C:** every five minutes for 14 days; six mice with regular cycles contributed 16 cycles. No cycle-stage difference in six-hour mean tail temperature or HLI was detected. This is not a direct serum-E2 analysis, and it does not test event counts.
- **Senktide, 23–24 °C:** an injected NK3R agonist evoked an acute tail-warming response sampled every minute, and E2 strongly attenuated it. The 2.64 °C value is therefore an evoked response amplitude, not the amplitude or frequency of natural dilation events.

In practical terms, a mouse tail dataset may contain the signal in **transitions, phase-specific effects, or conditional responses** even when its 24-hour average does not.

Two additional mouse studies help reconcile the mixed baseline results. Ding et al. gave 50 μg/kg E2 subcutaneously every day for five weeks starting at OVX in 3×TgAD mice.[19] OVX raised tail temperature by about 0.96 °C relative to sham at week 5, and E2 prevented the rise. This was intermittent infrared sampling without paired core or ambient temperature, so it establishes a chronic tail phenotype but not heat flux. Fleischer et al. found lower baseline tail temperature after 63 days of daily oral E2 at 0.2 mg/kg and a smaller senktide-evoked tail response.[22] Together with Krajewski-Hall, the mouse evidence suggests that exposure duration and whether the tail is challenged may matter as much as nominal dose.

Shuto et al. provide the clearest evoked mouse effect size.[20] Estradiol valerate, 1 mg/kg intramuscularly once a week for three weeks, reduced the integrated tail-temperature rise over the first six minutes after a 10-minute treadmill bout by 81.3 ± 11.9%. This does not mean E2 reduced total heat loss by 81%; it means it markedly suppressed this particular transient tail-temperature endpoint.

### 2.5 What is actually known about natural dilation events at approximately 22 °C

The literature supports the existence and physiological relevance of spontaneous tail-temperature transitions, but it does **not** provide a clean E2-dependent event rate for mice at 22 °C.

The broader, non-E2 literature changes the conclusion from “not maximally constricted” to “usually strongly constricted, with experiment-dependent exceptions”:

| Study | Species and conditions | Temporal resolution | What it says about spontaneous events | Applicability to Ivy's mice |
|---|---|---:|---|---|
| Pincedé et al., 2012[27] | Male Swiss CD-1 mice; 20–32 °C; two hours in a Plexiglas restrainer | 1 second | Below 24 °C, tail temperature was close to ambient in all 15 mice, interpreted as complete vasoconstriction. Above 25 °C, spontaneous fluctuations occurred at about 2–3 cycles/h and reached 8 °C. | Direct mouse evidence that events may be absent at 22 °C, but restraint, male sex, and a changing ambient-temperature protocol limit transfer to freely moving females. |
| Škop et al., 2020[9] | Male C57BL/6J mice; open cage with bedding at 22 °C after one-hour acclimation | Approximately 30-second IR videos collected intermittently | Baseline tail temperature was 23.5 ± 0.9 °C at 22 °C; α1 blockade with prazosin raised it by up to 4.6 ± 0.8 °C. This confirms strong baseline vasoconstrictor tone and preserved dilation capacity. | Strong direct evidence for a low baseline at 22 °C, but the intermittent clips cannot determine event frequency or prove continuous maximal constriction. |
| Krajewski-Hall et al., 2018[1] | Freely moving female CD-1 mice; exact 22 °C; OVX ± continuous E2 | 5 minutes for 3 hours | One-hour HLI range remained appreciable; mean HLI, mean tail temperature, and range did not differ with E2. Events were not counted. | Closest endocrine/sex match. It indicates dynamic variability in this setup but does not establish discrete event frequency. |
| Williams et al., 2010[3] | Freely moving female rats; 21.1–22.5 °C; OVX ± continuous E2 | 5 minutes for 21 days | Spontaneous large fluctuations were visible and E2 lowered dark-phase mean tail temperature, but six-hour means/SD were analyzed and Fourier analysis found no stable ultradian rhythm. | Direct room-temperature precedent that irregular events can exist and carry an E2 effect, but it is rat rather than mouse. |
| Hosono et al., 2001[6] | Freely moving female rats; exact 22 °C; OVX ± continuous estradiol benzoate | 1 minute for 3 hours | More than half of E2-treated rats had at least one >5 °C session excursion. The paper reports animal-level incidence, not episodes per hour. | High-resolution exact-22 °C evidence for events in rats, but not a transferable mouse rate. |

Published numerical event rates come from warmer, thermoneutral conditions—not from freely moving mice at 22 °C. Pincedé et al. observed about 2–3 cycles/h in mice only above 25 °C.[27] Young and Dawson observed an approximately 20-minute on/off cycle (about 3/h) in rats at 29–33 °C.[28] El Bitar et al. measured 2–7 cycles/h, with a mean near 3–4/h, in thermoneutral rats.[11] These studies validate the existence of a switch-like event process and provide useful expectations for event duration, but their rates should **not** be imposed as a prior rate for Ivy's cold-housed mice.

**Krajewski-Hall et al., mouse, exact 22 °C.[1]** The chamber experiment is the closest match to Ivy's animals. OVX and OVX+E2 mice were exposed to each ambient temperature for three hours, with tail, core, ambient temperature, and activity sampled every five minutes. The published analysis discarded the first two hours for acclimation and summarized the third hour as mean `Ttail`, mean HLI, and HLI range (`maximum HLI - minimum HLI`). At 22 °C, none of those endpoints differed with E2. Crucially, the HLI range was not near zero: across the low-temperature conditions, the authors did not observe the reduced range expected from constant maximal vasoconstriction. Therefore this paper argues against the claim that mouse tail vasomotion is already fully saturated at 22 °C. It cannot answer whether E2 changes the number, duration, or shape of events because those quantities were not computed.

**Hosono et al., rat, exact 22 °C.[6]** Freely moving OVX rats with continuous estradiol-benzoate or vehicle capsules were sampled every minute for three hours after a two-hour acclimation. Tail vasomotion occurred in an on/off fashion. E2 lowered the three-hour mean `Ttail-Ta` at 19, 22, and 25 °C. For an event-like endpoint, the investigators classified a rat as having a “large vasodilation” when its session maximum-minus-minimum tail temperature exceeded 5 °C. At 19 °C, 11/16 low-E2 rats versus 5/18 E2-treated rats met that criterion. At 22 °C, the text reports that more than half of the E2-treated rats showed vasodilation. The figure suggests that both groups commonly had a large excursion at 22 °C, but the paper did not report a number of episodes per rat. This pattern is best understood as an E2-dependent rightward shift of the switching window: E2 can reduce mean time spent dilated while still allowing, or even positioning an animal to show, transitions at 22 °C.

**Williams et al., rat, 21.1–22.5 °C.[3]** Tail temperature was sampled every five minutes for 21 days. Natural large-amplitude fluctuations were visible, and continuous E2 lowered dark-phase tail temperature beginning about 35–48 hours after implantation. However, the inferential analyses used six-hour means and six-hour standard deviations. The investigators also performed Fourier analysis over days 3–7 and found no significant ultradian rhythm across treatment groups. This is evidence that the natural fluctuations are irregular rather than a stable periodic oscillation with a literature-defined frequency. It does not establish zero events or equal event counts.

**Dacks et al., rat, 20 °C.[25]** Tail temperature was sampled every five seconds in freely moving OVX rats during one- to three-hour sessions. Continuous physiological E2 lowered the session-average HLI at 20 °C (approximately 0.17 with E2 versus 0.24 without E2 from the published figure). This is strong time-resolved acquisition below 22 °C, but the raw resolution was collapsed to a session mean; event frequency and dwell time were not reported.

**Closest explicit event-frequency precedent, but not at 22 °C or under controlled E2 replacement.[26]** Simpkins monitored spontaneous tail-temperature surges every five minutes for two hours at 24–25 or 29–30 °C in young cycling and aged rats of different reproductive states. Frequency and incidence differed by reproductive status, but age, estrogen, progesterone, restraint, and reproductive state were not experimentally separated. It shows that counting spontaneous tail events is an established phenotype, not that a specific E2 event-rate coefficient can be transferred to Ivy's mice.

The resulting evidence boundary is precise:

- **Known:** mouse tails are generally strongly constricted at 22 °C, yet at least one freely moving female-mouse protocol retained appreciable HLI variability; continuous E2 can lower average tail heat-loss proxies at 20–22 °C in rats; E2 does not alter the published one-hour mean or range endpoints in the closest mouse experiment.
- **Not known:** natural mouse dilation episodes per hour under E2 at exactly 22 °C, or whether E2 acts mainly through event frequency, duration, magnitude, or baseline between events.
- **Implication for Ivy's data:** first test whether the traces contain reproducible excursions above a strongly constricted baseline. If they do, event analysis is a justified secondary analysis with direct physiological precedent, but its numerical E2 coefficient must be estimated from Ivy's mice rather than borrowed from the literature. If the traces are flat and close to ambient, the event rate is effectively near zero under Ivy's conditions and group comparisons should focus on HLI level/occupancy rather than event counts.

### 2.6 Whole-body dry and evaporative heat loss

An older direct energy-balance study complicates any simple statement that E2 always decreases heat loss. In dark-phase OVX rats, chronic 17β-estradiol or estradiol-benzoate replacement through subcutaneous Silastic capsules increased both metabolic heat production and measured dry heat loss, while rectal temperature was unchanged.[23] Evaporative loss was unchanged. When replacement was delayed until 75 days after OVX, the increases in production and dry loss were absent. The accessible report does not provide a reliable, extractable numerical treatment difference, but its direction is important: whole-body dry loss can increase secondarily to a larger internally generated heat load even while E2 suppresses tail vasodilation at a matched thermal state.

At a much warmer 38 °C, E2 has a different action. In a within-rat, opposite-treatment-order experiment conducted over nine weeks, an E2 pellet increased evaporative water loss at all measured core temperatures, lowered the core-temperature threshold for saliva spreading, and allowed a lower core temperature during heat exposure.[24] The paper does not report a simple transferable percentage in its abstract. This route of heat loss will be invisible to a tail-only observation model.

The integrated interpretation is therefore:

- at ordinary subthermoneutral-to-midrange conditions, E2 often suppresses tail vasodilation or shifts its onset to warmer conditions;
- whole-body dry loss can nevertheless rise when E2 also raises heat production;
- at severe heat, E2 can facilitate evaporative defense;
- the net effect on core temperature is the balance of all three processes.

### 2.7 Numerical anchors that are useful - and those that are not

The following are the most useful numerical anchors for planning and sensitivity analysis:

| Quantity | Best available estimate | Appropriate use |
|---|---:|---|
| Rat lower threshold for tail vasodilation | +4.2 °C with chronic physiological E2 (21.9→26.1 °C)[2] | Supports a threshold-shift model; use only as a broad prior on mechanism, not a mouse parameter. |
| Rat dark-phase tail temperature at 25 °C | −3.8 °C with chronic high E2 (33.5→29.7 °C)[4] | Shows the possible size of a phase-specific surface-temperature effect. The serum E2 was high and the species is rat. |
| Intact rat proestrus vs diestrus dark tail temperature | −1.3 °C (32.7→34.0 °C comparison)[4] | A physiological-cycle benchmark in rat, confounded by other cycle hormones. |
| Mouse OVX-associated chronic tail rise | +0.96 °C relative to sham at week 5; daily E2 prevented it[19] | Rough scale for a slowly emerging mouse phenotype; not a paired core/ambient estimate. |
| Mouse exercise-evoked tail transient | 81.3 ± 11.9% suppression by weekly estradiol valerate[20] | Establishes that an event-level effect can be large. It is not percent whole-body heat loss. |
| WT mouse tamoxifen effect on 24-h HLI | +0.094 in the source-data reanalysis[8] | A useful HLI analysis benchmark and pipeline check, but not an E2 coefficient. |
| Maximum mouse tail contribution to whole-body loss | approximately 5-8% under severe heat stress[9] | Upper-bound sensitivity constraint when coupling tail state to core heat balance. |

A pooled “E2 effect size” across these studies would not be scientifically meaningful. The quantities differ in species, ambient temperature, endocrine state, route, duration, circadian phase, and whether the endpoint is a baseline, threshold, transient, or whole-body flux. The literature should constrain the **model form and plausible scale**, while the PI's paired mouse data estimate the cohort-specific parameters.

### 2.8 Human evidence is mechanistically supportive, not a mouse calibration

Human studies confirm that estrogen can alter thresholds and cutaneous vascular responses, but the direction depends on the population and endpoint. In a randomized, double-blind trial of 24 symptomatic postmenopausal women, oral 17β-E2 at 1 mg/day for 90 days raised the core-temperature threshold for sweating from 37.98 ± 0.09 to 38.14 ± 0.09 °C (+0.16 °C) and reduced hot-flash frequency.[15] This was chronic treatment, not a one-time dose.

After 10-12 days of ovarian-hormone suppression in 25 healthy young women, transdermal E2 at 0.1-0.2 mg/day for 3-4 days increased the initial local-heating cutaneous vascular conductance peak from 71 ± 2% to 79 ± 2% of maximum, while the later plateau, baseline, maximum flux, and reactive hyperemia were unchanged.[16] Thus, even within one study, the answer was positive for one phase of the vascular response and null for others.

In a chronic observational comparison, women receiving estrogen-only replacement for at least two years initiated cutaneous vasodilation at a body temperature about 0.4 °C below untreated postmenopausal women (36.5 ± 0.1 vs 36.9 ± 0.1 °C); combined estrogen-progestin therapy did not show the same pattern.[17] Differences among oral systemic therapy, local skin heating, exercise, symptomatic status, progesterone exposure, and the chosen threshold explain why these findings should not be compressed into a single sign. Their value for this project is conceptual: estrogen frequently changes the **activation threshold or early response**, rather than maximum vascular capacity.

---

## 3. Explicit audit of studies reporting no E2 effect

The null findings are not acute, one-time E2 tests. Most used continuous or repeated exposure.

| Null or partly null result | How E2 was administered | Exposure before/during endpoint | What was actually null | What was not null |
|---|---|---|---|---|
| Krajewski-Hall et al., mouse home cage[1] | Subcutaneous SILASTIC continuous-release capsule; initial 180 μg/mL E2 in capsule. The same animals later received 360 and 720 μg/mL capsule conditions sequentially. | Initial circadian endpoints on days 3-5; higher-dose conditions lasted 7 days each, separated by washout. Chamber tests occurred over days 13-21 after the initial treatment series. | Home-cage `Ttail`, HLI, and activity at the initial dose; tail/HLI also lacked a monotonic higher-dose response. Intact estrous-cycle phase did not alter tail/HLI in 6 mice across 16 cycles. | E2 lowered light-phase core temperature; selected chamber temperatures showed lower tail/HLI; E2 blocked the senktide-induced tail rise. |
| Dacks & Rance at 22 °C[2] | Two subcutaneous continuous-release SILASTIC capsules containing E2 in sesame oil. | Measurements 9-25 days after implantation. | Pointwise HLI difference at 22 °C. | Significant overall HLI effect, significant differences at 24.4 and 29 °C, and a +4.2 °C shift in the vasodilation boundary. |
| Marui et al. warm-challenge response[4] | Subcutaneous silicone implant containing 50-60 mg crystalline E2; constant exposure. | At least 14 days; plasma E2 approximately 142 ± 25 pg/mL in light and 159 ± 31 pg/mL in dark. | HLI and calculated whole-body conductance during 28, 31, and 34 °C challenges; only 5 animals per group/temperature condition. | At the 25 °C home-cage baseline, E2 reduced dark-phase tail temperature by 3.8 °C and increased the day-night amplitude. |
| Laudenslager et al., delayed-replacement condition[23] | 17β-estradiol or estradiol benzoate through a subcutaneous Silastic capsule; chronic replacement, not a single injection. | Replacement was delayed until 75 days after OVX; the accessible abstract does not state the subsequent implant duration precisely. | The increased heat production and dry heat loss seen with earlier replacement were absent after delayed initiation. | Earlier chronic replacement increased both production and dry loss while core temperature remained unchanged. |
| Brunt et al., human local-heating plateau[16] | Transdermal 17β-E2, 0.1-0.2 mg/day, after pharmacologic suppression of ovarian hormones. | 3-4 days of patch treatment following 10-12 days of GnRH-antagonist suppression. | Baseline skin blood flow, local-heating plateau, maximum flux, and reactive hyperemia. | The initial local-heating peak rose from 71 ± 2% to 79 ± 2% of maximum cutaneous vascular conductance; the nadir rose from 57 ± 2% to 65 ± 2%. |
| Late-pregnancy rat whole-body conductance[13] | No administered E2; endogenous pregnancy physiology. | Rats tested at gestational days 10/11, 15/16, and 20/21. | No overall main effect of pregnancy on calculated thermal conductance across the full test. | Late pregnancy shifted the threshold for increasing conductance from about 30 to about 28 °C and lowered core temperature by about 1 °C. |

The key point is that a reported null often means “no difference in this endpoint, in this phase, at this ambient temperature, after days or weeks of E2.” It does not mean that E2 was given once and had no thermoregulatory action.

---

## 4. What the estrogen-receptor studies imply for the inducible *Esr1* KO

The exact definition of the PI's “Esr1i KO” is essential. At least three biologically different comparisons could be meant:

- a global inducible deletion of *Esr1*;
- an adult brain- or hypothalamus-specific inducible deletion;
- a cell-type-specific inducible deletion in a thermoregulatory neuronal population.

They do not have the same prediction.

Opas et al. reported that estradiol cypionate lowered basal tail skin temperature in both ERαKO and ERβKO mice.[7] The authors interpreted this to mean that either receptor subtype can support the estrogenic suppression of tail temperature. Because these were constitutive single-receptor knockouts, developmental compensation is also possible. Therefore, “E2 still affects the tail in an ERα knockout” would not, by itself, invalidate an estrogen mechanism.

Zhang et al. used *Esr1* floxed mice crossed to Nkx2-1-Cre, which removes ERα from much of the hypothalamus/preoptic region.[8] Chronic low-dose tamoxifen increased HLI in wild-type controls, but the treatment effect was not significant in conditional knockouts. This supports a central ERα pathway for the tamoxifen response, but it does not establish that physiological E2 acts identically.

### Reanalysis of Zhang et al. source data

The published source-data workbooks were reanalyzed descriptively. These are group differences in HLI after four weeks of tamoxifen versus oil, not E2 effects.

| Genotype and phase | Oil mean ± SD | Tamoxifen mean ± SD | Difference | Descriptive Hedges' g |
|---|---:|---:|---:|---:|
| WT, light | 0.500 ± 0.025 | 0.597 ± 0.063 | +0.097 | +1.90 |
| WT, dark | 0.346 ± 0.082 | 0.439 ± 0.056 | +0.093 | +1.25 |
| WT, 24 h | 0.422 ± 0.044 | 0.516 ± 0.052 | +0.094 | +1.85 |
| cKO, light | 0.657 ± 0.056 | 0.611 ± 0.058 | −0.047 | −0.76 |
| cKO, dark | 0.473 ± 0.129 | 0.377 ± 0.099 | −0.096 | −0.79 |
| cKO, 24 h | 0.564 ± 0.089 | 0.492 ± 0.068 | −0.072 | −0.86 |

The WT treatment effect was statistically supported in the original repeated-measures analysis; the cKO treatment effect was not. The opposite-signed cKO point estimates should not be overinterpreted because the sample is small (7-8 mice per group), the original omnibus test was nonsignificant, and the experiment tests tamoxifen rather than E2.

For the PI's dataset, the most defensible test is not simply “does KO lower the mean tail temperature?” It is whether the relationship among tail state, core temperature, ambient temperature, phase, and reproductive/E2 state is altered by the deletion. A central *Esr1* deletion may alter the switching threshold or state occupancy without changing the physical maximum or minimum tail temperature.

---

## 5. Pregnancy is not an E2-only intervention

No primary study was found that causally isolates E2 and simultaneously measures tail and core temperature across mouse pregnancy in a way that could directly calibrate this model. That is an important evidence gap.

In pregnant rats, Eliason and Fewell found that core temperature fell by roughly 1 °C during the second half of gestation without a corresponding overall change in selected ambient temperature, oxygen consumption, or calculated thermal conductance in the thermocline experiment.[13] In a metabolic-chamber protocol, late pregnancy shifted the approximate ambient-temperature threshold for increasing conductance from 30 to 28 °C, even though the overall conductance effect was not significant. This again supports analyzing **thresholds and control states**, rather than only daily mean conductance.

Pregnancy cannot be used as a clean surrogate for high E2. It simultaneously changes progesterone, prolactin, glucocorticoids, body mass and geometry, fetal/placental metabolism, insulation, activity, food intake, cardiovascular output, and often the circadian pattern. A pregnancy versus cycling comparison estimates the combined reproductive state. It identifies an E2 mechanism only if additional design features separate E2 from those covariates.

Recent mouse work links reproductive experience to altered thermal preference and changes in preoptic ERα neurons, but does not provide a tail heat-loss calibration during pregnancy.[14] It strengthens the biological rationale for an ERα/reproductive-state interaction while leaving the quantitative heat-loss question open.

Recommended language for Meeting 8 is therefore:

> Pregnancy-associated changes in HLI will be treated as reproductive-state effects. Their concordance with cycling E2 and their dependence on inducible *Esr1* loss can test an E2/ERα-consistent mechanism, but pregnancy alone cannot quantify an E2 coefficient.

---

## 6. Can tail temperature quantify vasoconstriction?

### Short answer

**Yes, as a relative index or latent state; no, as absolute percent vasoconstriction without calibration.**

Thermographic/surface temperature is physiologically linked to skin blood flow. In simultaneous laser-Doppler and thermal measurements in rat skin, the relationship was strong but temperature was smoothed and delayed relative to flow, with a mean lag of about 35 seconds and a mean correlation of about 0.73.[11] That supports tail temperature as a vasomotor observation while warning against interpreting every rapid temperature fluctuation as an instantaneous flow change.

### Recommended derived variable

For every synchronized observation, calculate

`HLI(t) = [Ttail(t) - Ta(t)] / [Tcore(t) - Ta(t)]`.

Then optionally define a communication-friendly **vasoconstriction index**

`VCI(t) = 1 - HLI(t)`.

VCI near 1 means relatively constricted and VCI near 0 means relatively dilated. It must be labeled a **dimensionless index**, not “percent vasoconstriction.”

Do not clip HLI to 0-1 before quality control. Values outside that interval can reveal sensor detachment, tail contact with a warm/cold surface, measurement mismatch, local heating, or an incorrectly assumed ambient temperature. A robust downstream model may winsorize or use an outlier component only after those episodes are identified.

HLI also becomes noisy when `Tcore - Ta` is small. Predefine a minimum core-to-ambient gradient, inspect sensitivity to that cutoff, and report warm-condition observations separately rather than allowing a near-zero denominator to create extreme values. This matters most during 34-36 °C challenges, not ordinary 22-25 °C housing.

### A useful quasi-steady relative-conductance proxy

A two-node steady-state tail balance gives

`Gblood (Tcore - Ttail) = Htail (Ttail - Ta)`.

Therefore,

`Gblood / Htail ≈ (Ttail - Ta) / (Tcore - Ttail) = HLI / (1 - HLI)`.

This ratio is potentially more mechanistic than HLI: it reflects internal blood-to-tail conductance relative to tail-to-environment conductance. But it is only valid approximately when:

- the tail is near thermal steady state;
- local tail metabolism is negligible;
- sensor position and tail geometry are fixed;
- convection, radiation, posture, and bedding contact are stable;
- the tail can be treated as a single well-mixed compartment.

It becomes numerically unstable as HLI approaches 1. It should be used for within-setup relative comparisons and sensitivity analysis, not reported as absolute perfusion.

### Why HLI is not whole-body heat loss

Removing the mouse tail did not materially change core temperature, energy expenditure, or whole-body conductance at 22 °C.[9] Prazosin raised tail temperature by as much as 4.8 ± 0.8 °C, yet the estimated tail contribution to heat dissipation was only 3.4% and was not statistically significant in that condition. Under severe heat stress at 30 °C with β3-adrenergic stimulation, tailless mice had core temperature 0.72 ± 0.14 °C higher and whole-body conductance 0.0058 ± 0.0021 kcal/h/°C lower. The authors estimated that the tail accounted for about 5-8% of maximal whole-body heat dissipation.[9]

Thus, the tail can be a sensitive **control signal** while being a modest **heat exchanger**. A large change in HLI need not produce a proportionally large change in core temperature.

Absolute rat-tail heat-flow models likewise require blood flow and experimental heat-loss information, and are sensitive to venous anatomy and countercurrent exchange.[12] A temperature-only dataset does not close that inverse problem.

---

## 7. Proposed model architecture

### 7.1 Keep production and peripheral control separable

The Meeting 7 model treated E2 primarily as a thermogenic input and removed a small E2-dependent whole-body conductance multiplier. The new tail data justify adding a peripheral observation model, but do not yet justify reinstating that fixed multiplier.

A useful causal organization is:

`E2 / reproductive state / Esr1 / clock / activity`

`→ latent tail vasomotor state v(t)`

`→ tail blood-to-surface heat transfer → Ttail and HLI`

`→ a bounded contribution to whole-body heat loss → Tcore`.

In parallel, E2, pregnancy, activity, progesterone, and circadian signals may alter heat production directly. This parallel structure matters because core-temperature changes alone cannot distinguish increased heat production from reduced heat loss.

### 7.2 Two-node physical model

A minimal mechanistic extension is

`Ccore dTcore/dt = Qmet(t) - Grest[Tcore - Ta] - Gbt[v(t)] [Tcore - Ttail]`

`Ctail dTtail/dt = Gbt[v(t)] [Tcore - Ttail] - Htail(t) [Ttail - Ta]`.

Here:

- `Grest` is non-tail whole-body conductance;
- `Gbt` is blood/body-to-tail conductance controlled by vasomotor state `v`;
- `Htail` includes tail convection, radiation, and contact loss;
- `v = 0` is constricted and `v = 1` is dilated.

Temperatures alone will not separately identify `Gbt` and `Htail` without strong priors or calibration. The model is nevertheless valuable as a conceptual observation model and as a future framework if indirect calorimetry or blood-flow data become available.

### 7.3 Statistical state model for the data available now

The more identifiable first analysis is a hierarchical model of HLI or a latent dilated/constricted state:

`logit(HLI_i,t) = α_i + f(Tcore_i,t) + βTa Ta_i,t + βphase phase_t`

`                 + βE2 E2_i,t−τ + βKO KO_i + βpreg pregnancy_i,t`

`                 + interactions + βact activity_i,t + circadian harmonics + ε_i,t`.

`α_i` is a mouse-specific random intercept. Use a bounded likelihood such as beta regression only after dealing explicitly with HLI values at or outside 0 and 1. A transformed Gaussian model with robust residuals is another defensible option.

The interactions should include at minimum:

- reproductive/E2 state × light/dark phase;
- reproductive/E2 state × ambient temperature;
- genotype × reproductive/E2 state;
- pregnancy day × phase;
- genotype × core temperature or thermal load.

If the raw series is visibly bimodal or contains abrupt switches, a two-state hidden Markov model is preferable. It can estimate:

- probability of occupying the dilated state;
- constricted→dilated and dilated→constricted transition rates;
- state-specific HLI means;
- how E2/reproductive state and *Esr1* deletion move each transition threshold.

This directly matches the literature's threshold/switching interpretation.

### 7.4 A conservative coupling back to the core model

Only after the observation model is stable should the tail state enter the core-temperature equation, for example:

`Ccore dTcore/dt = Qmet(t) - [Grest + Gtail,max v(t)] [Tcore - Ta]`.

Use the mouse-tail literature to constrain, not estimate freely, the magnitude of `Gtail,max`. A reasonable sensitivity analysis is to let maximal tail loss account for 0%, 3%, 5%, and 8% of maximal whole-body dry conductance.[9] Compare predictive performance and parameter stability. Do not assume that E2 changes total conductance by the same percentage that it changes HLI.

The tail term deserves inclusion in the final core model only if it:

1. improves held-out core-temperature prediction;
2. reduces systematic residuals around tail-dilation events;
3. yields stable parameters across mice or cohorts;
4. does not force implausible compensation in heat-production parameters.

---

## 8. Analysis plan for the PI's dataset

### Stage 1: synchronize and audit

For every mouse, retain the original time resolution and align:

- `Tcore`;
- `Ttail`;
- measured local `Ta`, if available;
- light/dark phase and zeitgeber time;
- activity;
- treatment and exact time since treatment;
- pregnancy day;
- genotype, Cre, induction treatment, and control type;
- E2 measurement or modeled E2 waveform, if available.

Record sensor position along the tail and whether it was ventral, dorsal, implanted, taped, thermographic, or telemetry-based. Rodent studies place sensors anywhere from about 1 to 6 cm from the tail base; proximal and distal temperatures are not interchangeable. Sensor attachment and bedding contact can also change local heat transfer.

Audit for:

- impossible jumps and flat lines;
- `Ttail` below/above plausible environmental/core limits;
- HLI outside 0-1;
- core-tail clock drift;
- missingness linked to movement or phase;
- pregnancy-related changes in sensor fit or tail geometry.

### Stage 2: make per-mouse summaries before inferential testing

The mouse, not the minute, is the experimental unit. Compute within each mouse and defined biological window:

- 24-hour, light-phase, and dark-phase HLI mean and median;
- circadian HLI amplitude and phase;
- 10th and 90th percentiles;
- fraction of time in a high-HLI/dilated state;
- number and duration of dilation episodes;
- constriction→dilation and dilation→constriction transition rates;
- HLI conditional on core-temperature bins and activity bins;
- tail-to-core lead/lag around spontaneous events.

Use hierarchical uncertainty or bootstrap mice, rather than treating thousands of autocorrelated measurements as independent replicates. Girbig et al. specifically showed that rat tail temperature has high-frequency fluctuations, phase-dependent variance, and artifacts from animals sitting on the tail; they modeled circadian amplitude separately from baseline trend.[18]

### Stage 3: test the biologically discriminating contrasts

1. **Cycling nonpregnant mice:** Align observations to measured E2 if available. If only cycle stage is available, do not assume stage labels perfectly reconstruct E2. Test within-mouse phase-specific HLI, amplitude, and state occupancy.

2. **Pregnancy:** Model gestational day continuously or with biologically justified windows. Test whether the HLI/core relationship shifts before the decline in core temperature. A tail-state change that precedes the core change supports altered heat-loss control; a core change without HLI change points more strongly toward heat production, other vascular beds, evaporative loss, or altered set point.

3. **Inducible *Esr1* KO:** Compare induced mutants to induction-matched littermate controls. Tamoxifen itself can alter thermoregulation, so genotype and induction treatment must not be conflated. Test genotype × reproductive/E2 state and genotype × phase effects on HLI/state transitions.

4. **Cross-cohort concordance:** Ask whether high endogenous E2 in cycling mice and pregnancy moves HLI in the same direction, and whether the movement is attenuated by *Esr1* loss. Concordance supports, but does not prove, an E2/ERα mechanism. Discordance during pregnancy is biologically plausible because pregnancy is multihormonal.

### Stage 4: identify whether the tail explains core dynamics

Fit the existing core model first without the tail state, then add lagged `v(t)` or HLI. Compare:

- held-out prediction error by mouse;
- residuals around dilation/constriction transitions;
- phase-specific bias;
- inferred E2 thermogenesis coefficient;
- whether the tail term remains within the 0-8% physiological sensitivity bound.

Because tail temperature lags local flow, test modest lags rather than forcing simultaneity. The approximately 35-second value from rat skin is a qualitative guide, not a fixed mouse-tail constant.[11] If the sampling interval is several minutes, a one-bin lag may be all the data can resolve.

### Stage 5: optional calibration for future experiments

If future data collection is possible, the greatest gains would come from:

- a controlled ambient-temperature ramp spanning roughly 22-30 °C, where tail switching is most informative;
- simultaneous indirect calorimetry to constrain whole-body conductance;
- tail blood-flow measurement or a vasodilator calibration in a subset;
- standardized imaging/sensor position and tail posture;
- a factorial comparison of E2 exposure and inducible *Esr1* status with induction-matched controls.

At a typical mouse room temperature near 22 °C, mice are below thermoneutrality and several studies show strong or apparently complete tail constriction.[9,27] The closest freely moving female E2 study nevertheless retained appreciable HLI variability at the same nominal temperature.[1] A warm ramp or mild thermal challenge would identify the switching threshold more cleanly than a single cold-housed daily mean.

---

## 9. Predictions for the three experimental conditions

These are hypotheses to test, not assumptions to encode.

### Cycling E2, nonpregnant

The rat literature predicts lower dark-phase tail temperature/HLI at higher E2, but the intact mouse study found no clear cycle-stage effect.[1,3,4] Therefore, expect a modest, phase-dependent signal that may appear more clearly in state occupancy or conditional HLI than in a 24-hour mean. Within-mouse E2 measurements are substantially more informative than categorical estrous stage.

### Pregnancy

Late pregnancy may lower core temperature and shift thermoregulatory thresholds, but the direction of the tail effect cannot be inferred from E2 alone.[13] If HLI rises before or during the core-temperature fall, increased peripheral dissipation is plausible. If HLI falls while core also falls, the dominant mechanism is more likely reduced heat production or a lower regulated set point rather than tail heat loss.

### Inducible *Esr1* KO

If the deleted ERα population mediates central vasomotor control, expect altered HLI state occupancy, phase dependence, or switching thresholds. Complete loss of all E2 effects is not guaranteed because ERβ, peripheral pathways, and developmental/induction effects can contribute.[7,8] The exact Cre line, tissue specificity, timing, and induction control determine the interpretation.

---

## 10. Decision for the current model

### Recommended now

- Calculate and model HLI from paired tail, core, and ambient temperatures.
- Treat HLI or `1-HLI` as a relative vasomotor index.
- Model E2/reproductive state as shifting the probability or threshold of dilation, with phase and ambient-temperature interactions.
- Keep the existing E2 thermogenesis pathway separable from the tail pathway.
- Use mouse-level hierarchical inference and event/state metrics.
- Test whether adding the tail state improves out-of-sample core-temperature predictions.

### Do not do yet

- Do not label HLI as percent vasoconstriction.
- Do not convert tail temperature directly to blood flow using a coefficient from another experiment.
- Do not interpret pregnancy as a dose of E2.
- Do not infer an E2-specific whole-body conductance coefficient from a tail-temperature group difference.
- Do not treat all high-frequency time points as independent observations.
- Do not hard-code the earlier 2.75% E2 conductance modification. The evidence supports a nonlinear, context-dependent state effect, and the tail is a small fraction of whole-body exchange.

### Practical Meeting 8 conclusion

Because Ivy's animals were housed at **22 °C**, interpret the tail data as a cold-housed measurement of **tail vasomotor state**, not as a calorimetric measurement of total heat loss. Strong baseline constriction is likely, and some mouse protocols find essentially complete constriction at this temperature.[9,27] However, this cannot be assumed to be constant in Ivy's freely moving females because Krajewski-Hall's mice retained appreciable HLI range at 22 °C.[1] First inspect whether Ivy's HLI traces show a genuine second/high state or reproducible excursions coupled to a later core-temperature fall. Only then estimate event occupancy, rate, duration, magnitude, and core coupling. A null mean difference between groups would reproduce the closest published mouse summary result, but it would neither prove physiological saturation nor establish equal event dynamics.

---

## 11. Metadata needed before fitting

The following details will materially change the model:

1. Exact *Esr1i* genotype, Cre/driver, targeted tissue or cell type, induction agent, dose, and timing.
2. Whether controls received the same tamoxifen/induction protocol.
3. Tail-temperature technology and exact sensor position.
4. Sampling interval and clock synchronization with core temperature.
5. Actual ambient-temperature time series at cage level.
6. Housing, bedding, cage type, airflow, single/group housing, and handling schedule.
7. Pregnancy days represented and whether data are longitudinal within mouse.
8. Estrous-stage method and any serum/plasma E2 measurements.
9. Activity, body mass, food intake, progesterone, and metabolic-rate data, if available.
10. Whether tail temperature was measured continuously or only during handling sessions.

---

## 12. Limitations of this review

- Direct, quantitative mouse E2-tail studies are sparse; much of the mechanistic threshold evidence comes from rats.
- Implant concentration is not equivalent to delivered dose, and serum levels vary over time.
- Published studies differ in ambient temperature, phase, sensor location, assay duration, and definition of conductance.
- Several results are phase-specific or nonlinear, making simple pooled means misleading.
- The pregnancy literature rarely isolates E2 from the rest of reproductive physiology.
- Human cutaneous-circulation studies help establish that estrogen can shift vascular thresholds, but their vascular bed, thermal environment, and hormone preparations do not justify direct mouse parameter transfer.

---

## Sources

1. Krajewski-Hall SJ, Blackmore EM, McMinn JR, Rance NE. “Estradiol alters body temperature regulation in the female mouse.” *Temperature* (2018). [DOI: 10.1080/23328940.2017.1384090](https://doi.org/10.1080/23328940.2017.1384090)
2. Dacks PA, Rance NE. “Effects of estradiol on the thermoneutral zone and core temperature in ovariectomized rats.” *Endocrinology* (2010). [DOI: 10.1210/en.2009-1112](https://doi.org/10.1210/en.2009-1112)
3. Williams H, Dacks PA, Rance NE. “An improved method for recording tail skin temperature in the rat reveals changes during the estrous cycle and effects of ovarian steroids.” *Endocrinology* (2010). [DOI: 10.1210/en.2010-0630](https://doi.org/10.1210/en.2010-0630)
4. Marui S, et al. “Influence of exogenous and endogenous estrogen on thermoregulatory responses to mild heat and interaction with light and dark phases.” *Journal of Physiological Sciences* (2020). [DOI: 10.1186/s12576-020-00782-x](https://doi.org/10.1186/s12576-020-00782-x)
5. Okada M, Hayashi N, Kometani M, Nakao K, Inukai T. “Influences of ovariectomy and continuous replacement of 17β-estradiol on the tail skin temperature and behavior in the forced swimming test in rats.” *Japanese Journal of Pharmacology* (1997). [DOI: 10.1254/jjp.73.93](https://doi.org/10.1254/jjp.73.93)
6. Hosono T, Chen XM, Miyatsuji A, Yoda T, Yoshida K, Yanase-Fujiwara M, Kanosue K. “Effects of estrogen on thermoregulatory tail vasomotion and heat-escape behavior in freely moving female rats.” *American Journal of Physiology-Regulatory, Integrative and Comparative Physiology* (2001). [DOI: 10.1152/ajpregu.2001.280.5.R1341](https://doi.org/10.1152/ajpregu.2001.280.5.R1341)
7. Opas EE, Rutledge SJ, Vogel RL, Rodan GA, Schmidt A. “Estrogenic control of thermoregulation in ERαKO and ERβKO mice.” *Maturitas* (2006). [DOI: 10.1016/j.maturitas.2005.04.006](https://doi.org/10.1016/j.maturitas.2005.04.006)
8. Zhang Z, et al. “Estrogen receptor alpha in the brain mediates tamoxifen-induced changes in physiology in mice.” *eLife* (2021). [DOI: 10.7554/eLife.63333](https://doi.org/10.7554/eLife.63333)
9. Škop V, et al. “The contribution of the mouse tail to thermoregulation is modest.” *American Journal of Physiology-Endocrinology and Metabolism* (2020). [DOI: 10.1152/ajpendo.00133.2020](https://doi.org/10.1152/ajpendo.00133.2020)
10. Romanovsky AA, Ivanov AI, Shimansky YP. “Ambient temperature for experiments in rats: a new method for determining the zone of thermal neutrality.” *Journal of Applied Physiology* (2002). [DOI: 10.1152/japplphysiol.01173.2001](https://doi.org/10.1152/japplphysiol.01173.2001)
11. El Bitar N, Pollin B, Karroum EG, Pinçon C, Mouraux A, Le Bars D. “Thermoregulatory vasomotor tone of the rat tail and paws in thermoneutral conditions and its impact on a behavioral model of acute pain.” *Journal of Neurophysiology* (2014). [DOI: 10.1152/jn.00721.2013](https://doi.org/10.1152/jn.00721.2013)
12. Raman ER, Vanhuyse VJ, Roberts MF. “Mathematical circulation model for the blood-flow-heat-loss relationship in the rat tail.” *Physics in Medicine and Biology* (1987). [DOI: 10.1088/0031-9155/32/7/006](https://doi.org/10.1088/0031-9155/32/7/006)
13. Eliason HL, Fewell JE. “Thermoregulatory control during pregnancy and lactation in rats.” *Journal of Applied Physiology* (1997). [DOI: 10.1152/jappl.1997.83.3.837](https://doi.org/10.1152/jappl.1997.83.3.837)
14. Zhang C, et al. “Altered thermal preference by preoptic ERα neurons in postpartum females.” *Molecular Metabolism* (2025). [DOI: 10.1016/j.molmet.2025.102108](https://doi.org/10.1016/j.molmet.2025.102108)
15. Freedman RR, Blacker CM. “Estrogen raises the sweating threshold in postmenopausal women with hot flashes.” *Fertility and Sterility* (2002). [DOI: 10.1016/S0015-0282(01)03009-6](https://doi.org/10.1016/S0015-0282(01)03009-6)
16. Brunt VE, Miner JA, Meendering JR, Kaplan PF, Minson CT. “17β-estradiol and progesterone independently augment cutaneous thermal hyperemia but not reactive hyperemia.” *Microcirculation* (2011). [DOI: 10.1111/j.1549-8719.2011.00095.x](https://doi.org/10.1111/j.1549-8719.2011.00095.x)
17. Brooks EM, Morgan AL, Pierzga JM, Wladkowski SL, O'Gorman JT, Derr JA, Kenney WL. “Chronic hormone replacement therapy alters thermoregulatory and vasomotor function in postmenopausal women.” *Journal of Applied Physiology* (1997). [DOI: 10.1152/jappl.1997.83.2.477](https://doi.org/10.1152/jappl.1997.83.2.477)
18. Girbig D, Keller K, Prelle K, Patchev V, Vonk R, Igl B-W. “A dynamic model of circadian rhythms in rodent tail skin temperature for comparison of drug effects.” *Journal of Circadian Rhythms* (2012). [DOI: 10.1186/1740-3391-10-1](https://doi.org/10.1186/1740-3391-10-1)
19. Ding F, Yao J, Zhao L, Mao Z, Chen S, Brinton RD. “Ovariectomy induces a shift in fuel availability and metabolism in the hippocampus of the female transgenic model of familial Alzheimer's.” *PLOS ONE* (2013). [DOI: 10.1371/journal.pone.0059825](https://doi.org/10.1371/journal.pone.0059825)
20. Shuto H, Yamauchi A, Ikeda M, Sohda Y, Koga A, Tominaga K, Egawa T, Kataoka Y. “Forced exercise-induced flushing of tail skin in ovariectomized mice, as a new experimental model of menopausal hot flushes.” *Journal of Pharmacological Sciences* (2005). [DOI: 10.1254/jphs.SC0050046](https://doi.org/10.1254/jphs.SC0050046)
21. Bowe J, Li XF, Kinsey-Jones J, Heyerick A, Brain S, Milligan S, O'Byrne K. “The hop phytoestrogen, 8-prenylnaringenin, reverses the ovariectomy-induced rise in skin temperature in an animal model of menopausal hot flushes.” *Journal of Endocrinology* (2006). [DOI: 10.1677/joe.1.06919](https://doi.org/10.1677/joe.1.06919)
22. Fleischer AW, Schalk JC, Wetzel EA, Hanson AM, Sem DS, Donaldson WA, Frick KM. “Long-term oral administration of a novel estrogen receptor beta agonist enhances memory and alleviates drug-induced vasodilation in young ovariectomized mice.” *Hormones and Behavior* (2021). [DOI: 10.1016/j.yhbeh.2021.104948](https://doi.org/10.1016/j.yhbeh.2021.104948)
23. Laudenslager ML, Wilkinson CW, Carlisle HJ, Hammel HT. “Energy balance in ovariectomized rats with and without estrogen replacement.” *American Journal of Physiology* (1980). [DOI: 10.1152/ajpregu.1980.238.5.R400](https://doi.org/10.1152/ajpregu.1980.238.5.R400)
24. Baker MA, Dawson DD, Peters CE, Walker AM. “Effects of estrogen on thermoregulatory evaporation in rats exposed to heat.” *American Journal of Physiology-Regulatory, Integrative and Comparative Physiology* (1994). [DOI: 10.1152/ajpregu.1994.267.3.R673](https://doi.org/10.1152/ajpregu.1994.267.3.R673)
25. Dacks PA, Krajewski SJ, Rance NE. “Ambient temperature and 17β-estradiol modify Fos immunoreactivity in the median preoptic nucleus, a putative regulator of skin vasomotion.” *Endocrinology* (2011). [DOI: 10.1210/en.2010-1230](https://doi.org/10.1210/en.2010-1230)
26. Simpkins JW. “Effects of age, reproductive status and ambient temperature on skin temperature regulation in the female rat.” *Maturitas* (1995). [DOI: 10.1016/0378-5122(94)00877-A](https://doi.org/10.1016/0378-5122(94)00877-A)
27. Pincedé I, Pollin B, Meert T, Plaghki L, Le Bars D. “Psychophysics of a nociceptive test in the mouse: ambient temperature as a key factor for variation.” *PLOS ONE* (2012). [DOI: 10.1371/journal.pone.0036699](https://doi.org/10.1371/journal.pone.0036699)
28. Young AA, Dawson NJ. “Evidence for on-off control of heat dissipation from the tail of the rat.” *Canadian Journal of Physiology and Pharmacology* (1982). [DOI: 10.1139/y82-057](https://doi.org/10.1139/y82-057)

### Local source material retained with this memo

- `source_data/Zhang2021_fig1_HLI_summary.xlsx`: published WT animal-level HLI data used for the descriptive reanalysis.
- `source_data/Zhang2021_fig6_Esr1cKO_HLI_summary.xlsx`: published conditional-KO animal-level HLI data used for the descriptive reanalysis.
- `source_data/Okada1997.pdf`: primary paper used to verify continuous minipump dose and duration.
- `source_data/PMC2954718_full.html`, `PMC3115603_full.html`, and `PMC5902204_full.html`: primary-source full text used to verify sampling intervals, analysis windows, and natural-event limitations.
- `source_data/PMC2954718_Fig2.jpg`, `PMC2954718_Fig3.jpg`, `PMC3115603_Fig2.jpg`, and `PMC5902204_Fig4.jpg`: retained primary figures used to inspect time-series and HLI summaries.
- `source_data/PMC10717522.xml`, `PMC3398330.xml`, and `PMC7924955.xml`: additional machine-readable primary-source full text retained for reproducibility.
