# Weight-only temperature-model exploration memory

## Scope and non-negotiable constraints

This exploration changed only coefficient weights in copies of the model. The original `MODEL.py` was not edited. Model equations, pathway topology, phase timing, six-second integration, warm-up procedure, experimental conditions, and the cycling-E2 input waveform were held fixed.

The desired temperature behavior was:

- preserve legitimate bimodality in cycling E2 and pregnancy;
- make the first cycling lobe clearly taller than the second;
- keep cycling troughs low with little trough ripple;
- make constant-low-E2 and constant-high-E2 oscillation amplitudes about half the cycling-E2 amplitude;
- retain close visual alignment with `references/mousebodytemp.png`;
- never improve the temperature result by changing the E2 time series.

## Coefficient groups explored

The searches varied only the ten existing weights below, singly and in groups:

- baseline heat production (`B`)
- circadian gain (`CG`)
- heat-loss gain (`HLG`)
- clock/conservation coupling (`CC`)
- conservation Hill exponent (`N`)
- progesterone weight (`P`)
- loss/conservation coupling (`LC`)
- E2 loss weight (`E2W`)
- medium self-loss (`SLM`)
- high self-loss (`SLH`)

## Main candidates and outcomes

### Balanced reference-aligned candidate

Weights:

| Weight | Value |
|---|---:|
| B | 4.287827872594785 |
| CG | 1.4101630616719283 |
| HLG | 1.4177032021540943 |
| CC | 0.64 |
| N | 4.6 |
| P | 0.033 |
| LC | 1.30 |
| E2W | 1.0000044169174835e-6 |
| SLM | 356.4444972841674 |
| SLH | 87.99998323181106 |

Key metrics:

- cycling amplitude: 2.2496
- constant-low/cycling amplitude ratio: 0.6585
- pregnancy/cycling amplitude ratio: 0.3720
- cycling trough: 36.299
- trough ripple: 0.3365
- strict first-minus-second peak gap: only 0.123
- cycling and pregnancy passed the true-bimodality test

This was the best balanced/reference-like family, but it did not make the first lobe sufficiently taller than the second.

### Dominant-first-lobe, high-N candidate

Weights:

| Weight | Value |
|---|---:|
| B | 73.72380651341913 |
| CG | 21.995415552140607 |
| HLG | 21.01695558624496 |
| CC | 2.2011135676269085 |
| N | 3.883984140661388 |
| P | 0.0068244890394569425 |
| LC | 3.277440246607972 |
| E2W | 0.001794464307916788 |
| SLM | 772.2995775423863 |
| SLH | 56.986388319861135 |

Key metrics:

- cycling amplitude: 2.7868
- constant-low/cycling amplitude ratio: 0.7399
- pregnancy/cycling amplitude ratio: 0.4082
- first and second cycling peaks: 38.416 and 37.216 (gap 1.20)
- cycling peak prominences: 0.259 and 0.817
- post-peak trough: 36.399
- trough ripple: 0.272
- cycling and pregnancy passed the true-bimodality test
- the E2 waveform matched the reference model exactly

This produced a genuinely taller first lobe while retaining strict bimodality, but its temperature trace—especially pregnancy and constant-low E2—was not visually convincing enough to keep.

## Other informative endpoints

- An older restored-bimodality endpoint had a 1.236 first-lobe gap, constant-low ratio near 0.55, and pregnancy ratio near 0.46, but trough ripple was about 0.872 and therefore unacceptable.
- Several smooth/roughly-half-amplitude candidates reached constant-low ratios near 0.55 and trough ripple near 0.2. Their apparent second peak was only a shoulder under the strict prominence test, so they were not legitimately bimodal.
- A v2 dominant candidate retained true bimodality, a 0.472 first-lobe gap, and ripple of 0.163, but its constant-low ratio was about 0.994 rather than one-half.
- A current-family grouped search reached cycling amplitude 2.787, constant-low ratio 0.693, a 0.39 first-lobe gap, and ripple 0.446; it missed both dominance and smooth-trough requirements.
- One endpoint passed the numerical hard gates with a 1.98 first-lobe gap, constant-low ratio 0.68, pregnancy ratio 0.322, ripple 0.33, and true bimodality in both conditions. Its cycling amplitude rose to about 3.99 and the waveform looked poor, so it was rejected.
- Scaling gain and damping around a dominant endpoint yielded cycling amplitude 2.787, constant-low ratio 0.704, ripple 0.220, and a 1.412 nominal gap, but the second-lobe prominence fell below the strict bimodality threshold.

## What the searches taught us

- The coefficient-only model has a strong Pareto tradeoff: first-lobe dominance plus strict bimodality and smooth troughs compete with half-amplitude constant conditions and close reference alignment.
- Raising the conservation exponent `N` alone did not solve the problem. It commonly increased the constant/cycling amplitude ratio or collapsed the first/second-lobe gap.
- Increasing E2 loss did not smooth the balanced bridge as hoped and tended to weaken pregnancy bimodality.
- Interpolating whole coefficient groups often crossed an unstable operating region, producing unrealistic temperatures or trough ripple.
- Co-scaling gains with self-damping could resize total amplitude, but it did not substantially improve the constant/cycling amplitude ratio and often weakened the second lobe.
- Every verified candidate kept the cycling-E2 time series unchanged; temperature differences came only from coefficient weights.

## Pivot note

The searches suggest that satisfying all waveform requirements simultaneously may require a structural mechanism—such as phase-specific or condition-dependent damping, an additional delayed pathway, or another state—rather than more coefficient-only tuning. That is a future design decision, not a change made here.

All temporary model copies, search scripts, metrics, and generated comparison figures created during this exploration outside `meeting6/POWERPOINT` were removed on 2026-08-27. The original model/reference files and everything under `meeting6/POWERPOINT` were preserved.
