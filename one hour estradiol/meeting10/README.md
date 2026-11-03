# Meeting 10: first-edition sleep and SCN modulation

This edition adds two distinct time-of-day effects to the HPG model:

1. sleep slows the GnRH pulse generator by multiplying its instantaneous
   period by `1.57` during the Licinio sleep window; and
2. a smooth SCN-associated gain between `1` and `2` multiplies GnRH pulse
   amplitude and peaks at 08:00.

The two drives are intentionally separate. Sleep changes **when** a pulse
occurs; the SCN hypothesis changes the **jump size** only when a pulse occurs.

## Clock convention

The shared project clock uses ZT0 = 07:00. Therefore:

- Licinio sleep from 23:00 to 07:00 is ZT16 to ZT24;
- the provisional 08:00 SCN-gain peak is ZT1.

For phase (z(t)=t\bmod 24), define the sleep-period multiplier

$$
M_{\mathrm{sleep}}(t)=
\begin{cases}
1.57, & 16\le z(t)<24,\\
1, & 0\le z(t)<16.
\end{cases}
$$

This is implemented as an editable bout list, in the same style as the
activity bouts.

## Pulse timing

Let

$$
H_+(x;K,n)=\frac{x^n}{K^n+x^n},
\qquad
H_-(x;K,n)=\frac{K^n}{K^n+x^n}.
$$

The neuropeptide state remains

$$
\frac{dN}{dt}=100H_-(E_2;11,5)-N.
$$

The waking pulse period and sleep-modulated period, in minutes, are

$$
P_0(N)=75-(75-51)H_+(N;1,5),
$$

$$
P(N,t)=P_0(N)M_{\mathrm{sleep}}(t).
$$

The continuous pulse counter evolves as

$$
\frac{d\phi}{dt}=\frac{60}{P(N,t)},
$$

and a pulse fires whenever \(\phi\) crosses an integer. Multiplication by
`1.57` is the rounded ratio of Licinio's mean nighttime to daytime LH
inter-pulse intervals, (118/75=1.573\). Licinio's participants were in bed,
in darkness, and asleep from 23:00 to 07:00. Their tabulated comparison used
23:00--08:00 as the nighttime analysis window; this model uses the actual
23:00--07:00 sleep interval.

## SCN-associated amplitude gain

The first-edition smooth gain is

$$
M_{\mathrm{SCN}}(t)
=1.5+0.5\cos\!\left(\frac{2\pi[z(t)-1]}{24}\right).
$$

It has

$$
M_{\mathrm{SCN}}(\mathrm{ZT1})=2,
\qquad
M_{\mathrm{SCN}}(\mathrm{ZT13})=1.
$$

At pulse time (t_j), the GnRH jump is

$$
\mathrm{GnRH}(t_j^+)
=\mathrm{GnRH}(t_j^-)
+100H_+(N(t_j);1,5)M_{\mathrm{SCN}}(t_j).
$$

Thus the SCN factor does not create pulses or change their timing.

## Hormone cascade

Between pulses, the existing cascade is unchanged:

$$
\begin{aligned}
\frac{d\mathrm{GnRH}}{dt}&=-k_{\mathrm{GnRH}}\mathrm{GnRH},\\
\frac{d\mathrm{LH}}{dt}&=\mathrm{GnRH}-k_{\mathrm{LH}}\mathrm{LH},\\
\frac{dE_2}{dt}&=15\,\mathrm{LH}-k_{E_2}E_2,
\end{aligned}
$$

with half-lives of 4, 20, and 22 minutes, respectively.

## Why the provisional peak is 08:00 / ZT1

There is no human experiment that measured a twofold circadian SCN gain on
GnRH and located its maximum. The 08:00 phase is therefore a transparent
first-edition reproductive-gating hypothesis. Licinio explicitly exposed
participants to light from 07:00 to 23:00, making 08:00 ZT1 in this model.
Kerdelhué sampled every four hours and detected LH-surge initiation together
with the cortisol maximum at 04:00 in 20% and 08:00 in 80% of participants.
That study did not report a standardized lights-on time, so its 08:00 result
supports early-morning civil-clock timing but cannot independently establish
ZT1. The four-hour sampling grid also means that an event detected at 08:00
could have begun during the preceding 04:00--08:00 interval.

The gain is interpreted as an SCN-linked reproductive gate rather than total
SCN firing. Its maximum can coincide with surge initiation even though plasma
LH continues rising afterward: human midcycle data show an approximately
14-hour ascending limb, followed by an approximately 14-hour peak plateau.
Estradiol-dependent hypothalamic and pituitary dynamics can therefore sustain
the surge after the daily SCN gain has started to fall.

Direct mechanistic evidence comes from rodents and has a different phase.
In estradiol-treated mice, VIP excited GnRH neurons near surge onset at about
ZT10; at surge peak, blocking VIP receptors reduced GnRH-neuron firing,
consistent with an endogenous VIP contribution. In intact mice on a 12:12
light-dark schedule, a later study measured mean LH-surge onset near ZT11.9
and peak near ZT13. These data support a smooth circadian gain near the surge,
but not a human nighttime maximum or the numerical range `1` to `2`.

Nor should this gain be read as total human SCN neural activity. Human
perfusion imaging found the highest SCN activity at noon and the lowest at
06:00. The modeled early-morning gain is therefore an output-specific
reproductive gate attributed to the SCN pathway, not a waveform for whole-SCN
activity.

Consequently, `SCN_ENDOCRINE_PEAK_ZT = 1.0` is a phase parameter intended for
sensitivity testing.

## Sources

- [Licinio et al. (1998): 24-hour LH and estradiol sampling](https://www.pnas.org/doi/10.1073/pnas.95.5.2541)
- [Cahill et al. (1998): human LH-surge onset is usually 00:00--08:00](https://pubmed.ncbi.nlm.nih.gov/9660421/)
- [Kerdelhué et al. (2002): early-morning LH-surge onset and cortisol coupling](https://pubmed.ncbi.nlm.nih.gov/11914587/)
- [Hoff et al. (1983): approximately 14-hour ascending limb of the human LH surge](https://pubmed.ncbi.nlm.nih.gov/6411753/)
- [Christian and Moenter (2008): time-gated VIP effects on GnRH neurons](https://pubmed.ncbi.nlm.nih.gov/18326000/)
- [Bahougne et al. (2020): mouse LH-surge onset and peak in ZT](https://pubmed.ncbi.nlm.nih.gov/32100021/)
- [Oka et al. (2024): human SCN activity is highest at noon and lowest near 06:00](https://pubmed.ncbi.nlm.nih.gov/38238074/)
