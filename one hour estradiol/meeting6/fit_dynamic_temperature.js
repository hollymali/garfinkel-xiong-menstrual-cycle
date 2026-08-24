"use strict";

/*
Fit the minimal dynamic temperature model

  dT/dt = b + aC*C(t) - aL*L(t) - aT*T^n/(KT^n + T^n)

to the established Meeting 6 algebraic traces.  The endocrine equations and
the upstream conservation/loss pathways match MODEL.py and
final_temp_model.py.  This script has no third-party dependencies.
*/

const MODEL_DT_MIN = 0.1;
const TOTAL_DAYS = 9;
const SAMPLE_MIN = 5;
const ANALYSIS_DAYS = 5;
const HOURS_DT = SAMPLE_MIN / 60;

const constants = {
  endocrine: {
    kG: Math.log(2) / 4,
    kL: Math.log(2) / 20,
    kE: Math.log(2) / 30,
    cL: 0.5,
    cE: 0.1145,
    K: 20.1,
    A: 6.8,
    nHalf: 0.28,
    clockPeakMin: 14 * 60,
    clockAmplitude: 0.7,
  },
  thermal: {
    e2Half: 50,
    e2Exponent: 2,
    clockPeakHours: 11.71,
    armExponent: 4,
    conservationClockWeight: 2.27,
    conservationNWeight: 0.71,
    progesteroneWeight: 0.80,
    progesteroneCoefficient: 0.0276,
    progesteroneHalf: 10,
    progesteroneExponent: 2,
    lossClockWeight: 0.71,
    directE2Weight: 0.02,
    e2Coefficient: 0.000088,
    linearHormones: false,
    targetBase: 35.94,
    targetGain: 2.53,
  },
  conditions: {
    cyclingP: 0.2,
    lowE2: 5,
    lowP: 0.2,
    pregnancyE2: 215,
    pregnancyP: 25,
  },
};

function endocrineClock(timeMin) {
  const tod = ((timeMin % 1440) + 1440) % 1440;
  return Math.cos(2 * Math.PI * (tod - constants.endocrine.clockPeakMin) / 1440);
}

function simulateEndocrine() {
  const e = constants.endocrine;
  const steps = Math.round(TOTAL_DAYS * 1440 / MODEL_DT_MIN);
  const sampleStride = Math.round(SAMPLE_MIN / MODEL_DT_MIN);
  const samples = Math.floor(steps / sampleStride) + 1;
  const nOut = new Float64Array(samples);
  const e2Out = new Float64Array(samples);

  let N = 1 / (1 + Math.pow(15 / e.K, 5));
  let gnrh = 0;
  let lh = 5;
  let e2 = 15;
  let phase = 0;
  let sampleIndex = 0;
  nOut[0] = N;
  e2Out[0] = e2;

  for (let i = 0; i < steps; i += 1) {
    const time = i * MODEL_DT_MIN;
    const n2 = N * N;
    const n4 = n2 * n2;
    const ampDrive = n4 / (Math.pow(e.nHalf, 4) + n4);
    const freqDrive = n2 / (Math.pow(e.nHalf, 2) + n2);
    const period = 63 + (114 - 63) * freqDrive;
    const newPhase = phase + MODEL_DT_MIN / period;
    const fire = Math.floor(newPhase) > Math.floor(phase);
    const amplitudeMultiplier = 1 + e.clockAmplitude * (1 + endocrineClock(time)) / 2;
    const pulse = fire ? e.A * ampDrive * amplitudeMultiplier : 0;
    phase = newPhase;

    const nextN = N + MODEL_DT_MIN * (0.25 / (1 + Math.pow(e2 / e.K, 5)) - 0.25 * N);
    const nextGnrh = gnrh + MODEL_DT_MIN * (-e.kG * gnrh) + pulse;
    const nextLh = lh + MODEL_DT_MIN * (e.cL * gnrh - e.kL * lh);
    const nextE2 = e2 + MODEL_DT_MIN * (e.cE * lh - e.kE * e2);
    N = nextN;
    gnrh = nextGnrh;
    lh = nextLh;
    e2 = nextE2;

    if ((i + 1) % sampleStride === 0) {
      sampleIndex += 1;
      nOut[sampleIndex] = N;
      e2Out[sampleIndex] = e2;
    }
  }
  return { N: nOut, E2: e2Out };
}

function hill(value, exponent, half = 1) {
  const x = Math.max(value, 0) / half;
  const power = Math.pow(x, exponent);
  return power / (1 + power);
}

function thermalClock(hour, thermal = constants.thermal) {
  const t = thermal;
  const hod = ((hour % 24) + 24) % 24;
  return 0.5 * (1 + Math.cos(2 * Math.PI * (hod - t.clockPeakHours) / 24));
}

function buildCondition(N, E2, progesterone, thermalOverrides = {}) {
  const t = { ...constants.thermal, ...thermalOverrides };
  const count = E2.length;
  const conservation = new Float64Array(count);
  const loss = new Float64Array(count);
  const target = new Float64Array(count);
  const pArm = hill(progesterone, t.progesteroneExponent, t.progesteroneHalf);
  const progesteroneInput = t.linearHormones
    ? t.progesteroneCoefficient * progesterone
    : t.progesteroneWeight * pArm;
  for (let i = 0; i < count; i += 1) {
    const clock = thermalClock(i * HOURS_DT, t);
    const direct = hill(E2[i], t.e2Exponent, t.e2Half);
    conservation[i] = hill(
      t.conservationClockWeight * clock
        + t.conservationNWeight * N[i]
        + progesteroneInput,
      t.armExponent,
    );
    loss[i] = hill(
      t.lossClockWeight * clock
        + (t.linearHormones ? t.e2Coefficient * E2[i] : t.directE2Weight * direct),
      t.armExponent,
    );
    target[i] = t.targetBase + t.targetGain * (conservation[i] - loss[i]);
  }
  return { conservation, loss, target };
}

function constantArray(length, value) {
  const output = new Float64Array(length);
  output.fill(value);
  return output;
}

function buildInputs(thermalOverrides = {}) {
  const trace = simulateEndocrine();
  const n = trace.N.length;
  const c = constants.conditions;
  const retainedPoints = Math.round(ANALYSIS_DAYS * 24 / HOURS_DT) + 1;
  const retainedE2 = trace.E2.subarray(trace.E2.length - retainedPoints);
  let lowE2 = Infinity;
  for (const value of retainedE2) lowE2 = Math.min(lowE2, value);
  const lowN = 1 / (1 + Math.pow(lowE2 / constants.endocrine.K, 5));
  const pregnancyN = 1 / (1 + Math.pow(c.pregnancyE2 / constants.endocrine.K, 5));
  return {
    cycling: buildCondition(trace.N, trace.E2, c.cyclingP, thermalOverrides),
    low: buildCondition(
      constantArray(n, lowN),
      constantArray(n, lowE2),
      c.lowP,
      thermalOverrides,
    ),
    pregnancy: buildCondition(
      constantArray(n, pregnancyN),
      constantArray(n, c.pregnancyE2),
      c.pregnancyP,
      thermalOverrides,
    ),
  };
}

function selfHill(T, K, exponent) {
  const ratio = Math.max(T, 1e-9) / K;
  const power = Math.pow(ratio, exponent);
  return power / (1 + power);
}

function rhs(T, conservation, loss, p) {
  return p.b + p.aC * conservation - p.aL * loss
    - p.aT * selfHill(T, p.KT, p.n);
}

function simulateDynamic(condition, p) {
  const n = condition.target.length;
  const output = new Float64Array(n);
  let T = 37.0;
  output[0] = T;
  for (let i = 0; i < n - 1; i += 1) {
    const c0 = condition.conservation[i];
    const c1 = condition.conservation[i + 1];
    const l0 = condition.loss[i];
    const l1 = condition.loss[i + 1];
    const cm = 0.5 * (c0 + c1);
    const lm = 0.5 * (l0 + l1);
    const k1 = rhs(T, c0, l0, p);
    const k2 = rhs(T + 0.5 * HOURS_DT * k1, cm, lm, p);
    const k3 = rhs(T + 0.5 * HOURS_DT * k2, cm, lm, p);
    const k4 = rhs(T + HOURS_DT * k3, c1, l1, p);
    T += HOURS_DT * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    if (!Number.isFinite(T) || T < 20 || T > 50) return null;
    output[i + 1] = T;
  }
  return output;
}

function dailyMetrics(values, start) {
  const perDay = Math.round(24 / HOURS_DT);
  let peakSum = 0;
  let troughSum = 0;
  for (let day = 0; day < ANALYSIS_DAYS; day += 1) {
    const first = start + day * perDay;
    const last = Math.min(first + perDay, values.length - 1);
    let peak = -Infinity;
    let trough = Infinity;
    for (let i = first; i <= last; i += 1) {
      peak = Math.max(peak, values[i]);
      trough = Math.min(trough, values[i]);
    }
    peakSum += peak;
    troughSum += trough;
  }
  return {
    peak: peakSum / ANALYSIS_DAYS,
    trough: troughSum / ANALYSIS_DAYS,
    amplitude: (peakSum - troughSum) / ANALYSIS_DAYS,
  };
}

function score(inputs, p, includeSeries = false) {
  const start = Math.round((TOTAL_DAYS - ANALYSIS_DAYS) * 24 / HOURS_DT);
  const names = ["cycling", "low", "pregnancy"];
  const conditionWeights = { cycling: 3, low: 1, pregnancy: 1 };
  let totalMse = 0;
  let totalWeight = 0;
  const outputs = {};
  const metrics = {};

  for (const name of names) {
    const actual = simulateDynamic(inputs[name], p);
    if (actual === null) return { value: 1e12 };
    let mse = 0;
    let count = 0;
    for (let i = start; i < actual.length; i += 1) {
      const error = actual[i] - inputs[name].target[i];
      mse += error * error;
      count += 1;
    }
    mse /= count;
    totalMse += conditionWeights[name] * mse;
    totalWeight += conditionWeights[name];

    const got = dailyMetrics(actual, start);
    const wanted = dailyMetrics(inputs[name].target, start);
    metrics[name] = { got, wanted, rmse: Math.sqrt(mse) };
    const peakError = got.peak - wanted.peak;
    const troughError = got.trough - wanted.trough;
    const amplitudeError = got.amplitude - wanted.amplitude;
    totalMse += 2.5 * peakError * peakError
      + 2.5 * troughError * troughError
      + amplitudeError * amplitudeError;
    totalWeight += 6;
    if (includeSeries) outputs[name] = actual;
  }

  // Preserve the defining phenotype: cycling troughs lower than both constant
  // conditions while peak heights remain aligned.
  const desiredLowGap = metrics.low.wanted.trough - metrics.cycling.wanted.trough;
  const gotLowGap = metrics.low.got.trough - metrics.cycling.got.trough;
  const desiredPregGap = metrics.pregnancy.wanted.trough - metrics.cycling.wanted.trough;
  const gotPregGap = metrics.pregnancy.got.trough - metrics.cycling.got.trough;
  const gapPenalty = Math.pow(gotLowGap - desiredLowGap, 2)
    + Math.pow(gotPregGap - desiredPregGap, 2);
  const peakSpread = Math.max(
    metrics.cycling.got.peak,
    metrics.low.got.peak,
    metrics.pregnancy.got.peak,
  ) - Math.min(
    metrics.cycling.got.peak,
    metrics.low.got.peak,
    metrics.pregnancy.got.peak,
  );
  totalMse += 4 * gapPenalty + 3 * peakSpread * peakSpread;
  totalWeight += 11;

  return {
    value: totalMse / totalWeight,
    outputs,
    metrics,
    start,
  };
}

function makeRng(seed) {
  let state = seed >>> 0;
  return function random() {
    state += 0x6D2B79F5;
    let t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const bounds = [
  [0, 40],    // basal heat input b
  [0, 25],    // conservation gain aC
  [0, 25],    // hormone/clock loss gain aL
  [1, 100],   // maximum self-loss aT
  [30, 44],   // self-loss Hill half saturation KT
];

function vectorToParams(vector, exponent) {
  return {
    b: vector[0],
    aC: vector[1],
    aL: vector[2],
    aT: vector[3],
    KT: vector[4],
    n: exponent,
  };
}

function fitExponent(inputs, exponent, seed) {
  const random = makeRng(seed);
  const dimension = bounds.length;
  const populationSize = 54;
  const generations = 240;
  const F = 0.72;
  const CR = 0.88;
  const pop = [];
  const values = [];

  for (let i = 0; i < populationSize; i += 1) {
    const v = bounds.map(([lo, hi]) => lo + random() * (hi - lo));
    pop.push(v);
    values.push(score(inputs, vectorToParams(v, exponent)).value);
  }

  for (let generation = 0; generation < generations; generation += 1) {
    for (let i = 0; i < populationSize; i += 1) {
      let a;
      let b;
      let c;
      do { a = Math.floor(random() * populationSize); } while (a === i);
      do { b = Math.floor(random() * populationSize); } while (b === i || b === a);
      do { c = Math.floor(random() * populationSize); } while (c === i || c === a || c === b);
      const forced = Math.floor(random() * dimension);
      const trial = new Array(dimension);
      for (let j = 0; j < dimension; j += 1) {
        let value = (random() < CR || j === forced)
          ? pop[a][j] + F * (pop[b][j] - pop[c][j])
          : pop[i][j];
        const [lo, hi] = bounds[j];
        if (value < lo || value > hi) value = lo + random() * (hi - lo);
        trial[j] = value;
      }
      const trialValue = score(inputs, vectorToParams(trial, exponent)).value;
      if (trialValue < values[i]) {
        pop[i] = trial;
        values[i] = trialValue;
      }
    }
  }

  let bestIndex = 0;
  for (let i = 1; i < populationSize; i += 1) {
    if (values[i] < values[bestIndex]) bestIndex = i;
  }
  return {
    params: vectorToParams(pop[bestIndex], exponent),
    score: values[bestIndex],
  };
}

function foldSeries(values, start) {
  const perDay = Math.round(24 / HOURS_DT);
  const folded = [];
  for (let index = 0; index <= perDay; index += 1) {
    let sum = 0;
    let count = 0;
    for (let day = 0; day < ANALYSIS_DAYS; day += 1) {
      const i = start + day * perDay + index;
      if (i < values.length) {
        sum += values[i];
        count += 1;
      }
    }
    folded.push(sum / count);
  }
  return folded;
}

function main() {
  const inputs = buildInputs();
  const candidates = [];
  for (let exponent = 1; exponent <= 5; exponent += 1) {
    const fit = fitExponent(inputs, exponent, 20260824 + exponent * 7919);
    candidates.push(fit);
    process.stderr.write(`n=${exponent} score=${fit.score.toFixed(7)}\n`);
  }
  candidates.sort((a, b) => a.score - b.score);
  const best = candidates[0];
  const result = score(inputs, best.params, true);
  const output = {
    model: "dT/dt = b + aC*Conservation - aL*HormonalLoss - aT*T^n/(KT^n + T^n)",
    units: {
      b_aC_aL_aT: "degC/hour after whole-body heat-capacity scaling",
      KT: "degC in the fitted temperature coordinate",
      n: "dimensionless",
    },
    best_parameters: best.params,
    objective: result.value,
    exponent_candidates: candidates.map((candidate) => ({
      n: candidate.params.n,
      score: candidate.score,
      parameters: candidate.params,
    })),
    metrics: result.metrics,
    folded: {},
  };
  for (const name of ["cycling", "low", "pregnancy"]) {
    output.folded[name] = {
      hour: Array.from({ length: Math.round(24 / HOURS_DT) + 1 }, (_, i) => i * HOURS_DT),
      dynamic: foldSeries(result.outputs[name], result.start),
      target: foldSeries(inputs[name].target, result.start),
    };
  }
  if (process.argv.includes("--summary")) delete output.folded;
  process.stdout.write(`${JSON.stringify(output)}\n`);
}

if (require.main === module) main();

module.exports = {
  buildInputs,
  buildCondition,
  constantArray,
  constants,
  dailyMetrics,
  score,
  simulateEndocrine,
  simulateDynamic,
};
