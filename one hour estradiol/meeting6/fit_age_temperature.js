"use strict";

/*
Fit simple temperature-balance weights after replacing the hormone-specific
progesterone and E2/MPA Hill transforms with one linear coefficient each.
The outer pathway Hill functions, the self-loss Hill function
(K=39.08 C, n=5), and the endocrine model remain fixed.
*/

const fs = require("fs");
const path = require("path");
const base = require("./fit_dynamic_temperature.js");

const DT_H = 5 / 60;
const TOTAL_DAYS = 9;
const ANALYSIS_DAYS = 5;
const START = Math.round((TOTAL_DAYS - ANALYSIS_DAYS) * 24 / DT_H);
const CURRENT = {
  b: 49.42,
  aC: 9.76,
  aL: 15.02,
  aT: 125.0,
  KT: 39.08,
  n: 5,
  clockPeakHours: 18.16,
  conservationNWeight: 0.61,
  progesteroneCoefficient: 0.0276,
  e2Coefficient: 0.000088,
};

function parseMouseData() {
  const csvPath = path.join(__dirname, "mice_age_temperature_digitized.csv");
  const rows = fs.readFileSync(csvPath, "utf8").trim().split(/\r?\n/).slice(1);
  return rows.map((row) => {
    const values = row.split(",").map(value => Number(value.replaceAll('"', "")));
    return { day: values[0], young: values[1], aged: values[2] };
  });
}

function conditionInputs() {
  const trace = base.simulateEndocrine();
  const retainedPoints = Math.round(ANALYSIS_DAYS * 24 / DT_H) + 1;
  const retainedE2 = trace.E2.subarray(trace.E2.length - retainedPoints);
  let lowE2 = Infinity;
  for (const value of retainedE2) lowE2 = Math.min(lowE2, value);
  const endocrineK = base.constants.endocrine.K;
  const lowN = 1 / (1 + Math.pow(lowE2 / endocrineK, 5));
  const pregnancyE2 = 215;
  const pregnancyN = 1 / (1 + Math.pow(pregnancyE2 / endocrineK, 5));
  return {
    cycling: { N: trace.N, E2: trace.E2, progesterone: 0.2 },
    low: {
      N: base.constantArray(trace.N.length, lowN),
      E2: base.constantArray(trace.N.length, lowE2),
      progesterone: 0.2,
    },
    pregnancy: {
      N: base.constantArray(trace.N.length, pregnancyN),
      E2: base.constantArray(trace.N.length, pregnancyE2),
      progesterone: 25,
    },
    lowE2,
  };
}

function upstreamFor(
  raw,
  clockPeakHours,
  conservationNWeight,
  progesteroneCoefficient,
  e2Coefficient,
  linearHormones,
) {
  const override = {
    clockPeakHours,
    conservationNWeight,
    lossClockWeight: 0.78,
    progesteroneCoefficient,
    e2Coefficient,
    linearHormones,
  };
  return {
    cycling: base.buildCondition(
      raw.cycling.N,
      raw.cycling.E2,
      raw.cycling.progesterone,
      override,
    ),
    low: base.buildCondition(
      raw.low.N,
      raw.low.E2,
      raw.low.progesterone,
      override,
    ),
    pregnancy: base.buildCondition(
      raw.pregnancy.N,
      raw.pregnancy.E2,
      raw.pregnancy.progesterone,
      override,
    ),
  };
}

function interpolateModel(values, day) {
  const position = START + day * 24 / DT_H;
  const lower = Math.max(START, Math.min(values.length - 1, Math.floor(position)));
  const upper = Math.max(START, Math.min(values.length - 1, lower + 1));
  const fraction = position - lower;
  return values[lower] + fraction * (values[upper] - values[lower]);
}

function traceMetrics(values) {
  let peak = -Infinity;
  let trough = Infinity;
  for (let i = START; i < values.length; i += 1) {
    peak = Math.max(peak, values[i]);
    trough = Math.min(trough, values[i]);
  }
  return { peak, trough, amplitude: peak - trough };
}

function mouseMetrics(mouse, key) {
  let peak = -Infinity;
  let trough = Infinity;
  for (const point of mouse) {
    peak = Math.max(peak, point[key]);
    trough = Math.min(trough, point[key]);
  }
  return { peak, trough, amplitude: peak - trough };
}

function objective(
  mouse,
  raw,
  vector,
  includeSeries = false,
  linearHormones = true,
) {
  const p = {
    b: vector[0],
    aC: vector[1],
    aL: vector[2],
    aT: CURRENT.aT,
    KT: CURRENT.KT,
    n: CURRENT.n,
  };
  const clockPeakHours = vector[3];
  const conservationNWeight = vector[4];
  const progesteroneCoefficient = vector[5];
  const e2Coefficient = vector[6];
  const upstream = upstreamFor(
    raw,
    clockPeakHours,
    conservationNWeight,
    progesteroneCoefficient,
    e2Coefficient,
    linearHormones,
  );
  const cycling = base.simulateDynamic(upstream.cycling, p);
  const low = base.simulateDynamic(upstream.low, p);
  const pregnancy = base.simulateDynamic(upstream.pregnancy, p);
  if (cycling === null || low === null || pregnancy === null) {
    return { value: 1e12 };
  }

  let cyclingMse = 0;
  let lowMse = 0;
  for (const point of mouse) {
    cyclingMse += Math.pow(interpolateModel(cycling, point.day) - point.young, 2);
    lowMse += Math.pow(interpolateModel(low, point.day) - point.aged, 2);
  }
  cyclingMse /= mouse.length;
  lowMse /= mouse.length;

  const cyclingMetrics = traceMetrics(cycling);
  const lowMetrics = traceMetrics(low);
  const pregnancyMetrics = traceMetrics(pregnancy);
  const wantedCycling = mouseMetrics(mouse, "young");
  const wantedLow = mouseMetrics(mouse, "aged");
  // Keep the defining lower cycling nadir and similar peak heights.
  const troughGap = lowMetrics.trough - cyclingMetrics.trough;
  const peakGap = lowMetrics.peak - cyclingMetrics.peak;
  let value = cyclingMse + lowMse;
  value += 4.0 * Math.pow(cyclingMetrics.peak - wantedCycling.peak, 2);
  value += 4.0 * Math.pow(cyclingMetrics.trough - wantedCycling.trough, 2);
  value += 4.0 * Math.pow(lowMetrics.peak - wantedLow.peak, 2);
  value += 4.0 * Math.pow(lowMetrics.trough - wantedLow.trough, 2);
  value += 1.0 * Math.pow(troughGap - 0.32, 2);
  value += 1.0 * Math.pow(peakGap + 0.32, 2);
  value += 3.0 * Math.pow(pregnancyMetrics.peak - 37.96, 2);
  value += 3.0 * Math.pow(pregnancyMetrics.trough - 36.33, 2);
  value += 2.0 * Math.pow(pregnancyMetrics.peak - cyclingMetrics.peak, 2);

  // A light regularizer keeps the requested adjustment genuinely simple.
  value += 0.015 * Math.pow((p.b - CURRENT.b) / 4.0, 2);
  value += 0.015 * Math.pow((p.aC - CURRENT.aC) / 3.0, 2);
  value += 0.015 * Math.pow((p.aL - CURRENT.aL) / 3.0, 2);
  value += 0.015 * Math.pow(
    (conservationNWeight - CURRENT.conservationNWeight) / 0.25,
    2,
  );
  value += 0.005 * Math.pow(
    (progesteroneCoefficient - CURRENT.progesteroneCoefficient) / 0.02,
    2,
  );
  value += 0.005 * Math.pow(
    (e2Coefficient - CURRENT.e2Coefficient) / 0.0001,
    2,
  );

  return {
    value,
    p,
    clockPeakHours,
    conservationNWeight,
    progesteroneCoefficient,
    e2Coefficient,
    cyclingMse,
    lowMse,
    metrics: {
      cycling: cyclingMetrics,
      low: lowMetrics,
      pregnancy: pregnancyMetrics,
    },
    ...(includeSeries ? { cycling, low, pregnancy } : {}),
  };
}

function rng(seed) {
  let state = seed >>> 0;
  return () => {
    state += 0x6D2B79F5;
    let t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const BOUNDS = [
  [47.0, 52.0],       // basal input b
  [7.0, 13.0],        // conservation coefficient aC
  [11.0, 18.0],       // hormonal-loss coefficient aL
  [15.0, 21.0],       // shared thermal-clock peak hour
  [0.35, 0.85],       // N contribution to conservation
  [0.0, 0.055],       // one raw-progesterone coefficient
  [0.0, 0.00030],     // one raw-E2 coefficient
];

function fit(mouse, raw) {
  const random = rng(20260824);
  const populationSize = 40;
  const generations = 120;
  const dimension = BOUNDS.length;
  const population = [];
  const scores = [];

  population.push([
    CURRENT.b,
    CURRENT.aC,
    CURRENT.aL,
    CURRENT.clockPeakHours,
    CURRENT.conservationNWeight,
    CURRENT.progesteroneCoefficient,
    CURRENT.e2Coefficient,
  ]);
  scores.push(objective(mouse, raw, population[0]).value);
  for (let i = 1; i < populationSize; i += 1) {
    const vector = BOUNDS.map(([lo, hi]) => lo + random() * (hi - lo));
    population.push(vector);
    scores.push(objective(mouse, raw, vector).value);
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
      const trial = [];
      for (let j = 0; j < dimension; j += 1) {
        let candidate = (random() < 0.88 || j === forced)
          ? population[a][j] + 0.70 * (population[b][j] - population[c][j])
          : population[i][j];
        const [lo, hi] = BOUNDS[j];
        if (candidate < lo || candidate > hi) candidate = lo + random() * (hi - lo);
        trial.push(candidate);
      }
      const trialScore = objective(mouse, raw, trial).value;
      if (trialScore < scores[i]) {
        population[i] = trial;
        scores[i] = trialScore;
      }
    }
    if ((generation + 1) % 20 === 0) {
      process.stderr.write(`generation ${generation + 1}\n`);
    }
  }

  let best = 0;
  for (let i = 1; i < populationSize; i += 1) {
    if (scores[i] < scores[best]) best = i;
  }
  return objective(mouse, raw, population[best], true);
}

function sampledSeries(values) {
  const output = [];
  for (let i = START; i < values.length; i += 3) {
    output.push({
      day: (i - START) * DT_H / 24,
      temperature_c: values[i],
    });
  }
  return output;
}

function main() {
  const mouse = parseMouseData();
  const raw = conditionInputs();
  const current = objective(
    mouse,
    raw,
    [
      CURRENT.b,
      CURRENT.aC,
      CURRENT.aL,
      CURRENT.clockPeakHours,
      CURRENT.conservationNWeight,
      CURRENT.progesteroneCoefficient,
      CURRENT.e2Coefficient,
    ],
    true,
    false,
  );
  const fitted = fit(mouse, raw);
  const result = {
    fixed: {
      aT: CURRENT.aT,
      KT: CURRENT.KT,
      n: CURRENT.n,
      lossClockWeight: 0.78,
    },
    constantLowE2PgMl: raw.lowE2,
    before: {
      b: current.p.b,
      aC: current.p.aC,
      aL: current.p.aL,
      clockPeakHours: current.clockPeakHours,
      conservationNWeight: current.conservationNWeight,
      hormoneInputs: "Hill-transformed",
      cyclingRmse: Math.sqrt(current.cyclingMse),
      lowRmse: Math.sqrt(current.lowMse),
      metrics: current.metrics,
    },
    after: {
      b: fitted.p.b,
      aC: fitted.p.aC,
      aL: fitted.p.aL,
      clockPeakHours: fitted.clockPeakHours,
      conservationNWeight: fitted.conservationNWeight,
      progesteroneCoefficient: fitted.progesteroneCoefficient,
      e2Coefficient: fitted.e2Coefficient,
      hormoneInputs: "linear raw concentrations",
      cyclingRmse: Math.sqrt(fitted.cyclingMse),
      lowRmse: Math.sqrt(fitted.lowMse),
      metrics: fitted.metrics,
    },
    series: {
      beforeCycling: sampledSeries(current.cycling),
      afterCycling: sampledSeries(fitted.cycling),
      beforeLow: sampledSeries(current.low),
      afterLow: sampledSeries(fitted.low),
      beforePregnancy: sampledSeries(current.pregnancy),
      afterPregnancy: sampledSeries(fitted.pregnancy),
    },
  };
  const outputPath = path.join(__dirname, "age_temperature_fit_results.json");
  fs.writeFileSync(outputPath, `${JSON.stringify(result, null, 2)}\n`);
  process.stdout.write(`${JSON.stringify(result, null, 2)}\n`);
}

if (require.main === module) main();
