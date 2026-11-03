"use strict";

/*
Fit the thermal-pathway weights to the daily temperature amplitudes in the
two digitized reference datasets.  The endocrine trace is generated once at
the model's exact six-second step; temperature is integrated at five-minute
resolution during optimization so the fit has no third-party dependencies.

Run:
  node fit_temperature_amplitudes.js
*/

const base = require(
  "../_archive/meeting6_development_2026-08-24/fit_dynamic_temperature.js"
);

const HOURS_DT = 5 / 60;
const ANALYSIS_DAYS = 5;
const START = Math.round(4 * 24 / HOURS_DT);
let CYCLING_SHAPE_REFERENCE = null;

const TARGETS = {
  cycling: {
    mean: 36.9730819640,
    peak: 38.1834832471,
    trough: 35.7342702946,
    amplitude: 2.4492129526,
  },
  low: {
    mean: 36.9629038967,
    peak: 37.4629038967,
    trough: 36.4629038967,
    amplitude: 1.0000000000,
  },
  pregnancy: {
    mean: 36.8730845449,
    peak: 37.3523354536,
    trough: 36.5288116623,
    amplitude: 0.8235237913,
  },
};

const CURRENT = {
  basal: 49.44,
  conservationGain: 9.76,
  lossGain: 15.29,
  conservationClockWeight: 2.27,
  conservationNWeight: 0.60,
  progesteroneWeight: 0.0274,
  lossClockWeight: 0.78,
  e2LossWeight: 0.0000869,
};

const BOUNDS = [
  [42.0, 56.0],       // basal heat rate
  [2.0, 18.0],        // conservation gain
  [2.0, 25.0],        // hormonal-loss gain
  [0.3, 3.0],         // conservation clock weight
  [0.0, 2.0],         // N contribution to conservation
  [0.0, 0.12],        // progesterone contribution to conservation
  [0.1, 2.5],         // loss clock weight
  [0.0, 0.012],       // E2 contribution to hormonal loss
];

function constantArray(length, value) {
  const output = new Float64Array(length);
  output.fill(value);
  return output;
}

function makeRawInputs() {
  const trace = base.simulateEndocrine();
  const retained = Math.round(ANALYSIS_DAYS * 24 / HOURS_DT) + 1;
  let lowE2 = Infinity;
  for (const value of trace.E2.subarray(trace.E2.length - retained)) {
    lowE2 = Math.min(lowE2, value);
  }
  const feedbackHalf = base.constants.endocrine.K;
  const lowN = 1 / (1 + Math.pow(lowE2 / feedbackHalf, 5));
  const pregnancyE2 = 215.0;
  const pregnancyN = 1 / (
    1 + Math.pow(pregnancyE2 / feedbackHalf, 5)
  );
  return {
    cycling: { N: trace.N, E2: trace.E2, progesterone: 0.2 },
    low: {
      N: constantArray(trace.N.length, lowN),
      E2: constantArray(trace.E2.length, lowE2),
      progesterone: 0.2,
    },
    pregnancy: {
      N: constantArray(trace.N.length, pregnancyN),
      E2: constantArray(trace.E2.length, pregnancyE2),
      progesterone: 25.0,
    },
    lowE2,
  };
}

function vectorToParameters(vector) {
  return {
    basal: vector[0],
    conservationGain: vector[1],
    lossGain: vector[2],
    conservationClockWeight: vector[3],
    conservationNWeight: vector[4],
    progesteroneWeight: vector[5],
    lossClockWeight: vector[6],
    e2LossWeight: vector[7],
  };
}

function buildConditions(raw, p) {
  const override = {
    clockPeakHours: 14.0,
    conservationClockWeight: p.conservationClockWeight,
    conservationNWeight: p.conservationNWeight,
    progesteroneCoefficient: p.progesteroneWeight,
    lossClockWeight: p.lossClockWeight,
    e2Coefficient: p.e2LossWeight,
    linearHormones: true,
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

function temperatureParameters(p) {
  return {
    b: p.basal,
    aC: p.conservationGain,
    aL: p.lossGain,
    aT: 125.0,
    KT: 39.08,
    n: 5,
  };
}

function metrics(values) {
  const daily = base.dailyMetrics(values, START);
  let sum = 0;
  let count = 0;
  for (let index = START; index < values.length; index += 1) {
    sum += values[index];
    count += 1;
  }
  return { ...daily, mean: sum / count };
}

function standardizedShapeMse(values, reference) {
  let valueSum = 0;
  let referenceSum = 0;
  let count = 0;
  for (let index = START; index < values.length; index += 1) {
    valueSum += values[index];
    referenceSum += reference[index];
    count += 1;
  }
  const valueMean = valueSum / count;
  const referenceMean = referenceSum / count;
  let valueVariance = 0;
  let referenceVariance = 0;
  for (let index = START; index < values.length; index += 1) {
    valueVariance += Math.pow(values[index] - valueMean, 2);
    referenceVariance += Math.pow(reference[index] - referenceMean, 2);
  }
  const valueSd = Math.sqrt(valueVariance / count);
  const referenceSd = Math.sqrt(referenceVariance / count);
  let mse = 0;
  for (let index = START; index < values.length; index += 1) {
    const valueZ = (values[index] - valueMean) / valueSd;
    const referenceZ = (reference[index] - referenceMean) / referenceSd;
    mse += Math.pow(valueZ - referenceZ, 2);
  }
  return mse / count;
}

function objective(raw, vector, includeDetails = false) {
  const p = vectorToParameters(vector);
  const conditions = buildConditions(raw, p);
  const temperatureP = temperatureParameters(p);
  const got = {};
  let cyclingShapeMse = 0;
  let value = 0;

  for (const name of ["cycling", "low", "pregnancy"]) {
    const values = base.simulateDynamic(conditions[name], temperatureP);
    if (values === null) return { value: 1e12 };
    got[name] = metrics(values);
    if (name === "cycling" && CYCLING_SHAPE_REFERENCE !== null) {
      cyclingShapeMse = standardizedShapeMse(
        values,
        CYCLING_SHAPE_REFERENCE,
      );
      value += 100.0 * cyclingShapeMse;
      value += 1000.0 * Math.pow(
        Math.max(0.0, cyclingShapeMse - 0.04),
        2,
      );
    }
    const wanted = TARGETS[name];
    const amplitudeWeight = name === "cycling"
      ? 5.0
      : (name === "low" ? 30.0 : 20.0);
    value += amplitudeWeight * Math.pow(got[name].amplitude - wanted.amplitude, 2);
    value += 1.5 * Math.pow(got[name].mean - wanted.mean, 2);
    value += Math.pow(got[name].peak - wanted.peak, 2);
    value += Math.pow(got[name].trough - wanted.trough, 2);
  }

  value += 200.0 * Math.pow(
    got.low.amplitude / got.cycling.amplitude - 0.5,
    2,
  );
  value += 200.0 * Math.pow(
    got.pregnancy.amplitude / got.cycling.amplitude - 0.5,
    2,
  );

  // Prefer a reweighting close to the present model when two candidates have
  // nearly identical temperature behavior.
  const scales = {
    basal: 8.0,
    conservationGain: 8.0,
    lossGain: 10.0,
    conservationClockWeight: 2.0,
    conservationNWeight: 1.0,
    progesteroneWeight: 0.08,
    lossClockWeight: 1.5,
    e2LossWeight: 0.008,
  };
  for (const key of Object.keys(CURRENT)) {
    value += 0.002 * Math.pow((p[key] - CURRENT[key]) / scales[key], 2);
  }

  return includeDetails
    ? { value, parameters: p, metrics: got, cyclingShapeMse }
    : { value };
}

function rng(seed) {
  let state = seed >>> 0;
  return () => {
    state += 0x6D2B79F5;
    let value = state;
    value = Math.imul(value ^ (value >>> 15), value | 1);
    value ^= value + Math.imul(value ^ (value >>> 7), value | 61);
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
}

function fit(raw) {
  const random = rng(20260825);
  const populationSize = 64;
  const generations = 260;
  const dimension = BOUNDS.length;
  const population = [];
  const scores = [];
  const currentVector = [
    CURRENT.basal,
    CURRENT.conservationGain,
    CURRENT.lossGain,
    CURRENT.conservationClockWeight,
    CURRENT.conservationNWeight,
    CURRENT.progesteroneWeight,
    CURRENT.lossClockWeight,
    CURRENT.e2LossWeight,
  ];
  population.push(currentVector);
  scores.push(objective(raw, currentVector).value);

  for (let i = 1; i < populationSize; i += 1) {
    const vector = BOUNDS.map(([lo, hi]) => lo + random() * (hi - lo));
    population.push(vector);
    scores.push(objective(raw, vector).value);
  }

  for (let generation = 0; generation < generations; generation += 1) {
    for (let i = 0; i < populationSize; i += 1) {
      let a;
      let b;
      let c;
      do { a = Math.floor(random() * populationSize); } while (a === i);
      do { b = Math.floor(random() * populationSize); } while (b === i || b === a);
      do { c = Math.floor(random() * populationSize); } while (
        c === i || c === a || c === b
      );
      const forced = Math.floor(random() * dimension);
      const trial = [];
      for (let j = 0; j < dimension; j += 1) {
        let candidate = (random() < 0.90 || j === forced)
          ? population[a][j] + 0.70 * (population[b][j] - population[c][j])
          : population[i][j];
        const [lo, hi] = BOUNDS[j];
        if (candidate < lo || candidate > hi) {
          candidate = lo + random() * (hi - lo);
        }
        trial.push(candidate);
      }
      const trialScore = objective(raw, trial).value;
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
  return objective(raw, population[best], true);
}

function main() {
  const raw = makeRawInputs();
  const currentVector = [
    CURRENT.basal,
    CURRENT.conservationGain,
    CURRENT.lossGain,
    CURRENT.conservationClockWeight,
    CURRENT.conservationNWeight,
    CURRENT.progesteroneWeight,
    CURRENT.lossClockWeight,
    CURRENT.e2LossWeight,
  ];
  const result = {
    integration_note: "six-second endocrine; five-minute RK4 temperature fit",
    constant_low_e2_pg_ml: raw.lowE2,
    targets: TARGETS,
    before: null,
    after: null,
  };
  const currentParameters = vectorToParameters(currentVector);
  const currentConditions = buildConditions(raw, currentParameters);
  CYCLING_SHAPE_REFERENCE = base.simulateDynamic(
    currentConditions.cycling,
    temperatureParameters(currentParameters),
  );
  result.before = objective(raw, currentVector, true);
  result.after = fit(raw);
  process.stdout.write(`${JSON.stringify(result, null, 2)}\n`);
}

if (require.main === module) main();
