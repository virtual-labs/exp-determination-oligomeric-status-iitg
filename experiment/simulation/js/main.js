"use strict";

/* ================================================================
 * Determination of Oligomeric Status of a Protein — Virtual Lab
 * script.js
 * ----------------------------------------------------------------
 * Two-level structure:
 *
 *   STEP 1  SDS-PAGE ................. subunit MW (placeholder for now)
 *   STEP 2  Gel Filtration (GPC) ..... native MW, with sub-steps
 *             2.1 Column Packing
 *             2.2 Sample Preparation
 *             2.3 Running of Samples   (5 sequential runs, one curve each)
 *             2.4 Determination of MW  (Kav → calibration curve → MW)
 *             2.5 Column Regeneration
 *   STEP 3  Oligomeric Status ........ n = native MW / subunit MW
 *
 * Each animated sub-step is a "phase": it owns an SVG scene and a
 * sub-step timeline (StateMachine) and renders as a pure function of
 * its master time T, so Start/Pause/Next/Prev/Reset/scrub all work.
 * ================================================================ */

const SVG_NS = "http://www.w3.org/2000/svg";

/* ---------------- helpers ---------------- */
const clamp = (v, lo, hi) => (v < lo ? lo : v > hi ? hi : v);
const lerp = (a, b, t) => a + (b - a) * t;
const invLerp = (a, b, v) => (b === a ? 0 : (v - a) / (b - a));
const TAU = Math.PI * 2;
const log10 = (x) => Math.log(x) / Math.LN10;
const $ = (id) => document.getElementById(id);
const svgEl = (name, attrs = {}) => {
  const el = document.createElementNS(SVG_NS, name);
  for (const k in attrs) el.setAttribute(k, attrs[k]);
  return el;
};
function makePRNG(seed) {
  let s = seed >>> 0;
  return function () {
    s |= 0;
    s = (s + 0x6d2b79f5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function mixHex(a, b, t) {
  const pa = [
    parseInt(a.slice(1, 3), 16),
    parseInt(a.slice(3, 5), 16),
    parseInt(a.slice(5, 7), 16),
  ];
  const pb = [
    parseInt(b.slice(1, 3), 16),
    parseInt(b.slice(3, 5), 16),
    parseInt(b.slice(5, 7), 16),
  ];
  const c = pa.map((v, i) => Math.round(lerp(v, pb[i], t)));
  return "#" + c.map((v) => v.toString(16).padStart(2, "0")).join("");
}
const fmtClock = (t) =>
  `${Math.floor(t / 60)}:${String(Math.floor(t % 60)).padStart(2, "0")}`;

/* ================================================================
 * GPC DATA MODEL  (real values from the lab manual)
 *   Kav = (Ve - Vo) / (Vt - Vo)
 *   calibration: Kav = -m·log(MW) + c   (linear regression on standards)
 * ================================================================ */
function computeGPC() {
  const Vo = 7.79,
    Vt = 23.562;
  const kav = (Ve) => (Ve - Vo) / (Vt - Vo);

  // running order: Blue Dextran, three standards, then the unknown
  const blueDextran = {
    key: "bd",
    label: "Blue Dextran (void marker)",
    short: "BD",
    color: "#2f6fdb",
    Ve: 7.79,
    role: "void",
  };
  const standards = [
    {
      key: "s1",
      label: "β-Amylase · 200 kDa",
      short: "S1",
      color: "#0e9aa7",
      Ve: 10.0,
      mw: 200,
    },
    {
      key: "s2",
      label: "BSA · 66.4 kDa",
      short: "S2",
      color: "#8e44ad",
      Ve: 12.5,
      mw: 66.4,
    },
    {
      key: "s3",
      label: "Lysozyme · 14.3 kDa",
      short: "S3",
      color: "#e08a2e",
      Ve: 20.0,
      mw: 14.3,
    },
  ];
  const unknown = {
    key: "unk",
    label: "Unknown sample",
    short: "Unk",
    color: "#e0313b",
    Ve: 18.6,
  };

  // fill Kav / logMW
  [blueDextran, ...standards, unknown].forEach((s) => {
    s.Kav = kav(s.Ve);
  });
  standards.forEach((s) => {
    s.logMW = log10(s.mw);
  });

  // linear regression of Kav (y) on logMW (x) over the standards
  const n = standards.length;
  const mx = standards.reduce((a, s) => a + s.logMW, 0) / n;
  const my = standards.reduce((a, s) => a + s.Kav, 0) / n;
  let sxy = 0,
    sxx = 0;
  for (const s of standards) {
    sxy += (s.logMW - mx) * (s.Kav - my);
    sxx += (s.logMW - mx) ** 2;
  }
  const slope = sxy / sxx; // negative
  const c = my - slope * mx;
  const m = -slope; // report as positive m in "Kav = -m·logMW + c"

  // unknown MW from its Kav
  unknown.logMW = (c - unknown.Kav) / m;
  unknown.mw = Math.pow(10, unknown.logMW);

  const runs = [blueDextran, ...standards, unknown];
  // travel speed proxy: bigger Kav (smaller molecule) arrives later
  runs.forEach((r) => {
    r.arrivalP = clamp(0.45 + r.Kav * 0.4, 0.45, 0.82);
    r.wiggle = 1.5 + r.Kav * 6;
  });

  const sdsPageMW = 18; // subunit MW from Step 1 (doc value / placeholder)
  const oligomeric = unknown.mw / sdsPageMW;

  return {
    Vo,
    Vt,
    blueDextran,
    standards,
    unknown,
    runs,
    m,
    c,
    sdsPageMW,
    oligomeric,
  };
}
const GPC = computeGPC();

/* ================================================================
 * SDS-PAGE DATA MODEL
 *   Rf = distance(protein) / distance(dye front)
 *   calibration: log(MW) = slope·Rf + intercept   (slope < 0)
 * ================================================================ */
function computeSDS() {
  const markers = [
    { name: "Phosphorylase b", mw: 97, Rf: 0.15, color: "#0e9aa7" },
    { name: "BSA", mw: 66, Rf: 0.28, color: "#2f6fdb" },
    { name: "Ovalbumin", mw: 45, Rf: 0.42, color: "#8e44ad" },
    { name: "Carbonic anhyd.", mw: 30, Rf: 0.58, color: "#2faf6b" },
    { name: "Trypsin inhib.", mw: 20.1, Rf: 0.72, color: "#e0a32e" },
    { name: "Lysozyme", mw: 14.4, Rf: 0.85, color: "#e08a2e" },
  ];
  markers.forEach((m) => {
    m.logMW = log10(m.mw);
  });
  const n = markers.length;
  const mx = markers.reduce((a, m) => a + m.Rf, 0) / n;
  const my = markers.reduce((a, m) => a + m.logMW, 0) / n;
  let sxy = 0,
    sxx = 0;
  for (const m of markers) {
    sxy += (m.Rf - mx) * (m.logMW - my);
    sxx += (m.Rf - mx) ** 2;
  }
  const slope = sxy / sxx; // negative
  const intercept = my - slope * mx;
  const unknown = { Rf: 0.764, color: "#e0313b" };
  unknown.logMW = slope * unknown.Rf + intercept;
  unknown.mw = Math.pow(10, unknown.logMW);
  return { markers, slope, intercept, m: -slope, c: intercept, unknown };
}
const SDS = computeSDS();

/* ================================================================
 * CONFIG
 * ================================================================ */
const CONFIG = {
  geometry: {
    column: { travelTopY: 176, outletY: 397, innerL: 408, innerR: 446 },
    collector: {
      tubeCount: 6,
      tubeTopY: 502,
      tubeH: 32,
      tubeW: 30,
      tubeBottomY: 534,
      nozzleY: 486,
      centers: [524, 572, 621, 670, 719, 768],
    },
  },
  chart: {
    x0: 38,
    x1: 292,
    yBase: 150,
    yTop: 16,
    veMin: 6,
    veMax: 24,
    gridRows: 4,
  },
  transport: { defaultSpeed: 1, minSpeed: 0.25, maxSpeed: 3 },
  equipment: { valveInjectAngle: 62, plungerTravel: 16 },

  phases: {
    sdsprep: [
      {
        key: "mix",
        name: "Add Dye to Protein",
        duration: 6,
        instruction:
          "The protein and SDS loading dye (with bromophenol blue and a reducing agent) are poured into the tube and mixed. SDS unfolds the protein and gives it a uniform negative charge.",
        note: "Protein + SDS loading dye added to the tube.",
        hint: {
          text: "Prepare the protein sample for SDS-PAGE",
          action: "play",
        },
      },
      {
        key: "heat",
        name: "Heat (95 °C)",
        duration: 6,
        instruction:
          "The tube is held in the heating block at 95 °C. Heating with SDS and the reducing agent fully denatures the protein into its subunits.",
        note: "Heated to 95 °C — protein fully denatured.",
        hint: null,
      },
      {
        key: "spin",
        name: "Centrifuge",
        duration: 5,
        instruction:
          "The tube is briefly centrifuged to pellet any insoluble debris, leaving a clear blue supernatant.",
        note: "Centrifuged — insoluble debris pelleted.",
        hint: null,
      },
      {
        key: "ready",
        name: "Sample Ready",
        duration: 3,
        instruction:
          "The denatured, negatively-charged blue sample is ready to load onto the gel.",
        note: "Blue SDS sample ready to load.",
        hint: { text: "Go to next step: Run on SDS-PAGE", action: "nextsub" },
      },
    ],
    sdsrun: [
      {
        key: "load",
        name: "Load the Gel",
        duration: 5,
        instruction:
          "The molecular-weight marker (a protein ladder) and the unknown sample are loaded into the wells at the top of the gel.",
        note: "Marker and unknown loaded into the wells.",
        hint: null,
      },
      {
        key: "run",
        name: "Apply Voltage",
        duration: 8,
        instruction:
          "A voltage is applied. SDS-coated proteins are pulled toward the anode (+) and separate by size — smaller proteins move faster and farther.",
        note: "Proteins migrate and separate by size.",
        hint: null,
      },
      {
        key: "stain",
        name: "Stain the Bands",
        duration: 5,
        instruction:
          "The gel is stained with Coomassie Brilliant Blue, revealing the separated protein bands. Each band's migration distance is measured for the next step.",
        note: "Bands revealed by Coomassie; distances measured.",
        hint: {
          text: "Go to next step: Determination of MW",
          action: "nextsub",
        },
      },
    ],
    sdsdetermine: [
      {
        key: "rf",
        name: "Plot Rf vs log(MW)",
        duration: 5,
        instruction:
          "For each band the relative mobility Rf = d(protein)/d(dye front) is measured and plotted against log(MW).",
        note: "Standards plotted — Rf vs log(MW).",
        hint: { text: "Draw the calibration curve", action: "play" },
      },
      {
        key: "line",
        name: "Fit the Line",
        duration: 5,
        instruction:
          "A straight line is fitted through the standards: Rf = −m·log(MW) + c.",
        note: "Calibration line fitted (y = −mx + c).",
        hint: null,
      },
      {
        key: "mw",
        name: "Read the Unknown",
        duration: 6,
        instruction:
          "The unknown's Rf is located on the line (★) to read its log(MW), giving the subunit molecular weight.",
        note: "Subunit (denatured) MW read from the calibration curve.",
        hint: {
          text: "Go to next step: Determination of native Mol. Wt.",
          action: "nextstep",
        },
      },
    ],
    packing: [
      {
        key: "swell",
        name: "Matrix Swelling",
        duration: 6,
        instruction:
          "The dry gel matrix is allowed to swell in the mobile phase, forming a uniform, pourable slurry.",
        note: "Matrix swollen in mobile phase.",
        hint: { text: "▶ Pack the column", action: "play" },
      },
      {
        key: "pour",
        name: "Pouring the Slurry",
        duration: 6,
        instruction:
          "The slurry is poured into the column down a glass rod in one continuous motion, so that no air is trapped within the bed.",
        note: "Slurry poured into the column.",
        hint: null,
      },
      {
        key: "settle",
        name: "Bed Settling",
        duration: 7,
        instruction:
          "The beads settle under buffer flow into an even packed bed. Trapped air bubbles rise and escape, leaving a uniform column ready for use.",
        note: "Beads settled — no air channels; column packed.",
        hint: {
          text: "Go to next step: Sample Preparation",
          action: "nextsub",
        },
      },
    ],
    prep: [
      {
        key: "dissolve",
        name: "Dissolve in Mobile Phase",
        duration: 6,
        instruction:
          "The protein sample is dissolved in the mobile phase — the same buffer used to elute — so sample and eluent are matched.",
        note: "Sample dissolved in mobile phase.",
        hint: { text: "▶ Prepare the sample", action: "play" },
      },
      {
        key: "clarify",
        name: "Remove Particulates",
        duration: 7,
        instruction:
          "The solution is pushed through a syringe filter to remove suspended particles that would clog the column.",
        note: "Suspended particles removed by filtration.",
        hint: null,
      },
      {
        key: "load",
        name: "Draw into Syringe",
        duration: 6,
        instruction:
          "The clarified sample is drawn into a syringe — the recommended way to apply it — ready to inject onto the column.",
        note: "Clarified sample loaded into the injection syringe.",
        hint: {
          text: "Go to next step: Running of Samples",
          action: "nextsub",
        },
      },
    ],
    regen: [
      {
        key: "wash",
        name: "Salt Wash",
        duration: 7,
        instruction:
          "The column is washed with a salt-containing mobile phase to strip off protein adsorbed non-specifically to the matrix.",
        note: "Salt wash removes non-specifically adsorbed protein.",
        hint: { text: "▶ Regenerate the column", action: "play" },
      },
      {
        key: "equilibrate",
        name: "Re-equilibration",
        duration: 6,
        instruction:
          "The column is re-equilibrated with mobile phase, restoring the starting conditions for the next separation.",
        note: "Column re-equilibrated with mobile phase.",
        hint: null,
      },
      {
        key: "store",
        name: "Storage",
        duration: 6,
        instruction:
          "For storage the column is kept at 4 °C in 20% alcohol with 0.05% sodium azide, preventing microbial growth.",
        note: "Stored at 4 °C in 20% ethanol + 0.05% sodium azide.",
        hint: { text: "→ Next: Oligomeric Status", action: "nextstep" },
      },
    ],
    determination: [
      {
        key: "kav",
        name: "Calculate Kₐᵥ",
        duration: 5,
        instruction:
          "For every sample the distribution coefficient Kₐᵥ is found from its elution volume: Kₐᵥ = (Vₑ − Vₒ)/(Vₜ − Vₒ).",
        note: "Kₐᵥ computed for each protein from its Vₑ.",
        hint: { text: "▶ Draw the calibration curve", action: "play" },
      },
      {
        key: "plot",
        name: "Plot the Standards",
        duration: 5,
        instruction: "Each standard is plotted as Kₐᵥ against log(MW).",
        note: "Standard proteins plotted (Kₐᵥ vs log MW).",
        hint: null,
      },
      {
        key: "line",
        name: "Fit the Line",
        duration: 5,
        instruction:
          "A straight line is fitted through the standards: Kₐᵥ = −m·log(MW) + c.",
        note: "Calibration line fitted through the standards.",
        hint: null,
      },
      {
        key: "read",
        name: "Read the Unknown",
        duration: 6,
        instruction:
          "The unknown's Kₐᵥ is located on the line to read its log(MW), giving the native molecular weight.",
        note: "Native MW of the unknown read from the calibration curve.",
        hint: {
          text: "Go to next step: Oligomeric Status",
          action: "nextstep",
        },
      },
    ],
  },
};

/* ================================================================
 * STATE MACHINE
 * ================================================================ */
class StateMachine {
  constructor(steps) {
    let t = 0;
    this.stages = steps.map((s) => {
      const st = { ...s, startT: t, endT: t + s.duration };
      t += s.duration;
      return st;
    });
    this.totalT = t;
  }
  indexAt(T) {
    if (T >= this.totalT) return this.stages.length - 1;
    for (let i = 0; i < this.stages.length; i++)
      if (T < this.stages[i].endT) return i;
    return this.stages.length - 1;
  }
  stageAt(T) {
    return this.stages[this.indexAt(T)];
  }
  localAt(T) {
    const i = this.indexAt(T);
    const s = this.stages[i];
    return { i, s, p: clamp(invLerp(s.startT, s.endT, T), 0, 1) };
  }
}

/* ================================================================
 * PHASE 2.1 — COLUMN PACKING
 * ================================================================ */
class ColumnPackingPhase {
  constructor(sceneEl) {
    this.key = "packing";
    this.code = "2.1";
    this.label = "Column Packing";
    this.scene = sceneEl;
    this.usesChart = false;
    this.sm = new StateMachine(CONFIG.phases.packing);
    this.T = 0;
    this.geo = { bedBottom: 452, bedFullTop: 160 };
    this.el = {
      bed: $("p1Bed"),
      liquid: $("p1Liquid"),
      band: $("p1Band"),
      slurry: $("p1Slurry"),
      slurryBeads: $("p1SlurryBeads"),
      slurrySurf: $("p1SlurrySurf"),
      pour: $("p1PourGroup"),
      bubbles: $("p1Bubbles"),
      waste: $("p1WasteLiquid"),
      result: $("p1Result"),
    };
    const rng = makePRNG(0xb0bb1e);
    this.bubbles = [];
    for (let i = 0; i < 7; i++) {
      const c = svgEl("circle", {
        r: (1.5 + rng() * 2.4).toFixed(1),
        fill: "#eafaff",
        opacity: "0",
      });
      this.el.bubbles.appendChild(c);
      this.bubbles.push({
        el: c,
        x: 552 + rng() * 68,
        off: rng(),
        speed: 0.7 + rng() * 0.6,
      });
    }
  }
  render(T) {
    const { s, p } = this.sm.localAt(T);
    const { bedBottom, bedFullTop } = this.geo;
    const bedFullH = bedBottom - bedFullTop;
    let bedFrac = 0,
      slurryFrac = 1,
      liquidShow = 0,
      pouring = false,
      bubbleP = -1,
      bandP = -1,
      wasteFrac = 0,
      resultO = 0,
      testReveal = 0;
    switch (s.key) {
      case "swell":
        bedFrac = 0;
        slurryFrac = 1;
        break;
      case "pour":
        pouring = true;
        slurryFrac = 1 - p;
        bedFrac = 0.6 * p;
        liquidShow = 0.35 * p;
        break;
      case "settle":
        slurryFrac = 0;
        bedFrac = lerp(0.6, 1, p);
        liquidShow = lerp(0.35, 1, p);
        bubbleP = p;
        wasteFrac = 0.14 * p;
        break;
      case "test":
        slurryFrac = 0;
        bedFrac = 1;
        liquidShow = 1;
        bandP = p;
        testReveal = p;
        wasteFrac = 0.14 + 0.55 * clamp((p - 0.75) / 0.25, 0, 1);
        resultO = clamp((p - 0.85) / 0.12, 0, 1);
        break;
    }
    const bedH = bedFullH * bedFrac,
      bedTop = bedBottom - bedH;
    this.el.bed.setAttribute("y", bedTop.toFixed(1));
    this.el.bed.setAttribute("height", bedH.toFixed(1));
    const liqTop = 130;
    this.el.liquid.setAttribute("y", liqTop);
    this.el.liquid.setAttribute(
      "height",
      Math.max(0, bedTop - liqTop).toFixed(1),
    );
    this.el.liquid.setAttribute("opacity", (0.45 * liquidShow).toFixed(2));
    const slurryH = 140 * slurryFrac,
      slurryY = 300 + (140 - slurryH);
    for (const r of [this.el.slurry, this.el.slurryBeads]) {
      r.setAttribute("y", slurryY.toFixed(1));
      r.setAttribute("height", slurryH.toFixed(1));
    }
    const bob = s.key === "swell" ? Math.sin(T * 3) * 1.5 : 0;
    this.el.slurrySurf.setAttribute("cy", (slurryY + bob).toFixed(1));
    this.el.slurrySurf.setAttribute("opacity", (0.85 * slurryFrac).toFixed(2));
    this.el.pour.classList.toggle("is-flowing", pouring);
    for (const b of this.bubbles) {
      if (bubbleP < 0) {
        b.el.setAttribute("opacity", "0");
        continue;
      }
      const prog = (bubbleP * b.speed + b.off) % 1;
      b.el.setAttribute("cx", b.x.toFixed(1));
      b.el.setAttribute(
        "cy",
        lerp(bedBottom - 10, bedFullTop + 6, prog).toFixed(1),
      );
      b.el.setAttribute("opacity", ((1 - bubbleP) * 0.6).toFixed(2));
    }
    if (bandP < 0) {
      this.el.band.setAttribute("opacity", "0");
      this.el.band.setAttribute("height", "0");
    } else {
      const bandH = 16,
        by = lerp(bedFullTop, bedBottom - bandH, bandP);
      this.el.band.setAttribute("y", by.toFixed(1));
      this.el.band.setAttribute("height", bandH);
      const fade = bandP < 0.94 ? 1 : 1 - (bandP - 0.94) / 0.06;
      this.el.band.setAttribute(
        "opacity",
        (0.9 * clamp(fade, 0, 1)).toFixed(2),
      );
    }
    const wasteH = 40 * wasteFrac;
    this.el.waste.setAttribute("height", wasteH.toFixed(1));
    this.el.waste.setAttribute("y", (520 - wasteH).toFixed(1));
    this.el.result.setAttribute("opacity", resultO.toFixed(2));
  }
  events() {}
}

/* ================================================================
 * PHASE 2.2 — SAMPLE PREPARATION
 * ================================================================ */
class SamplePrepPhase {
  constructor(sceneEl) {
    this.key = "prep";
    this.code = "2.2";
    this.label = "Sample Preparation";
    this.scene = sceneEl;
    this.usesChart = false;
    this.sm = new StateMachine(CONFIG.phases.prep);
    this.T = 0;
    this.el = {
      vialLiquid: $("p2VialLiquid"),
      cloud: $("p2Cloud"),
      drop: $("p2DropperDrop"),
      barrel: $("p2BarrelSample"),
      plunger: $("p2Plunger"),
      cake: $("p2FilterCake"),
      fdrop: $("p2Drop"),
      filtrate: $("p2Filtrate"),
      loadSyr: $("p2LoadSyringe"),
      loadSample: $("p2LoadSample"),
      loadPlunger: $("p2LoadPlunger"),
    };
    const rng = makePRNG(0xc10d);
    this.specks = [];
    for (let i = 0; i < 12; i++) {
      const c = svgEl("circle", {
        r: (1.3 + rng() * 1.6).toFixed(1),
        fill: "#8a5a3a",
        opacity: "0",
      });
      this.el.cloud.appendChild(c);
      this.specks.push({
        el: c,
        x: 128 + rng() * 54,
        y: 278 + rng() * 140,
        ph: rng() * TAU,
      });
    }
  }
  render(T) {
    const { s, p } = this.sm.localAt(T);
    let vialFrac = 1,
      cloudO = 0,
      barrelFrac = 0,
      plungerY = 0,
      cakeO = 0,
      dropOn = false,
      filtrateFrac = 0,
      loadO = 0,
      loadFrac = 0;
    switch (s.key) {
      case "dissolve":
        vialFrac = p;
        cloudO = p;
        break;
      case "clarify":
        vialFrac = lerp(1, 0.15, p);
        cloudO = 1 - p;
        barrelFrac = 1 - p;
        plungerY = 130 * p;
        cakeO = clamp(p * 1.3, 0, 0.9);
        filtrateFrac = p;
        dropOn = true;
        break;
      case "load":
        vialFrac = 0.15;
        cloudO = 0;
        cakeO = 0.9;
        barrelFrac = 0;
        loadO = 1;
        loadFrac = p;
        filtrateFrac = 1 - p;
        break;
    }
    // vial liquid (bottom at 430)
    const vialH = 172 * vialFrac;
    this.el.vialLiquid.setAttribute("height", vialH.toFixed(1));
    this.el.vialLiquid.setAttribute("y", (430 - vialH).toFixed(1));
    for (const sp of this.specks) {
      sp.el.setAttribute(
        "cx",
        (sp.x + Math.sin(sp.ph + T * 2.2) * 2.4).toFixed(1),
      );
      sp.el.setAttribute(
        "cy",
        (sp.y + Math.cos(sp.ph + T * 1.8) * 2.4).toFixed(1),
      );
      sp.el.setAttribute("opacity", (cloudO * 0.85).toFixed(2));
    }
    this.el.drop.setAttribute(
      "r",
      (dropOn
        ? 0
        : s.key === "dissolve"
          ? Math.max(0, 2.4 + Math.sin(T * 8) * 1.6)
          : 0
      ).toFixed(1),
    );
    // syringe barrel (drains from bottom 315 as plunger pushes)
    const barH = 140 * barrelFrac;
    this.el.barrel.setAttribute("height", barH.toFixed(1));
    this.el.barrel.setAttribute("y", (315 - barH).toFixed(1));
    this.el.plunger.setAttribute(
      "transform",
      `translate(0,${plungerY.toFixed(1)})`,
    );
    this.el.cake.setAttribute("opacity", cakeO.toFixed(2));
    this.el.fdrop.setAttribute(
      "r",
      (dropOn ? Math.max(0, 2.2 + Math.sin(T * 10) * 1.4) : 0).toFixed(1),
    );
    // clarified filtrate tube (bottom at 450)
    const filH = 90 * filtrateFrac;
    this.el.filtrate.setAttribute("height", filH.toFixed(1));
    this.el.filtrate.setAttribute("y", (450 - filH).toFixed(1));
    // injection syringe: plunger is pulled up as the sample is drawn in from
    // the bottom (400). Its tip tracks the rising liquid surface.
    this.el.loadSyr.setAttribute("opacity", loadO.toFixed(2));
    const loadH = 140 * loadFrac;
    this.el.loadSample.setAttribute("height", loadH.toFixed(1));
    this.el.loadSample.setAttribute("y", (400 - loadH).toFixed(1));
    this.el.loadPlunger.setAttribute(
      "transform",
      `translate(0,${(146 - 140 * loadFrac).toFixed(1)})`,
    );
  }
  events() {}
}

/* ================================================================
 * PHASE 2.3 — RUNNING OF SAMPLES  (5 sequential accumulating runs)
 * ================================================================ */
class RunningPhase {
  constructor(sceneEl) {
    this.key = "running";
    this.code = "2.3";
    this.label = "Running of Samples";
    this.scene = sceneEl;
    this.usesChart = true;
    this.stepwise = true; // each sample is run by a separate click
    this.chartCaption =
      "Each sample run in turn — peaks accumulate by elution volume";
    this.runs = GPC.runs;
    this.chartLegend = this.runs
      .map(
        (r) =>
          `<span><i class="dot" style="background:${r.color}"></i>${r.short}</span>`,
      )
      .join("");
    this.geo = CONFIG.geometry;

    // descriptive run-button text; standards are numbered by their run order
    const runHint = (r, i) =>
      r.role === "void"
        ? "Run blue dextran to calculate the void volume"
        : r.key === "unk"
          ? "Run unknown protein sample"
          : `Run standard ${i} ${r.label.split(" ·")[0]}`;

    // one stage per run (each fired by its own click) + a final summary
    const stages = this.runs.map((r, i) => ({
      key: `run${i}`,
      name: `Run ${i + 1}: ${r.short}`,
      duration: 7,
      instruction: `Running <b>${r.label}</b>. ${
        r.role === "void"
          ? "Being the largest, it is fully excluded and elutes first, marking the void volume Vₒ."
          : r.key === "unk"
            ? "Its peak appears between the standards — this fixes the unknown's elution volume."
            : "It elutes at a volume set by its size — larger proteins come off earlier."
      }`,
      note: `${r.label} eluted at Vₑ = ${r.Ve.toFixed(2)} mL (Kₐᵥ = ${r.Kav.toFixed(3)}).`,
      hint: { text: runHint(r, i), action: "play" },
    }));
    stages.push({
      key: "done",
      name: "All Runs Complete",
      duration: 3,
      instruction:
        "All five chromatograms are recorded. The unknown's peak sits between the standards. Proceed to determine its molecular weight.",
      note: "Blue Dextran, three standards and the unknown all run.",
      hint: { text: "Go to next step: Determination of MW", action: "nextsub" },
    });
    this.sm = new StateMachine(stages);
    this.T = 0;

    this._refs();
    this._buildTubes();
    this._buildMolecules();
    this._buildPeaks();
  }
  _refs() {
    this.el = {
      lab: $("lab"),
      tubings: Array.from($("lab").querySelectorAll(".tubing")),
      valveRotor: $("valveRotor"),
      plunger: $("syringePlunger"),
      syrSample: $("syringeSample"),
      band: $("sampleBand"),
      glow: $("detectorGlow"),
      lamp: $("detectorLamp"),
      nozzle: $("nozzle"),
      drop: $("nozzleDrop"),
      layer: $("moleculeLayer"),
      rack: $("tubeRack"),
      chip: $("runChip"),
      chipDot: $("runChipDot"),
      chipText: $("runChipText"),
      chipBg: $("runChipBg"),
      runsG: $("chartRuns"),
      cursor: $("chartCursor"),
    };
  }
  _buildTubes() {
    const g = this.geo.collector;
    this.tubes = [];
    for (let i = 0; i < g.tubeCount; i++) {
      const x = g.centers[i] - g.tubeW / 2;
      const grp = svgEl("g", { class: "tube" });
      const glass = svgEl("rect", {
        class: "tube-glass",
        x,
        y: g.tubeTopY,
        width: g.tubeW,
        height: g.tubeH,
        rx: 7,
        fill: "url(#gradGlass)",
        stroke: "#9fb4c2",
        "stroke-width": "1.5",
      });
      const liquid = svgEl("rect", {
        class: "tube-liquid",
        x: x + 2,
        width: g.tubeW - 4,
        y: g.tubeBottomY,
        height: 0,
        rx: 4,
        fill: this.runs[i] ? this.runs[i].color : "#c9d6de",
      });
      const label = svgEl("text", {
        class: "svg-tag",
        x: g.centers[i],
        y: g.tubeTopY - 4,
        "text-anchor": "middle",
      });
      label.textContent = this.runs[i] ? this.runs[i].short : `F${i + 1}`;
      grp.append(glass, liquid, label);
      this.el.rack.appendChild(grp);
      this.tubes.push({ liquid });
    }
  }
  _buildMolecules() {
    const rng = makePRNG(0x123abc);
    const c = this.geo.column;
    this.mols = [];
    for (let i = 0; i < 12; i++) {
      const el = svgEl("circle", {
        r: 5,
        opacity: "0",
        stroke: "rgba(255,255,255,.5)",
        "stroke-width": "0.8",
      });
      this.el.layer.appendChild(el);
      this.mols.push({
        el,
        lane: lerp(c.innerL + 8, c.innerR - 8, rng()),
        ph: rng() * TAU,
      });
    }
  }
  _buildPeaks() {
    this.peaks = this.runs.map((r) => {
      const path = svgEl("path", {
        fill: r.color,
        "fill-opacity": "0.32",
        stroke: r.color,
        "stroke-width": "2",
        d: "",
      });
      const label = svgEl("text", {
        class: "svg-tag",
        fill: r.color,
        "text-anchor": "middle",
        opacity: "0",
      });
      label.textContent = r.short;
      this.el.runsG.append(path, label);
      return { path, label, run: r };
    });
  }
  _xOf(ve) {
    const { x0, x1, veMin, veMax } = CONFIG.chart;
    return lerp(x0, x1, clamp((ve - veMin) / (veMax - veMin), 0, 1));
  }

  render(T, playing) {
    const { i: si, s: stage, p } = this.sm.localAt(T);
    const runs = this.runs,
      isDone = stage.key === "done";
    const cur = isDone ? null : runs[si];
    const c = this.geo.column,
      col = this.geo.collector,
      eq = CONFIG.equipment;

    // ---- accumulating chromatogram ----
    const { yBase, yTop } = CONFIG.chart,
      amp = (yBase - yTop) * 0.82,
      sig = 0.5;
    const elutedOf = (j) => {
      if (isDone || j < si) return 1;
      if (j === si) return clamp((p - runs[j].arrivalP + 0.05) / 0.22, 0, 1);
      return 0;
    };
    this.peaks.forEach((pk, j) => {
      const e = elutedOf(j);
      if (e <= 0) {
        pk.path.setAttribute("d", "");
        pk.label.setAttribute("opacity", "0");
        return;
      }
      let d = "",
        first = null,
        last = 0;
      for (let k = 0; k <= 70; k++) {
        const ve = lerp(CONFIG.chart.veMin, CONFIG.chart.veMax, k / 70);
        const y =
          yBase - amp * e * Math.exp(-0.5 * ((ve - pk.run.Ve) / sig) ** 2);
        const x = this._xOf(ve);
        d += `${k === 0 ? "M" : "L"}${x.toFixed(1)},${y.toFixed(1)} `;
        if (first === null) first = x;
        last = x;
      }
      pk.path.setAttribute(
        "d",
        d + `L${last.toFixed(1)},${yBase} L${first.toFixed(1)},${yBase} Z`,
      );
      const lx = this._xOf(pk.run.Ve),
        ly = yBase - amp * e - 5;
      pk.label.setAttribute("x", lx.toFixed(1));
      pk.label.setAttribute("y", Math.max(yTop + 8, ly).toFixed(1));
      pk.label.setAttribute("opacity", e > 0.25 ? "1" : "0");
    });

    // ---- apparatus for the current run ----
    const injectEnd = 0.12;
    let flowing = false,
      plungerP = 0,
      bandO = 0,
      glow = 0,
      curReveal = 0;
    const tubeIdx = isDone ? -1 : si;

    if (!isDone) {
      const arr = cur.arrivalP;
      curReveal = elutedOf(si);
      glow = curReveal;
      // syringe + injection
      plungerP = clamp(p / injectEnd, 0, 1);
      bandO =
        p < injectEnd ? plungerP * 0.5 : Math.max(0, 0.5 - (p - injectEnd) * 3);
      flowing = p >= injectEnd && p < arr + 0.2;
      // molecules
      for (const m of this.mols) {
        let x = m.lane,
          y = c.travelTopY,
          o = 0;
        if (p < injectEnd) {
          y = c.travelTopY;
          o = plungerP;
          x = m.lane + Math.sin(m.ph + T) * 1.5;
        } else if (p <= arr) {
          const u = invLerp(injectEnd, arr, p);
          y = lerp(c.travelTopY, c.outletY, u);
          x = m.lane + Math.sin(m.ph + T * 1.4) * cur.wiggle * (0.5 + u);
          o = 1;
        } else {
          const u2 = clamp(invLerp(arr, arr + 0.18, p), 0, 1);
          x = lerp(m.lane, col.centers[tubeIdx], u2 * u2);
          y = lerp(c.outletY, col.tubeTopY, u2);
          o = 1 - 0.4 * u2;
        }
        m.el.setAttribute("cx", x.toFixed(1));
        m.el.setAttribute("cy", y.toFixed(1));
        m.el.setAttribute("opacity", o.toFixed(2));
        m.el.setAttribute("fill", cur.color);
      }
    } else {
      for (const m of this.mols) m.el.setAttribute("opacity", "0");
    }

    // band, syringe
    this.el.band.setAttribute("opacity", bandO.toFixed(2));
    this.el.band.setAttribute("fill", cur ? cur.color : "#b5532f");
    this.el.band.setAttribute(
      "height",
      (13 * clamp(bandO / 0.5, 0, 1)).toFixed(1),
    );
    this.el.band.setAttribute("y", "170");
    this.el.plunger.setAttribute(
      "transform",
      `translate(0,${(plungerP * eq.plungerTravel).toFixed(1)})`,
    );
    this.el.syrSample.setAttribute("height", (14 * (1 - plungerP)).toFixed(1));
    this.el.syrSample.setAttribute(
      "opacity",
      (0.9 * (1 - plungerP)).toFixed(2),
    );
    this.el.syrSample.setAttribute("fill", cur ? cur.color : "#c0504d");

    // flow decorations
    for (const g of this.el.tubings) g.classList.toggle("is-flowing", flowing);
    this.el.lab.classList.toggle("is-spinning", flowing && playing);
    this.el.valveRotor.setAttribute(
      "transform",
      `rotate(${flowing ? eq.valveInjectAngle : 0} 355 104)`,
    );

    // detector
    this.el.glow.setAttribute("opacity", (glow * 0.9).toFixed(2));
    this.el.lamp.setAttribute("fill", glow > 0.15 ? "#ffd34d" : "#3a5666");

    // fraction tubes fill
    this.tubes.forEach((t, j) => {
      const e = elutedOf(j);
      const h = e * (col.tubeH - 6);
      t.liquid.setAttribute("height", h.toFixed(1));
      t.liquid.setAttribute("y", (col.tubeBottomY - h).toFixed(1));
    });

    // nozzle + drop
    const nx = col.centers[tubeIdx < 0 ? 0 : tubeIdx];
    this.el.nozzle.setAttribute("transform", `translate(${nx},${col.nozzleY})`);
    const dropping = !isDone && curReveal > 0 && curReveal < 1;
    this.el.drop.setAttribute(
      "r",
      dropping ? (2.4 + Math.sin(T * 12) * 1.2).toFixed(2) : "0",
    );
    this.el.drop.setAttribute("fill", cur ? cur.color : "#b5532f");

    // active-sample chip
    this.el.chipDot.setAttribute("fill", cur ? cur.color : "#2faf6b");
    this.el.chipText.textContent = cur ? cur.label : "All runs complete";
    this.el.chipBg.setAttribute(
      "width",
      cur ? Math.max(180, cur.label.length * 9 + 40) : 190,
    );

    // chromatogram cursor at the current forming peak
    this.el.cursor.setAttribute(
      "x1",
      this._xOf(cur ? cur.Ve : CONFIG.chart.veMin).toFixed(1),
    );
    this.el.cursor.setAttribute(
      "x2",
      this._xOf(cur ? cur.Ve : CONFIG.chart.veMin).toFixed(1),
    );
  }
  events(T, log) {
    const { i: si, s: stage, p } = this.sm.localAt(T);
    if (stage.key === "done") return;
    const r = this.runs[si];
    if (p > r.arrivalP + 0.1)
      log(
        `run-${r.key}`,
        `<strong>${r.short}</strong> (${r.label}) eluted at Vₑ = ${r.Ve.toFixed(2)} mL.`,
        "ok",
      );
  }
}

/* ================================================================
 * PHASE 2.4 — DETERMINATION OF MOLECULAR WEIGHT
 * ================================================================ */
class DeterminationPhase {
  constructor(sceneEl) {
    this.key = "determination";
    this.code = "2.4";
    this.label = "Determination of MW";
    this.scene = sceneEl;
    this.usesChart = false;
    this.sm = new StateMachine(CONFIG.phases.determination);
    this.T = 0;
    this.plot = {
      x0: 80,
      x1: 520,
      y0: 380,
      y1: 60,
      L0: 1.0,
      L1: 2.5,
      K0: 0,
      K1: 1.0,
    };
    this._buildGrid();
    this._buildPoints();
    this._buildLine();
    this._buildUnknown();
    this._buildFormulas();
  }
  _xOf(logMW) {
    const p = this.plot;
    return lerp(p.x0, p.x1, invLerp(p.L0, p.L1, logMW));
  }
  _yOf(kav) {
    const p = this.plot;
    return lerp(p.y0, p.y1, invLerp(p.K0, p.K1, kav));
  }
  _buildGrid() {
    const g = $("calibGrid");
    const p = this.plot;
    for (let k = 0; k <= 4; k++) {
      const kav = k / 4,
        y = this._yOf(kav);
      g.appendChild(
        svgEl("line", {
          class: "calib-grid-line",
          x1: p.x0,
          y1: y,
          x2: p.x1,
          y2: y,
        }),
      );
      const t = svgEl("text", {
        class: "calib-tick",
        x: p.x0 - 8,
        y: y + 4,
        "text-anchor": "end",
      });
      t.textContent = kav.toFixed(2);
      g.appendChild(t);
    }
    for (const lm of [1.0, 1.5, 2.0, 2.5]) {
      const x = this._xOf(lm);
      const t = svgEl("text", {
        class: "calib-tick",
        x,
        y: p.y0 + 18,
        "text-anchor": "middle",
      });
      t.textContent = lm.toFixed(1);
      g.appendChild(t);
    }
  }
  _buildPoints() {
    const g = $("calibPoints");
    this.ptEls = [];
    for (const s of GPC.standards) {
      const cxp = this._xOf(s.logMW),
        cyp = this._yOf(s.Kav);
      const dot = svgEl("circle", {
        class: "calib-pt",
        cx: cxp,
        cy: cyp,
        r: 6,
        fill: s.color,
        opacity: "0",
      });
      const lab = svgEl("text", {
        class: "calib-pt-label",
        x: cxp + 10,
        y: cyp - 6,
        opacity: "0",
      });
      lab.textContent = `${s.mw} kDa`;
      g.append(dot, lab);
      this.ptEls.push({ dot, lab });
    }
  }
  _buildLine() {
    // line spans from logMW at Kav=K1 side to where Kav hits 0, clamped to domain
    const p = this.plot;
    const lmAtK = (kav) => (GPC.c - kav) / GPC.m;
    this.lineA = { lm: clamp(lmAtK(p.K1), p.L0, p.L1) }; // upper-left
    this.lineB = { lm: clamp(lmAtK(0), p.L0, p.L1) }; // lower-right
    this.lineA.x = this._xOf(this.lineA.lm);
    this.lineA.y = this._yOf(GPC.c - GPC.m * this.lineA.lm);
    this.lineB.x = this._xOf(this.lineB.lm);
    this.lineB.y = this._yOf(GPC.c - GPC.m * this.lineB.lm);
    this.lineEl = svgEl("line", {
      class: "calib-line",
      x1: this.lineA.x,
      y1: this.lineA.y,
      x2: this.lineA.x,
      y2: this.lineA.y,
      opacity: "0",
    });
    $("calibLine").appendChild(this.lineEl);
  }
  _buildUnknown() {
    const g = $("calibUnknown");
    const u = GPC.unknown;
    const ux = this._xOf(u.logMW),
      uy = this._yOf(u.Kav);
    this.uGuideH = svgEl("line", {
      class: "calib-guide",
      x1: this.plot.x0,
      y1: uy,
      x2: ux,
      y2: uy,
      opacity: "0",
    });
    this.uGuideV = svgEl("line", {
      class: "calib-guide",
      x1: ux,
      y1: uy,
      x2: ux,
      y2: this.plot.y0,
      opacity: "0",
    });
    this.uDot = svgEl("circle", {
      cx: ux,
      cy: uy,
      r: 6,
      fill: u.color,
      stroke: "#fff",
      "stroke-width": "1.5",
      opacity: "0",
    });
    this.uRing = svgEl("circle", {
      class: "calib-unknown-ring",
      cx: ux,
      cy: uy,
      r: 11,
      opacity: "0",
    });
    this.uLab = svgEl("text", {
      class: "calib-pt-label",
      x: ux + 12,
      y: uy + 4,
      opacity: "0",
      fill: u.color,
    });
    this.uLab.textContent = `Unknown ≈ ${u.mw.toFixed(1)} kDa`;
    g.append(this.uGuideH, this.uGuideV, this.uRing, this.uDot, this.uLab);
  }
  _buildFormulas() {
    const wrap = $("detFormulas");
    const kavRows = GPC.runs
      .filter((r) => r.role !== "void")
      .map(
        (r) =>
          `<div>${r.short === "Unk" ? "<b>Unknown</b>" : r.short} · Vₑ ${r.Ve.toFixed(2)} → Kₐᵥ ${r.Kav.toFixed(3)}</div>`,
      )
      .join("");
    wrap.innerHTML = `
      <div class="formula-card" data-card="kav">
        <h3>1 · Distribution coefficient</h3>
        <div class="eq"><span class="var">K</span><sub>av</sub> =
          <span class="frac"><span class="num"><span class="var">V</span><sub>e</sub> − <span class="var">V</span><sub>o</sub></span>
          <span class="den"><span class="var">V</span><sub>t</sub> − <span class="var">V</span><sub>o</sub></span></span></div>
        <div class="res">Vₒ = ${GPC.Vo} mL&nbsp;·&nbsp;Vₜ = ${GPC.Vt} mL</div>
        <div style="font-size:.8rem;color:var(--c-ink-soft);margin-top:6px;line-height:1.6">${kavRows}</div>
      </div>
      <div class="formula-card" data-card="line">
        <h3>2 · Calibration line</h3>
        <div class="eq"><span class="var">K</span><sub>av</sub> = −${GPC.m.toFixed(3)}·log(MW) + ${GPC.c.toFixed(3)}</div>
        <div class="res">fitted through the standards (y = −mx + c)</div>
      </div>
      <div class="formula-card is-result" data-card="read">
        <h3>3 · Native molecular weight</h3>
        <div class="eq">log(MW) =
          <span class="frac"><span class="num"><span class="var">c</span> − <span class="var">K</span><sub>av</sub></span>
          <span class="den"><span class="var">m</span></span></span> = ${GPC.unknown.logMW.toFixed(3)}</div>
        <div class="res">MW (native) ≈ ${GPC.unknown.mw.toFixed(1)} kDa</div>
      </div>`;
    this.cards = {};
    wrap.querySelectorAll(".formula-card").forEach((c) => {
      this.cards[c.dataset.card] = c;
    });
  }
  render(T) {
    const { i: idx, s, p } = this.sm.localAt(T);
    // points
    const showPts = idx >= 1 ? 1 : 0;
    this.ptEls.forEach((pt) => {
      pt.dot.setAttribute("opacity", showPts);
      pt.lab.setAttribute("opacity", showPts);
    });
    // line
    if (idx >= 2) {
      const rp = idx === 2 ? p : 1;
      const x2 = lerp(this.lineA.x, this.lineB.x, rp),
        y2 = lerp(this.lineA.y, this.lineB.y, rp);
      this.lineEl.setAttribute("x2", x2.toFixed(1));
      this.lineEl.setAttribute("y2", y2.toFixed(1));
      this.lineEl.setAttribute("opacity", "1");
    } else {
      this.lineEl.setAttribute("opacity", "0");
      this.lineEl.setAttribute("x2", this.lineA.x);
      this.lineEl.setAttribute("y2", this.lineA.y);
    }
    // unknown
    const uo = idx >= 3 ? (idx === 3 ? clamp(p * 1.3, 0, 1) : 1) : 0;
    for (const e of [
      this.uGuideH,
      this.uGuideV,
      this.uDot,
      this.uRing,
      this.uLab,
    ])
      e.setAttribute("opacity", uo.toFixed(2));
    // formula cards
    this.cards.kav.classList.toggle("is-on", idx >= 0);
    this.cards.line.classList.toggle("is-on", idx >= 2);
    this.cards.read.classList.toggle("is-on", idx >= 3);
  }
  events() {}
}

/* ================================================================
 * PHASE 2.5 — COLUMN REGENERATION
 * ================================================================ */
class RegenerationPhase {
  constructor(sceneEl) {
    this.key = "regen";
    this.code = "2.5";
    this.label = "Column Regeneration";
    this.scene = sceneEl;
    this.usesChart = false;
    this.sm = new StateMachine(CONFIG.phases.regen);
    this.T = 0;
    this.el = {
      inTube: $("p4InTube"),
      outTube: $("p4OutTube"),
      liquid: $("p4Liquid"),
      adsorbed: $("p4Adsorbed"),
      waste: $("p4WasteLiquid"),
      resLiquid: $("p4ResLiquid"),
      resLabel: $("p4ResLabel"),
      fridge: $("p4Fridge"),
      cap: $("p4Cap"),
    };
    const rng = makePRNG(0xad5036);
    this.specks = [];
    for (let i = 0; i < 11; i++) {
      const c = svgEl("circle", {
        r: (1.6 + rng() * 1.8).toFixed(1),
        fill: "#d9673f",
        opacity: "0",
      });
      this.el.adsorbed.appendChild(c);
      this.specks.push({
        el: c,
        x: 398 + rng() * 54,
        y: 176 + rng() * 210,
        ph: rng() * TAU,
      });
    }
  }
  render(T) {
    const { s, p } = this.sm.localAt(T);
    let flowing = false,
      adsorbedShow = 1,
      adsorbedDrift = 0,
      wasteFrac = 0,
      wasteMix = 0,
      liquidO = 0,
      fridgeO = 0,
      capO = 0,
      resLabel = "Salt buffer",
      resFill = "url(#gradSalt)";
    switch (s.key) {
      case "wash":
        flowing = true;
        adsorbedShow = 1 - p;
        adsorbedDrift = p;
        wasteFrac = 0.5 * p;
        wasteMix = 1;
        liquidO = 0.4;
        break;
      case "equilibrate":
        flowing = true;
        adsorbedShow = 0;
        wasteFrac = 0.5 + 0.4 * p;
        wasteMix = 1 - p;
        liquidO = 0.4;
        resLabel = "Mobile phase";
        resFill = "url(#gradLiquid)";
        break;
      case "store":
        flowing = false;
        adsorbedShow = 0;
        wasteFrac = 0.9;
        wasteMix = 0;
        liquidO = 0.25;
        fridgeO = p;
        capO = clamp(p * 2, 0, 1);
        resLabel = "Mobile phase";
        resFill = "url(#gradLiquid)";
        break;
    }
    this.el.inTube.classList.toggle("is-flowing", flowing);
    this.el.outTube.classList.toggle("is-flowing", flowing);
    this.el.liquid.setAttribute("y", 158);
    this.el.liquid.setAttribute("height", 240);
    this.el.liquid.setAttribute("opacity", liquidO.toFixed(2));
    for (const sp of this.specks) {
      sp.el.setAttribute(
        "cx",
        (sp.x + Math.sin(sp.ph + T * 2) * 1.6).toFixed(1),
      );
      sp.el.setAttribute(
        "cy",
        (sp.y + adsorbedDrift * (430 - sp.y)).toFixed(1),
      );
      sp.el.setAttribute("opacity", (adsorbedShow * 0.9).toFixed(2));
    }
    const wasteH = 42 * wasteFrac;
    this.el.waste.setAttribute("height", wasteH.toFixed(1));
    this.el.waste.setAttribute("y", (512 - wasteH).toFixed(1));
    this.el.waste.setAttribute("fill", mixHex("#bfe0f2", "#d99a86", wasteMix));
    this.el.resLiquid.setAttribute("fill", resFill);
    this.el.resLabel.textContent = resLabel;
    this.el.fridge.setAttribute("opacity", fridgeO.toFixed(2));
    this.el.cap.setAttribute("opacity", capO.toFixed(2));
  }
  events() {}
}

/* ================================================================
 * PHASE 1.1 — SDS SAMPLE PREPARATION
 * ================================================================ */
class SdsPrepPhase {
  constructor(sceneEl) {
    this.key = "sdsprep";
    this.code = "1.1";
    this.label = "Sample Preparation";
    this.scene = sceneEl;
    this.usesChart = false;
    this.sm = new StateMachine(CONFIG.phases.sdsprep);
    this.T = 0;
    this.el = {
      liquid: $("sdsLiquid"),
      pellet: $("sdsPellet"),
      bubbles: $("sdsBubbles"),
      arrowA: $("sdsArrowA"),
      arrowB: $("sdsArrowB"),
      dropA: $("sdsDropA"),
      dropB: $("sdsDropB"),
      heat: $("sdsHeat"),
      spin: $("sdsSpin"),
      ready: $("sdsReady"),
    };
    const rng = makePRNG(0x5d5);
    this.bub = [];
    for (let i = 0; i < 6; i++) {
      const c = svgEl("circle", {
        r: (1.4 + rng() * 1.6).toFixed(1),
        fill: "#cfe0ff",
        opacity: "0",
      });
      this.el.bubbles.appendChild(c);
      this.bub.push({
        el: c,
        x: 398 + rng() * 34,
        off: rng(),
        speed: 0.6 + rng() * 0.7,
      });
    }
  }
  render(T) {
    const { s, p } = this.sm.localAt(T);
    let liquidFrac = 0,
      colorFrac = 0,
      pouring = false,
      heating = false,
      spinning = false,
      pelletO = 0,
      readyO = 0,
      bubbleP = -1;
    switch (s.key) {
      case "mix":
        liquidFrac = p;
        colorFrac = p;
        pouring = true;
        break;
      case "heat":
        liquidFrac = 1;
        colorFrac = 1;
        heating = true;
        bubbleP = p;
        break;
      case "spin":
        liquidFrac = 1;
        colorFrac = 1;
        spinning = true;
        pelletO = clamp((p - 0.4) / 0.6, 0, 1);
        break;
      case "ready":
        liquidFrac = 1;
        colorFrac = 1;
        pelletO = 1;
        readyO = clamp(p * 1.5, 0, 1);
        break;
    }
    // liquid level (bottom of cone at 398) + colour (pale protein -> blue dye)
    const topY = lerp(396, 250, liquidFrac);
    this.el.liquid.setAttribute("y", topY.toFixed(1));
    this.el.liquid.setAttribute("height", (398 - topY).toFixed(1));
    this.el.liquid.setAttribute(
      "fill",
      mixHex("#d9c79b", "#2b5fd0", colorFrac),
    );
    // reagent pour arrows + travelling drops
    this.el.arrowA.setAttribute("opacity", pouring ? "1" : "0");
    this.el.arrowB.setAttribute("opacity", pouring ? "1" : "0");
    if (pouring) {
      const tA = (T * 0.9) % 1,
        tB = (T * 0.9 + 0.5) % 1;
      this.el.dropA.setAttribute("cx", lerp(300, 396, tA).toFixed(1));
      this.el.dropA.setAttribute("cy", lerp(140, 236, tA).toFixed(1));
      this.el.dropA.setAttribute("r", "3.5");
      this.el.dropB.setAttribute("cx", lerp(530, 434, tB).toFixed(1));
      this.el.dropB.setAttribute("cy", lerp(140, 236, tB).toFixed(1));
      this.el.dropB.setAttribute("r", "3.5");
    } else {
      this.el.dropA.setAttribute("r", "0");
      this.el.dropB.setAttribute("r", "0");
    }
    // heat block glow + waves
    this.el.heat.setAttribute("opacity", heating ? "1" : "0");
    this.el.heat.classList.toggle("is-heating", heating);
    // centrifuge spin arcs
    this.el.spin.setAttribute("opacity", spinning ? "1" : "0");
    this.el.spin.classList.toggle("is-spinning", spinning);
    // bubbles (during heat)
    for (const b of this.bub) {
      if (bubbleP < 0) {
        b.el.setAttribute("opacity", "0");
        continue;
      }
      const prog = (bubbleP * b.speed + b.off) % 1;
      b.el.setAttribute("cx", b.x.toFixed(1));
      b.el.setAttribute("cy", lerp(390, topY + 6, prog).toFixed(1));
      b.el.setAttribute("opacity", "0.7");
    }
    // pellet + ready
    this.el.pellet.setAttribute("opacity", pelletO.toFixed(2));
    this.el.ready.setAttribute("opacity", readyO.toFixed(2));
  }
  events() {}
}

/* ================================================================
 * PHASE 1.2 — GEL RUN (SDS-PAGE)
 * ================================================================ */
class SdsRunPhase {
  constructor(sceneEl) {
    this.key = "sdsrun";
    this.code = "1.2";
    this.label = "Run on SDS-PAGE";
    this.scene = sceneEl;
    this.usesChart = false;
    this.sm = new StateMachine(CONFIG.phases.sdsrun);
    this.T = 0;
    this.gel = {
      wellY: 96,
      bottomY: 424,
      markerX: 106,
      sampleX: 226,
      bandW: 48,
    };
    this._buildGel();
    this.el = {
      dye: $("gelDye"),
      voltBadge: $("gelVoltBadge"),
      voltText: $("gelVoltText"),
    };
  }
  _buildGel() {
    const g = this.gel,
      mg = $("gelMarker"),
      sg = $("gelSample");
    this.markerBands = SDS.markers.map((m) => {
      const rect = svgEl("rect", {
        class: "gel-band",
        x: g.markerX,
        y: g.wellY,
        width: g.bandW,
        height: 5,
        rx: 1.5,
        fill: "#243f8f",
        opacity: "0",
      });
      const lab = svgEl("text", {
        class: "gel-band-label",
        x: g.markerX - 6,
        y: g.wellY,
        "text-anchor": "end",
        opacity: "0",
      });
      lab.textContent = m.mw;
      mg.append(rect, lab);
      return { rect, lab, Rf: m.Rf };
    });
    this.sampleBand = svgEl("rect", {
      class: "gel-band",
      x: g.sampleX,
      y: g.wellY,
      width: g.bandW,
      height: 6,
      rx: 1.5,
      fill: "#243f8f",
      opacity: "0",
    });
    sg.append(this.sampleBand);
  }
  render(T) {
    const { s, p } = this.sm.localAt(T);
    const g = this.gel,
      span = g.bottomY - g.wellY;
    let runP = 0,
      bandsO = 0,
      voltO = 0;
    switch (s.key) {
      case "load":
        runP = 0;
        bandsO = clamp(p, 0, 1);
        break;
      case "run":
        runP = p;
        bandsO = 1;
        voltO = 1;
        break;
      case "stain":
        runP = 1;
        bandsO = 1;
        voltO = 0;
        break;
    }
    this.markerBands.forEach((b) => {
      const y = g.wellY + b.Rf * runP * span;
      b.rect.setAttribute("y", y.toFixed(1));
      b.rect.setAttribute("opacity", bandsO.toFixed(2));
      b.lab.setAttribute("y", (y + 4).toFixed(1));
      b.lab.setAttribute(
        "opacity",
        (bandsO * (runP > 0.15 ? 1 : 0)).toFixed(2),
      );
    });
    const uy = g.wellY + SDS.unknown.Rf * runP * span;
    this.sampleBand.setAttribute("y", uy.toFixed(1));
    this.sampleBand.setAttribute("opacity", bandsO.toFixed(2));
    const dyeY = g.wellY + runP * span;
    this.el.dye.setAttribute("y1", dyeY.toFixed(1));
    this.el.dye.setAttribute("y2", dyeY.toFixed(1));
    this.el.voltBadge.setAttribute("opacity", voltO);
    this.el.voltText.setAttribute("opacity", voltO);
  }
  events() {}
}

/* ================================================================
 * PHASE 1.3 — DETERMINATION OF SUBUNIT MW (Rf curve + formulas)
 * ================================================================ */
class SdsDeterminePhase {
  constructor(sceneEl) {
    this.key = "sdsdetermine";
    this.code = "1.3";
    this.label = "Determination of MW";
    this.scene = sceneEl;
    this.usesChart = false;
    this.sm = new StateMachine(CONFIG.phases.sdsdetermine);
    this.T = 0;
    // graph: x = log(MW), y = Rf
    this.plot = {
      x0: 70,
      x1: 330,
      y0: 380,
      y1: 60,
      L0: 1.0,
      L1: 2.2,
      R0: 0,
      R1: 1,
    };
    this._buildGraph();
    this._buildFormulas();
  }
  _xOf(lm) {
    const p = this.plot;
    return lerp(p.x0, p.x1, invLerp(p.L0, p.L1, lm));
  }
  _yOf(rf) {
    const p = this.plot;
    return lerp(p.y0, p.y1, invLerp(p.R0, p.R1, rf));
  }
  _star(cx, cy, rO, rI) {
    let d = "";
    for (let i = 0; i < 10; i++) {
      const r = i % 2 ? rI : rO,
        a = -Math.PI / 2 + (i * Math.PI) / 5;
      d +=
        (i ? "L" : "M") +
        (cx + r * Math.cos(a)).toFixed(1) +
        "," +
        (cy + r * Math.sin(a)).toFixed(1) +
        " ";
    }
    return d + "Z";
  }
  _buildGraph() {
    const p = this.plot,
      grid = $("rfGrid");
    for (const lm of [1.0, 1.4, 1.8, 2.2]) {
      const x = this._xOf(lm);
      grid.appendChild(
        svgEl("line", {
          class: "calib-grid-line",
          x1: x,
          y1: p.y1,
          x2: x,
          y2: p.y0,
        }),
      );
      const t = svgEl("text", {
        class: "calib-tick",
        x,
        y: p.y0 + 16,
        "text-anchor": "middle",
      });
      t.textContent = lm.toFixed(1);
      grid.appendChild(t);
    }
    for (let k = 0; k <= 4; k++) {
      const rf = k / 4,
        y = this._yOf(rf);
      const t = svgEl("text", {
        class: "calib-tick",
        x: p.x0 - 8,
        y: y + 4,
        "text-anchor": "end",
      });
      t.textContent = rf.toFixed(2);
      grid.appendChild(t);
    }
    this.pts = SDS.markers.map((m) => {
      const dot = svgEl("circle", {
        class: "calib-pt",
        cx: this._xOf(m.logMW),
        cy: this._yOf(m.Rf),
        r: 5,
        fill: m.color,
        opacity: "0",
      });
      const lab = svgEl("text", {
        class: "calib-pt-label",
        x: this._xOf(m.logMW) + 8,
        y: this._yOf(m.Rf) - 6,
        opacity: "0",
      });
      lab.textContent = `${m.mw}`;
      $("rfPoints").append(dot, lab);
      return { dot, lab };
    });
    // calibration line (Rf as a function of log MW), clamped to the plot domain
    const rfAtLm = (lm) => (lm - SDS.intercept) / SDS.slope;
    const lmAtRf = (rf) => SDS.slope * rf + SDS.intercept;
    this.la = { lm: p.L0, rf: clamp(rfAtLm(p.L0), p.R0, p.R1) };
    this.lb = { rf: p.R0, lm: clamp(lmAtRf(p.R0), p.L0, p.L1) };
    this.la.x = this._xOf(this.la.lm);
    this.la.y = this._yOf(this.la.rf);
    this.lb.x = this._xOf(this.lb.lm);
    this.lb.y = this._yOf(this.lb.rf);
    this.lineEl = svgEl("line", {
      class: "calib-line",
      x1: this.la.x,
      y1: this.la.y,
      x2: this.la.x,
      y2: this.la.y,
      opacity: "0",
    });
    $("rfLine").appendChild(this.lineEl);
    // unknown — drawn as a star
    const u = SDS.unknown,
      ux = this._xOf(u.logMW),
      uy = this._yOf(u.Rf),
      g = $("rfUnknown");
    this.uGuideH = svgEl("line", {
      class: "rf-guide",
      x1: p.x0,
      y1: uy,
      x2: ux,
      y2: uy,
      opacity: "0",
    });
    this.uGuideV = svgEl("line", {
      class: "rf-guide",
      x1: ux,
      y1: uy,
      x2: ux,
      y2: p.y0,
      opacity: "0",
    });
    this.uStar = svgEl("path", {
      class: "calib-star",
      d: this._star(ux, uy, 10, 4.4),
      fill: u.color,
      opacity: "0",
    });
    this.uLab = svgEl("text", {
      class: "calib-pt-label",
      x: ux + 14,
      y: uy - 12,
      fill: u.color,
      "font-weight": "700",
      opacity: "0",
    });
    this.uLab.textContent = `Unknown ≈ ${u.mw.toFixed(1)} kDa`;
    g.append(this.uGuideH, this.uGuideV, this.uStar, this.uLab);
  }
  _buildFormulas() {
    const u = SDS.unknown;
    $("sdsFormula").innerHTML = `
      <div class="formula-card" data-card="rf">
        <h3>1 · Relative mobility</h3>
        <div class="eq"><span class="var">R</span><sub>f</sub> =
          <span class="frac"><span class="num">distance moved by protein</span>
          <span class="den">distance moved by dye front</span></span></div>
        <div class="res">Unknown: R<sub>f</sub> = ${u.Rf.toFixed(3)}</div>
      </div>
      <div class="formula-card is-result" data-card="mw">
        <h3>2 · Subunit molecular weight</h3>
        <div class="eq">log(MW) = (slope)·R<sub>f</sub> + intercept</div>
        <div class="eq">= (−${SDS.m.toFixed(3)})(${u.Rf.toFixed(3)}) + ${SDS.c.toFixed(3)} = ${u.logMW.toFixed(3)}</div>
        <div class="res">MW (subunit) ≈ ${u.mw.toFixed(1)} kDa</div>
      </div>`;
    this.cards = {};
    $("sdsFormula")
      .querySelectorAll(".formula-card")
      .forEach((c) => {
        this.cards[c.dataset.card] = c;
      });
  }
  render(T) {
    const { i: idx, s, p } = this.sm.localAt(T);
    // standard points (from the 'rf' stage on)
    const ptsO = idx >= 0 ? 1 : 0;
    this.pts.forEach((pt) => {
      pt.dot.setAttribute("opacity", ptsO);
      pt.lab.setAttribute("opacity", ptsO);
    });
    // calibration line (draw during 'line', full after)
    if (idx >= 1) {
      const lp = idx === 1 ? p : 1;
      this.lineEl.setAttribute("opacity", "1");
      this.lineEl.setAttribute("x2", lerp(this.la.x, this.lb.x, lp).toFixed(1));
      this.lineEl.setAttribute("y2", lerp(this.la.y, this.lb.y, lp).toFixed(1));
    } else {
      this.lineEl.setAttribute("opacity", "0");
      this.lineEl.setAttribute("x2", this.la.x);
      this.lineEl.setAttribute("y2", this.la.y);
    }
    // unknown star (during 'mw')
    const uo = idx >= 2 ? (idx === 2 ? clamp(p * 1.3, 0, 1) : 1) : 0;
    for (const e of [this.uGuideH, this.uGuideV, this.uStar, this.uLab])
      e.setAttribute("opacity", uo.toFixed(2));
    // formula cards
    this.cards.rf.classList.toggle("is-on", idx >= 0);
    this.cards.mw.classList.toggle("is-on", idx >= 2);
  }
  events() {}
}

/* ================================================================
 * UI CONTROLLER
 * ================================================================ */
class UIController {
  constructor() {
    this.dom = {
      progressFill: $("progressFill"),
      progressRoot: $("progressRoot"),
      ticks: $("progressTicks"),
      progress: document.querySelector(".progress"),
      instruction: $("instructionText"),
      hint: $("actionHint"),
      log: $("observationLog"),
      tooltip: $("tooltip"),
      sceneStack: document.querySelector(".scene-stack"),
      chartPanel: $("chartPanel"),
      chartCaption: $("chartCaption"),
      chartLegend: $("chartLegend"),
      stepTabs: Array.from(document.querySelectorAll(".step-tab")),
      subNav: $("subNav"),
      subTabs: Array.from(document.querySelectorAll(".sub-tab")),
      controls: document.querySelector(".controls"),
      btnStart: $("btnStart"),
      btnPause: $("btnPause"),
      btnPrev: $("btnPrev"),
      btnNext: $("btnNext"),
      btnReset: $("btnReset"),
      speed: $("speed"),
      speedOut: $("speedOut"),
      cTrace: $("chartTrace"),
      cLarge: $("chartLargeArea"),
      cNormal: $("chartNormalArea"),
      cSmall: $("chartSmallArea"),
      runsG: $("chartRuns"),
    };
    this._loggedKeys = new Set();
    this._lastT = 0;
    this._drawChartGrid();
  }
  _drawChartGrid() {
    const { x0, x1, yBase, yTop, gridRows } = CONFIG.chart;
    let g = "";
    for (let r = 0; r <= gridRows; r++) {
      const y = lerp(yBase, yTop, r / gridRows);
      g += `<line class="chart-grid-line" x1="${x0}" y1="${y}" x2="${x1}" y2="${y}"/>`;
    }
    $("chartGrid").innerHTML = g;
  }
  buildTicks(sm) {
    const frag = document.createDocumentFragment();
    sm.stages.forEach((s, i) => {
      const li = document.createElement("li");
      li.style.left = `${(s.startT / sm.totalT) * 100}%`;
      li.title = s.name;
      li.dataset.index = i;
      frag.appendChild(li);
    });
    this.dom.ticks.innerHTML = "";
    this.dom.ticks.appendChild(frag);
  }
  setStageInfo(stage) {
    this.dom.instruction.innerHTML = stage.instruction;
    Array.from(this.dom.ticks.children).forEach((li, i) => {
      li.classList.toggle("is-done", i < stage._idx);
      li.classList.toggle("is-active", i === stage._idx);
    });
  }
  setStaticInfo(instruction) {
    this.dom.instruction.innerHTML = instruction;
  }
  setProgress(T, sm) {
    const pct = clamp(T / sm.totalT, 0, 1) * 100;
    this.dom.progressFill.style.width = `${pct}%`;
    this.dom.progressRoot.setAttribute("aria-valuenow", Math.round(pct));
  }
  setStepTabs(active, doneSet) {
    this.dom.stepTabs.forEach((t, i) => {
      t.classList.toggle("is-active", i === active);
      t.classList.toggle("is-done", doneSet.has(i) && i !== active);
    });
  }
  setSubTabs(active, doneSet) {
    this.dom.subTabs.forEach((t, i) => {
      t.classList.toggle("is-active", i === active);
      t.classList.toggle("is-done", doneSet.has(i) && i !== active);
    });
  }
  buildSubNav(subs) {
    this.dom.subNav.innerHTML = subs
      .map(
        (ph, i) =>
          `<button class="sub-tab" type="button" data-sub="${i}"><b>${ph.code}</b> ${ph.label}</button>`,
      )
      .join("");
    this.dom.subTabs = Array.from(this.dom.subNav.querySelectorAll(".sub-tab"));
  }
  showSubNav(show) {
    this.dom.subNav.hidden = !show;
  }
  showProgress(show) {
    this.dom.progress.classList.toggle("is-hidden", !show);
  }
  setTransportEnabled(on) {
    this.dom.controls.classList.toggle("is-static", !on);
  }
  setChart(phase) {
    if (!phase || !phase.usesChart) {
      this.dom.chartPanel.hidden = true;
      return;
    }
    this.dom.chartPanel.hidden = false;
    this.dom.chartCaption.textContent = phase.chartCaption || "Detector trace";
    this.dom.chartLegend.innerHTML = phase.chartLegend || "";
  }
  clearChart() {
    for (const el of [
      this.dom.cTrace,
      this.dom.cLarge,
      this.dom.cNormal,
      this.dom.cSmall,
    ])
      el.setAttribute("d", "");
    this.dom.runsG
      .querySelectorAll("path")
      .forEach((p) => p.setAttribute("d", ""));
    this.dom.runsG
      .querySelectorAll("text")
      .forEach((t) => t.setAttribute("opacity", "0"));
  }
  setHint(stage) {
    const hint = stage && stage.hint;
    document
      .querySelectorAll(".apparatus.is-target")
      .forEach((el) => el.classList.remove("is-target"));
    if (!hint) {
      this.dom.hint.hidden = true;
      return;
    }
    this.dom.hint.hidden = false;
    this.dom.hint.textContent = hint.text;
    if (hint.target) {
      const el = $(hint.target);
      if (el) el.classList.add("is-target");
    }
  }
  suggestButton(name) {
    [this.dom.btnStart, this.dom.btnNext, this.dom.btnReset].forEach((b) =>
      b.classList.remove("is-suggested"),
    );
    if (name && this.dom[name]) this.dom[name].classList.add("is-suggested");
  }
  updateTransport({ playing, T, total }) {
    const completed = T >= total;
    this.dom.btnStart.disabled = playing;
    this.dom.btnPause.disabled = !playing;
    this.dom.btnStart.innerHTML = completed
      ? "&#9654; Replay"
      : T > 0
        ? "&#9654; Resume"
        : "&#9654; Start";
  }
  log(message, kind = "") {
    const li = document.createElement("li");
    if (kind) li.classList.add(kind === "warn" ? "is-warn" : "is-ok");
    li.innerHTML = message;
    this.dom.log.append(li);
    while (this.dom.log.children.length > 22) this.dom.log.firstChild.remove();
    this.dom.log.scrollTop = this.dom.log.scrollHeight;
  }
  logOnce(key, message, kind) {
    if (this._loggedKeys.has(key)) return;
    this._loggedKeys.add(key);
    this.log(message, kind);
  }
  resetLog() {
    this._loggedKeys.clear();
    this.dom.log.innerHTML = "";
  }
  showTooltip(el) {
    const title = el.dataset.label || "",
      info = el.dataset.info || "";
    if (!title && !info) return;
    const tip = this.dom.tooltip;
    tip.innerHTML = `<strong>${title}</strong>${info}`;
    tip.hidden = false;
    const box = el.getBoundingClientRect(),
      ref = this.dom.sceneStack.getBoundingClientRect();
    tip.style.left = `${clamp(box.left + box.width / 2 - ref.left, 70, ref.width - 70)}px`;
    tip.style.top = `${Math.max(54, box.top - ref.top - 8)}px`;
  }
  hideTooltip() {
    this.dom.tooltip.hidden = true;
  }
}

/* ================================================================
 * LAB SIMULATION — two-level orchestrator
 * ================================================================ */
class LabSimulation {
  constructor() {
    this.ui = new UIController();
    const sceneOf = (v) => document.querySelector(`.scene[data-view="${v}"]`);
    this.scenes = {
      sdsprep: sceneOf("sdsprep"),
      sdsrun: sceneOf("sdsrun"),
      sdsdetermine: sceneOf("sdsdetermine"),
      packing: sceneOf("packing"),
      prep: sceneOf("prep"),
      running: sceneOf("running"),
      determination: sceneOf("determination"),
      regen: sceneOf("regen"),
      oligo: sceneOf("oligo"),
    };
    // grouped steps (0 = SDS-PAGE, 1 = GPC); step 2 (oligomeric) is static
    this.groups = {
      0: [
        new SdsPrepPhase(this.scenes.sdsprep),
        new SdsRunPhase(this.scenes.sdsrun),
        new SdsDeterminePhase(this.scenes.sdsdetermine),
      ],
      1: [
        new ColumnPackingPhase(this.scenes.packing),
        new SamplePrepPhase(this.scenes.prep),
        new RunningPhase(this.scenes.running),
        new DeterminationPhase(this.scenes.determination),
      ],
    };
    this._buildOligo();

    this.activeStep = 0;
    this.subState = { 0: 0, 1: 0 };
    this.speed = CONFIG.transport.defaultSpeed;
    this.playing = false;
    this.doneSubs = { 0: new Set(), 1: new Set() };
    this.doneSteps = new Set();
    this._lastFrame = 0;
    this._stageIndex = -1;

    this._bindControls();
    this._bindNav();
    this._bindInteraction();
    this.selectStep(0);
  }

  get subs() {
    return this.groups[this.activeStep] || null;
  }
  get activeSub() {
    return this.subState[this.activeStep] || 0;
  }
  set activeSub(v) {
    this.subState[this.activeStep] = v;
  }
  get phase() {
    return this.subs ? this.subs[this.activeSub] : null;
  }

  /* ---------- step / sub selection ---------- */
  _showScene(view) {
    for (const k in this.scenes) this.scenes[k].hidden = k !== view;
  }

  selectStep(step) {
    this.pause();
    this.activeStep = step;
    this._applyTheme(step);
    this.ui.setStepTabs(step, this.doneSteps);
    const subs = this.groups[step];

    if (subs) {
      // grouped step (SDS or GPC)
      this.ui.buildSubNav(subs);
      this.ui.showSubNav(true);
      this.ui.showProgress(true);
      this.ui.setTransportEnabled(true);
      this.selectSub(this.activeSub, true);
      return;
    }
    // static step (Oligomeric status)
    this.ui.showSubNav(false);
    this.ui.showProgress(false);
    this.ui.setTransportEnabled(false);
    this.ui.setChart(null);
    this.ui.setHint(null);
    this.ui.suggestButton(null);
    this.ui.resetLog();
    this._showScene("oligo");
    this.doneSteps.add(2);
    this.ui.setStaticInfo(
      "The oligomeric number <strong>n</strong> is the native MW (gel filtration) divided by the subunit MW (SDS-PAGE).",
    );
    this.ui.log(
      "Step 3 — comparing native and subunit molecular weights.",
      "ok",
    );
  }

  /* colour theme per step: 1 blue · 2 red · 3 green */
  _applyTheme(step) {
    const THEME = [
      { color: "#2f6fdb", faint: "#e8f0fc" }, // Step 1 — blue
      { color: "#e0533b", faint: "#fdeee9" }, // Step 2 — red
      { color: "#2faf6b", faint: "#e8f7ef" }, // Step 3 — green
    ][step] || { color: "#2f9bdb", faint: "#eaf5fd" };
    document.body.style.setProperty("--step-color", THEME.color);
    document.body.style.setProperty("--step-faint", THEME.faint);
  }

  selectSub(sub, keepStep = false) {
    if (!keepStep) this.pause();
    this.activeSub = sub;
    const ph = this.phase;
    this._showScene(ph.key);
    this.ui.setSubTabs(sub, this.doneSubs[this.activeStep]);
    this.ui.clearChart();
    this.ui.setChart(ph);
    this.ui.buildTicks(ph.sm);
    this.ui.resetLog();
    this.ui.log(
      `<strong>Step ${ph.code} — ${ph.label}.</strong> Press Start to run this step.`,
      "ok",
    );
    this._stageIndex = -1;
    this.render();
    this._refreshChrome();
    // the gel run (1.2) starts on entry — no manual "Run the gel" click needed
    if (ph.key === "sdsrun" && ph.T === 0) this.play();
  }

  /* ---------- transport ---------- */
  play() {
    if (!this.subs) return;
    const ph = this.phase;
    if (ph.T >= ph.sm.totalT) ph.T = 0;
    // stepwise phases (e.g. running samples) advance one sub-step per click
    this._stopAt = ph.stepwise ? ph.sm.stages[ph.sm.indexAt(ph.T)].endT : null;
    this.playing = true;
    this._lastFrame = performance.now();
    this._loop();
    this._refreshChrome();
  }
  pause() {
    this.playing = false;
    this._stopAt = null;
    this._refreshChrome();
    if (this.phase) this.render();
  }
  togglePlay() {
    this.playing ? this.pause() : this.play();
  }
  setT(T) {
    const ph = this.phase;
    if (!ph) return;
    this._stopAt = null;
    ph.T = clamp(T, 0, ph.sm.totalT);
    this.render();
    this._refreshChrome();
  }

  nextStep() {
    if (!this.subs) return;
    const ph = this.phase,
      i = ph.sm.indexAt(ph.T);
    if (i < ph.sm.stages.length - 1) {
      this.pause();
      this.setT(ph.sm.stages[i + 1].startT);
    } else if (this.activeSub < this.subs.length - 1)
      this.selectSub(this.activeSub + 1);
    else this.selectStep(this.activeStep + 1);
  }
  prevStep() {
    if (!this.subs) return;
    const ph = this.phase,
      i = ph.sm.indexAt(ph.T),
      start = ph.sm.stages[i].startT;
    if (ph.T - start > 0.15) {
      this.pause();
      this.setT(start);
    } else if (i > 0) {
      this.pause();
      this.setT(ph.sm.stages[i - 1].startT);
    } else if (this.activeSub > 0) {
      const prev = this.subs[this.activeSub - 1];
      this.selectSub(this.activeSub - 1);
      this.setT(prev.sm.stages[prev.sm.stages.length - 1].startT);
    } else if (this.activeStep > 0 && this.groups[this.activeStep - 1]) {
      const g = this.groups[this.activeStep - 1];
      this.selectStep(this.activeStep - 1);
      this.selectSub(g.length - 1);
    }
  }
  reset() {
    if (!this.subs) return;
    this.pause();
    this.phase.T = 0;
    this._stageIndex = -1;
    this.ui.resetLog();
    this.ui.log(
      `<strong>Step ${this.phase.code} — ${this.phase.label}.</strong> Reset.`,
      "ok",
    );
    this.render();
    this._refreshChrome();
  }
  setSpeed(v) {
    this.speed = clamp(v, CONFIG.transport.minSpeed, CONFIG.transport.maxSpeed);
    this.ui.dom.speedOut.textContent = `${this.speed.toFixed(2).replace(/0$/, "")}×`;
  }

  _loop() {
    if (!this.playing) return;
    const now = performance.now(),
      dt = Math.min(0.05, (now - this._lastFrame) / 1000);
    this._lastFrame = now;
    const ph = this.phase;
    const limit =
      this._stopAt != null
        ? Math.min(this._stopAt, ph.sm.totalT)
        : ph.sm.totalT;
    ph.T += dt * this.speed;
    if (ph.T >= limit) {
      ph.T = limit;
      this.playing = false;
      this._stopAt = null;
      const atLast =
        ph.stepwise && ph.sm.indexAt(ph.T) === ph.sm.stages.length - 1;
      if (
        (ph.T >= ph.sm.totalT || atLast) &&
        !this.doneSubs[this.activeStep].has(this.activeSub)
      )
        this._onSubComplete();
    }
    this.render();
    if (this.playing) requestAnimationFrame(() => this._loop());
    else this._refreshChrome();
  }
  _onSubComplete() {
    const done = this.doneSubs[this.activeStep];
    done.add(this.activeSub);
    this.ui.setSubTabs(this.activeSub, done);
    if (done.size === this.subs.length) {
      this.doneSteps.add(this.activeStep);
      this.ui.setStepTabs(this.activeStep, this.doneSteps);
    }
    this.ui.log(`<strong>${this.phase.label} complete.</strong>`, "ok");
  }

  /* ---------- render ---------- */
  render() {
    const ph = this.phase;
    if (!ph) return;
    const T = ph.T;
    this.ui._lastT = T;
    ph.render(T, this.playing);
    this.ui.setProgress(T, ph.sm);
    const idx = ph.sm.indexAt(T);
    if (idx !== this._stageIndex) {
      this._stageIndex = idx;
      const stage = ph.sm.stages[idx];
      stage._idx = idx;
      this.ui.setStageInfo(stage);
      this.ui.setHint(stage);
      this.ui.suggestButton(
        stage.hint && stage.hint.action === "play" ? "btnStart" : null,
      );
      this.ui.logOnce(
        `${ph.key}-${stage.key}`,
        `<strong>${stage.name}.</strong> ${stage.note || ""}`,
      );
    }
    if (ph.events) ph.events(T, (k, m, kind) => this.ui.logOnce(k, m, kind));
  }
  _refreshChrome() {
    const ph = this.phase;
    if (!ph) {
      this.ui.updateTransport({ playing: false, T: 0, total: 1 });
      return;
    }
    this.ui.updateTransport({
      playing: this.playing,
      T: ph.T,
      total: ph.sm.totalT,
    });
    if (ph.T >= ph.sm.totalT) {
      const last = ph.sm.stages[ph.sm.stages.length - 1];
      last._idx = ph.sm.stages.length - 1;
      this.ui.setHint(last);
      this.ui.suggestButton("btnNext");
    }
  }

  /* ---------- Step 3 formula ---------- */
  _buildOligo() {
    const native = GPC.unknown.mw,
      subunit = SDS.unknown.mw;
    const n = native / subunit,
      rounded = Math.round(n);
    const name =
      rounded === 1
        ? "MONOMER"
        : rounded === 2
          ? "DIMER"
          : rounded === 3
            ? "TRIMER"
            : `${rounded}-MER`;
    $("oligoStage").innerHTML = `
      <div class="big-formula">
        <h2>Oligomeric status of the protein</h2>
        <p class="sub">How many subunits make up the native protein?</p>
        <div class="big-eq"><span class="var">n</span> =
          <span class="frac"><span class="num">MW&nbsp;(native, gel filtration)</span>
          <span class="den">MW&nbsp;(subunit, SDS-PAGE)</span></span></div>
        <div class="plug-row">
          <div class="plug"><b>${native.toFixed(1)} kDa</b><small>native · Step 2 (GPC)</small></div>
          <div class="plug"><b>${subunit.toFixed(1)} kDa</b><small>subunit · Step 1 (SDS-PAGE)</small></div>
        </div>
        <div class="big-eq"><span class="var">n</span> =
          <span class="frac"><span class="num">${native.toFixed(1)}</span><span class="den">${subunit.toFixed(1)}</span></span>
          = ${n.toFixed(2)} ≈ ${rounded}</div>
        <div class="verdict">The protein is a ${name}
          <small>n ≈ ${rounded} → ${rounded} polypeptide subunit${rounded === 1 ? "" : "s"} in the native protein</small></div>
      </div>`;
  }

  /* ---------- events ---------- */
  _bindControls() {
    const d = this.ui.dom;
    d.btnStart.addEventListener("click", () => this.play());
    d.btnPause.addEventListener("click", () => this.pause());
    d.btnPrev.addEventListener("click", () => this.prevStep());
    d.btnNext.addEventListener("click", () => this.nextStep());
    d.btnReset.addEventListener("click", () => this.reset());
    d.speed.addEventListener("input", (e) =>
      this.setSpeed(parseFloat(e.target.value)),
    );
    d.ticks.addEventListener("click", (e) => {
      const li = e.target.closest("li");
      if (!li || !this.phase) return;
      this.pause();
      this.setT(this.phase.sm.stages[+li.dataset.index].startT);
    });
    window.addEventListener("keydown", (e) => {
      if (["INPUT", "TEXTAREA"].includes(document.activeElement.tagName))
        return;
      if (e.code === "Space") {
        e.preventDefault();
        this.togglePlay();
      } else if (e.code === "ArrowRight") this.nextStep();
      else if (e.code === "ArrowLeft") this.prevStep();
      else if (e.key.toLowerCase() === "r") this.reset();
    });
  }
  _bindNav() {
    this.ui.dom.stepTabs.forEach((t) =>
      t.addEventListener("click", () => this.selectStep(+t.dataset.step)),
    );
    this.ui.dom.subNav.addEventListener("click", (e) => {
      const b = e.target.closest(".sub-tab");
      if (!b || !this.subs) return;
      this.selectSub(+b.dataset.sub);
    });
  }
  _bindInteraction() {
    Array.from(document.querySelectorAll(".apparatus")).forEach((el) => {
      el.addEventListener("click", () => this._onApparatusClick(el));
      el.addEventListener("mouseenter", () => this.ui.showTooltip(el));
      el.addEventListener("mouseleave", () => this.ui.hideTooltip());
      el.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.code === "Space") {
          e.preventDefault();
          this._onApparatusClick(el);
        }
      });
    });
    this.ui.dom.sceneStack.addEventListener("click", (e) => {
      if (!e.target.closest(".apparatus")) this.ui.hideTooltip();
    });
    this.ui.dom.hint.addEventListener("click", () => {
      if (this.phase) this._performHint(this.phase.sm.stageAt(this.phase.T));
    });
  }
  _onApparatusClick(el) {
    this.ui.showTooltip(el);
    if (!this.phase) return;
    const stage = this.phase.sm.stageAt(this.phase.T),
      hint = stage.hint;
    if (hint && hint.target === el.id && !this.playing)
      this._performHint(stage);
  }
  _performHint(stage) {
    if (!stage || !stage.hint) return;
    switch (stage.hint.action) {
      case "play":
        this.play();
        break;
      case "next":
        this.nextStep();
        break;
      case "nextsub":
        if (this.activeSub < this.subs.length - 1)
          this.selectSub(this.activeSub + 1);
        else this.selectStep(this.activeStep + 1);
        break;
      case "nextstep":
        this.selectStep(this.activeStep + 1);
        break;
      case "reset":
        this.reset();
        break;
    }
  }
}

/* ================================================================
 * BOOTSTRAP
 * ================================================================ */
window.addEventListener("DOMContentLoaded", () => {
  window.lab = new LabSimulation();
});
