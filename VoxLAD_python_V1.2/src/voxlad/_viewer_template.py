"""HTML page of the 3-D viewer (see viewer.py).  One self-contained file: no
network access, no other files.  __TITLE__ and __DATA__ are filled in by
viewer.build_html()."""

TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ - VoxLAD viewer</title>
<style>
:root { color-scheme: light dark; }
.viz-root {
  color-scheme: light;
  --surface-1: #fcfcfb; --page: #f9f9f7;
  --text-primary: #0b0b0b; --text-secondary: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --wash: rgba(11,11,11,0.05);
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) .viz-root {
    color-scheme: dark;
    --surface-1: #1a1a19; --page: #0d0d0d;
    --text-primary: #f0efec; --text-secondary: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --wash: rgba(255,255,255,0.07);
  }
}
:root[data-theme="dark"] .viz-root {
  color-scheme: dark;
  --surface-1: #1a1a19; --page: #0d0d0d;
  --text-primary: #f0efec; --text-secondary: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --wash: rgba(255,255,255,0.07);
}
html, body { height: 100%; margin: 0; }
body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; font-size: 13px; }
.viz-root { display: flex; flex-direction: column; height: 100%; background: var(--page); color: var(--text-primary); }
header { padding: 10px 14px 4px; }
h1 { font-size: 16px; font-weight: 600; margin: 0 0 2px; }
#stats { margin: 0; color: var(--text-secondary); }
#stats b { color: var(--text-primary); font-weight: 600; }
.toolbar { display: flex; flex-wrap: wrap; align-items: center; gap: 6px 16px; padding: 6px 14px 8px; color: var(--text-secondary); }
.toolbar label, .toolbar .grp { display: inline-flex; align-items: center; gap: 6px; white-space: nowrap; }
.toolbar input[type=range] { width: 92px; accent-color: #2a78d6; }
.toolbar input[type=number] { width: 58px; }
.toolbar select, .toolbar input[type=number], .toolbar button {
  font: inherit; color: var(--text-primary); background: var(--surface-1);
  border: 1px solid var(--border); border-radius: 6px; padding: 3px 7px;
}
.toolbar button { cursor: pointer; }
.toolbar button:hover { background: var(--wash); }
.toolbar button[aria-pressed=true] { border-color: var(--muted); font-weight: 600; }
.toolbar output { color: var(--text-primary); font-variant-numeric: tabular-nums; min-width: 2.6em; }
.stage { position: relative; flex: 1 1 auto; min-height: 240px; margin: 0 14px 12px; background: var(--surface-1);
  border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
canvas { position: absolute; inset: 0; width: 100%; height: 100%; display: block; touch-action: none; cursor: grab; outline-offset: -2px; }
canvas.dragging { cursor: grabbing; }
#labels { position: absolute; inset: 0; pointer-events: none; overflow: hidden; }
#labels span { position: absolute; transform: translate(-50%, -50%); color: var(--muted); font-size: 11px;
  font-variant-numeric: tabular-nums; white-space: nowrap; }
#labels span.t { color: var(--text-secondary); }
#legend { position: absolute; left: 10px; bottom: 10px; max-width: 46%; padding: 8px 10px; background: var(--surface-1);
  border: 1px solid var(--border); border-radius: 6px; color: var(--text-secondary); }
#legend .ttl { color: var(--text-primary); font-weight: 600; margin-bottom: 5px; }
#legend .bar { width: 190px; height: 10px; border-radius: 2px; }
#legend .ends { display: flex; justify-content: space-between; width: 190px; margin-top: 2px; font-variant-numeric: tabular-nums; }
#legend .row { display: flex; align-items: center; gap: 6px; margin-top: 3px; }
#legend button.row { font: inherit; color: inherit; background: none; border: 0; padding: 1px 2px; cursor: pointer; border-radius: 4px; width: 100%; text-align: left; }
#legend button.row:hover { background: var(--wash); }
#legend button.row[aria-pressed=false] { opacity: 0.45; text-decoration: line-through; }
#legend .sw { width: 11px; height: 11px; border-radius: 2px; flex: none; }
#legend .n { margin-left: auto; padding-left: 12px; font-variant-numeric: tabular-nums; }
#tooltip { position: absolute; display: none; pointer-events: none; min-width: 170px; padding: 8px 10px; background: var(--surface-1);
  border: 1px solid var(--border); border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.18); color: var(--text-secondary); }
#tooltip .v { color: var(--text-primary); font-size: 15px; font-weight: 600; }
#tooltip table { border-collapse: collapse; margin-top: 4px; }
#tooltip td { padding: 1px 0; }
#tooltip td:first-child { padding-right: 12px; }
#tooltip td:last-child { color: var(--text-primary); text-align: right; font-variant-numeric: tabular-nums; }
#help { position: absolute; right: 10px; bottom: 8px; color: var(--muted); font-size: 11px; text-align: right; }
@media (max-width: 760px) { #help { display: none; } .toolbar .grp { flex-wrap: wrap; } #legend { max-width: 80%; } }
#msg { position: absolute; inset: 0; display: none; align-items: center; justify-content: center; padding: 20px; text-align: center; color: var(--text-secondary); }
#tablewrap { margin: 0 14px 12px; max-height: 38vh; overflow: auto; background: var(--surface-1); border: 1px solid var(--border); border-radius: 8px; }
#tablewrap caption { caption-side: top; text-align: left; padding: 8px 10px; color: var(--text-secondary); }
#tablewrap table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }
#tablewrap th, #tablewrap td { padding: 3px 10px; text-align: right; border-top: 1px solid var(--grid); white-space: nowrap; }
#tablewrap th { position: sticky; top: 0; background: var(--surface-1); color: var(--text-secondary); font-weight: 600; }
#tablewrap td.l, #tablewrap th.l { text-align: left; }
</style>
</head>
<body>
<div class="viz-root" id="root">
  <header>
    <h1 id="title"></h1>
    <p id="stats"></p>
  </header>
  <div class="toolbar" role="group" aria-label="Display and filters">
    <label>Colour by
      <select id="mode">
        <option value="lad">Leaf area density</option>
        <option value="explored">Percent explored</option>
        <option value="flag">Quality flag</option>
        <option value="scan">Scan used</option>
      </select>
    </label>
    <span class="grp" id="ladgrp">LAD from <input type="range" id="ladmin" aria-label="Lowest LAD shown"><output id="ladminv"></output>
      to <input type="range" id="ladmax" aria-label="Highest LAD shown"><output id="ladmaxv"></output></span>
    <label id="cmaxlab">Colour scale up to <input type="number" id="cmax" min="0" step="0.1"></label>
    <span class="grp"><label>Cut along <select id="clipaxis"><option value="2">height (z)</option><option value="0">x</option><option value="1">y</option></select></label>
      <input type="range" id="cliplo" aria-label="Cut: lower limit"><input type="range" id="cliphi" aria-label="Cut: upper limit"><output id="clipv"></output>
      <button type="button" id="clipreset" hidden>Reset cuts</button></span>
    <span class="grp" id="flagtoggles">Show</span>
    <label><input type="checkbox" id="hidezero"> Hide LAD = 0</label>
    <label>Voxel size <input type="range" id="size" min="20" max="100" value="90"></label>
    <label><input type="checkbox" id="scanners"> Scanners</label>
    <span class="grp">View
      <button type="button" data-view="iso">3-D</button><button type="button" data-view="top">Top</button><button type="button" data-view="front">Front</button><button type="button" data-view="side">Side</button>
    </span>
    <label><input type="checkbox" id="ortho"> Parallel projection</label>
    <button type="button" id="tablebtn" aria-pressed="false">Table</button>
    <button type="button" id="png">Save image</button>
    <button type="button" id="theme">Dark</button>
  </div>
  <div class="stage" id="stage">
    <canvas id="gl" tabindex="0" aria-label="Three-dimensional view of the voxels. Drag to rotate, wheel to zoom. The Table button lists the same values."></canvas>
    <div id="labels"></div>
    <div id="legend"></div>
    <div id="help">Drag: rotate &middot; Shift-drag or right-drag: move &middot; Wheel: zoom &middot; Double-click: reset</div>
    <div id="tooltip" role="status"></div>
    <div id="msg"></div>
  </div>
  <div id="tablewrap" hidden></div>
</div>
<script type="application/json" id="voxlad-data">__DATA__</script>
<script>
(() => {
"use strict";
const D = JSON.parse(document.getElementById("voxlad-data").textContent);
const $ = id => document.getElementById(id);
const root = document.documentElement, viz = $("root"), stage = $("stage"), canvas = $("gl");

function decode(a) {
  const s = atob(a.b64), u = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) u[i] = s.charCodeAt(i);
  const T = {uint8: Uint8Array, uint16: Uint16Array, int32: Int32Array, float32: Float32Array}[a.dtype];
  return new T(u.buffer);
}
const A = {};
for (const k in D.arrays) A[k] = decode(D.arrays[k]);
const N = D.n, H = D.voxel_size, VOL = H * H * H;
const FLAG = ["ok", "low sampling", "unresolved", "occluded", "not seen"];
const FLAG_ICON = ["✓", "!", "?", "×", "–"];
const FLAG_COL = ["#0ca30c", "#fab219", "#ec835a", "#d03b3b", "#898781"];     // status palette + neutral
const CAT = {light: ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"],
             dark:  ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]};
const RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"];
const NODATA = "#898781";
const hex = h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
const RAMP_RGB = RAMP.map(hex);

// ---------------------------------------------------------------- scene extent
let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
for (let i = 0; i < N; i++) for (let a = 0; a < 3; a++) {
  const v = A.ijk[3 * i + a];
  if (v < lo[a]) lo[a] = v;
  if (v > hi[a]) hi[a] = v;
}
if (N === 0) { lo = [0, 0, 0]; hi = [0, 0, 0]; }
const wmin = [0, 1, 2].map(a => D.origin[a] + lo[a] * H), wmax = [0, 1, 2].map(a => D.origin[a] + (hi[a] + 1) * H);
const centre = [0, 1, 2].map(a => (wmin[a] + wmax[a]) / 2);
const half = [0, 1, 2].map(a => (wmax[a] - wmin[a]) / 2);
const radius = Math.hypot(half[0], half[1], half[2]) || H;
let ladMaxData = 0;
const finite = [];
for (let i = 0; i < N; i++) if (Number.isFinite(A.lad[i])) { finite.push(A.lad[i]); if (A.lad[i] > ladMaxData) ladMaxData = A.lad[i]; }
finite.sort((a, b) => a - b);
const p98 = finite.length ? finite[Math.min(finite.length - 1, Math.floor(0.98 * finite.length))] : 1;
const nice = v => { if (!(v > 0)) return 1; const e = Math.pow(10, Math.floor(Math.log10(v))); const m = v / e; return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10) * e; };
const up2 = v => { if (!(v > 0)) return 1; const e = Math.pow(10, Math.floor(Math.log10(v)) - 1); return Math.round(Math.ceil(v / e - 1e-9) * e * 1e6) / 1e6; };
const fmt = (v, d) => Number.isFinite(v) ? v.toFixed(d) : "no value";

// ---------------------------------------------------------------- state
const S = {
  mode: "lad", ladMin: 0, ladMax: ladMaxData, cmax: Math.max(up2(p98), 1e-6), clipAxis: 2,
  clip: [[lo[0], hi[0]], [lo[1], hi[1]], [lo[2], hi[2]]], flags: [true, true, true, true, true],
  scans: {}, hideZero: false, size: 0.9, ortho: false, scanners: false, table: false,
  az: -58, el: 22, dist: radius * 3.6, target: [0, 0, 0], hover: -1,
};

// ---------------------------------------------------------------- WebGL
const gl = canvas.getContext("webgl2", {antialias: true, alpha: false});
if (!gl) {
  $("msg").style.display = "flex";
  $("msg").textContent = "This viewer needs WebGL 2, which every current browser provides. The same values are in the voxel table (CSV file).";
  return;
}
function program(vs, fs) {
  const p = gl.createProgram();
  for (const [t, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
    const s = gl.createShader(t); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    gl.attachShader(p, s);
  }
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}
const cubeProg = program(`#version 300 es
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aNormal;
layout(location=2) in vec3 iCentre; layout(location=3) in vec4 iColor; layout(location=4) in vec4 iPick;
uniform mat4 uMVP; uniform float uSize; uniform int uPick; out vec4 vColor;
void main() {
  gl_Position = uMVP * vec4(iCentre + aPos * uSize, 1.0);
  if (uPick == 1) { vColor = iPick; }
  else {
    float s = aNormal.z > 0.5 ? 1.0 : (aNormal.z < -0.5 ? 0.72 : (abs(aNormal.x) > 0.5 ? 0.91 : 0.82));
    vColor = vec4(iColor.rgb * s, 1.0);
  }
}`, `#version 300 es
precision mediump float; in vec4 vColor; out vec4 o; void main() { o = vColor; }`);
const lineProg = program(`#version 300 es
layout(location=0) in vec3 aPos; uniform mat4 uMVP; void main() { gl_Position = uMVP * vec4(aPos, 1.0); }`,
`#version 300 es
precision mediump float; uniform vec4 uColor; out vec4 o; void main() { o = uColor; }`);
const U = {mvp: gl.getUniformLocation(cubeProg, "uMVP"), size: gl.getUniformLocation(cubeProg, "uSize"), pick: gl.getUniformLocation(cubeProg, "uPick"),
           lmvp: gl.getUniformLocation(lineProg, "uMVP"), lcol: gl.getUniformLocation(lineProg, "uColor")};

// unit cube, outward faces counter-clockwise
const FACES = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[-1, 0, 0], [0, 0, 1], [0, 1, 0]], [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
               [[0, -1, 0], [1, 0, 0], [0, 0, 1]], [[0, 0, 1], [1, 0, 0], [0, 1, 0]], [[0, 0, -1], [0, 1, 0], [1, 0, 0]]];
const cube = [];
for (const [n, u, v] of FACES) for (const [a, b] of [[-1, -1], [1, -1], [1, 1], [-1, -1], [1, 1], [-1, 1]])
  cube.push(...[0, 1, 2].map(k => 0.5 * (n[k] + a * u[k] + b * v[k])), ...n);
const vao = gl.createVertexArray();
gl.bindVertexArray(vao);
const cubeBuf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, cubeBuf);
gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(cube), gl.STATIC_DRAW);
gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 24, 0);
gl.enableVertexAttribArray(1); gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 24, 12);
const centreBuf = gl.createBuffer(), colourBuf = gl.createBuffer(), pickBuf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, centreBuf);
gl.enableVertexAttribArray(2); gl.vertexAttribPointer(2, 3, gl.FLOAT, false, 0, 0); gl.vertexAttribDivisor(2, 1);
gl.bindBuffer(gl.ARRAY_BUFFER, colourBuf);
gl.enableVertexAttribArray(3); gl.vertexAttribPointer(3, 4, gl.UNSIGNED_BYTE, true, 0, 0); gl.vertexAttribDivisor(3, 1);
gl.bindBuffer(gl.ARRAY_BUFFER, pickBuf);
gl.enableVertexAttribArray(4); gl.vertexAttribPointer(4, 4, gl.UNSIGNED_BYTE, true, 0, 0); gl.vertexAttribDivisor(4, 1);
const lineVao = gl.createVertexArray();
gl.bindVertexArray(lineVao);
const lineBuf = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, lineBuf);
gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
gl.bindVertexArray(null);
const pickFbo = gl.createFramebuffer(), pickTex = gl.createTexture(), pickDepth = gl.createRenderbuffer();
let pickW = 0, pickH = 0, pickDirty = true;

// ---------------------------------------------------------------- matrices (column-major)
const mul = (a, b) => { const o = new Float32Array(16);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) o[4 * c + r] = a[r] * b[4 * c] + a[4 + r] * b[4 * c + 1] + a[8 + r] * b[4 * c + 2] + a[12 + r] * b[4 * c + 3];
  return o; };
const persp = (fovy, asp, n, f) => { const t = 1 / Math.tan(fovy / 2), o = new Float32Array(16);
  o[0] = t / asp; o[5] = t; o[10] = (f + n) / (n - f); o[11] = -1; o[14] = 2 * f * n / (n - f); return o; };
const ortho = (w, h, n, f) => { const o = new Float32Array(16);
  o[0] = 1 / w; o[5] = 1 / h; o[10] = 2 / (n - f); o[14] = (f + n) / (n - f); o[15] = 1; return o; };
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const FOV = 32 * Math.PI / 180;
let MVP = new Float32Array(16), camRight = [1, 0, 0], camUp = [0, 0, 1];
function camAxes() {         // unit vectors: towards the eye, screen right, screen up
  const az = S.az * Math.PI / 180, el = S.el * Math.PI / 180;
  const dir = [Math.cos(el) * Math.cos(az), Math.cos(el) * Math.sin(az), Math.sin(el)], f = [-dir[0], -dir[1], -dir[2]];
  const up0 = Math.abs(S.el) > 89 ? [-Math.cos(az), -Math.sin(az), 0].map(v => v * Math.sign(S.el)) : [0, 0, 1];
  const r = norm(cross(f, up0));
  return {dir, f, r, u: cross(r, f)};
}
let autoFit = true;       // until the user zooms, moves or turns the scene, keep all of it in view
function fitDist() {
  const W = Math.max(stage.clientWidth, 1), Hs = Math.max(stage.clientHeight, 1), t = Math.tan(FOV / 2);
  const kx = t * (W / Hs) * Math.max(0.35, 1 - 150 / W), ky = t * Math.max(0.35, 1 - 130 / Hs);   // room left for the axis labels
  const {dir, r, u} = camAxes(), pts = [];
  for (let c = 0; c < 8; c++) pts.push([c & 1 ? half[0] : -half[0], c & 2 ? half[1] : -half[1], c & 4 ? half[2] : -half[2]]);
  if (S.scanners) for (const q of D.scans) pts.push(sub(q.position, centre));
  let d = radius * 0.2;
  for (const q of pts) { const z = S.ortho ? 0 : dot(q, dir); d = Math.max(d, z + Math.abs(dot(q, r)) / kx, z + Math.abs(dot(q, u)) / ky); }
  return d;
}
function camera() {
  const {dir, f, r, u} = camAxes();
  const eye = [0, 1, 2].map(k => S.target[k] + S.dist * dir[k]);
  camRight = r; camUp = u;
  const view = new Float32Array([r[0], u[0], -f[0], 0, r[1], u[1], -f[1], 0, r[2], u[2], -f[2], 0, -dot(r, eye), -dot(u, eye), dot(f, eye), 1]);
  const asp = canvas.width / Math.max(canvas.height, 1), n = Math.max(S.dist / 200, 1e-3), fr = S.dist * 60 + radius * 4;
  const hh = S.dist * Math.tan(FOV / 2);
  MVP = mul(S.ortho ? ortho(hh * asp, hh, -fr, fr) : persp(FOV, asp, n, fr), view);
}
function project(p) {     // scene point -> CSS pixels, or null when behind the eye
  const x = MVP[0] * p[0] + MVP[4] * p[1] + MVP[8] * p[2] + MVP[12], y = MVP[1] * p[0] + MVP[5] * p[1] + MVP[9] * p[2] + MVP[13],
        w = MVP[3] * p[0] + MVP[7] * p[1] + MVP[11] * p[2] + MVP[15];
  if (w <= 1e-6) return null;
  return [(x / w * 0.5 + 0.5) * stage.clientWidth, (1 - (y / w * 0.5 + 0.5)) * stage.clientHeight];
}

// ---------------------------------------------------------------- colours and filters
const dark = () => { const t = root.getAttribute("data-theme"); return t ? t === "dark" : matchMedia("(prefers-color-scheme: dark)").matches; };
function ramp(t) {          // sequential, one hue: light (low) to dark (high); the anchor flips on a dark surface
  if (!(t > 0)) t = 0; if (t > 1) t = 1;
  if (dark()) t = 1 - t;
  const x = t * (RAMP_RGB.length - 1), i = Math.min(Math.floor(x), RAMP_RGB.length - 2), f = x - i, a = RAMP_RGB[i], b = RAMP_RGB[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}
const OTHER = {light: "#4a4a45", dark: "#d2d1c8"};
const scanColour = s => s < 1 ? NODATA : s <= 8 ? CAT[dark() ? "dark" : "light"][s - 1] : OTHER[dark() ? "dark" : "light"];   // colour follows the scan, never its rank; scans after the eighth share one
function palette() {       // the colours of the current mode and theme, worked out once per redraw: [r, g, b] of voxel i
  const nd = hex(NODATA);
  if (S.mode === "flag") { const t = FLAG_COL.map(hex); return i => t[A.flag[i]] || nd; }
  if (S.mode === "scan") { const t = []; for (let k = 0; k < 256; k++) t.push(hex(scanColour(k))); return i => t[A.scan[i]]; }
  const lut = []; for (let k = 0; k < 256; k++) lut.push(ramp(k / 255));
  const src = S.mode === "lad" ? A.lad : A.explored, k = 255 / (S.mode === "lad" ? S.cmax : 100);
  return i => { const v = src[i] * k; return v !== v || v === Infinity || v === -Infinity ? nd : lut[v <= 0 ? 0 : v >= 255 ? 255 : Math.round(v)]; };
}
function shown(i) {
  if (!S.flags[A.flag[i]]) return false;
  if (S.scans[A.scan[i]] === false) return false;
  const l = A.lad[i];
  if (Number.isFinite(l)) { if (l < S.ladMin || l > S.ladMax || (S.hideZero && l === 0)) return false; }
  else if (S.ladMin > 0) return false;
  for (let a = 0; a < 3; a++) { const v = A.ijk[3 * i + a]; if (v < S.clip[a][0] || v > S.clip[a][1]) return false; }
  return true;
}
let nShown = 0, nDraw = 0, culled = false, vis = new Int32Array(0), occ = null;
const CULL_FROM = 120000, gdim = [0, 1, 2].map(a => hi[a] - lo[a] + 1), cells = N ? gdim[0] * gdim[1] * gdim[2] : 0;
function rebuild() {
  const idx = new Int32Array(N);
  let m = 0, area = 0, areaAll = 0;
  for (let i = 0; i < N; i++) {
    const l = A.lad[i], ok = Number.isFinite(l) && A.flag[i] <= 1;       // as the total of the summary: flags ok and low sampling
    if (ok) areaAll += l;
    if (!shown(i)) continue;
    if (ok) area += l;
    idx[m++] = i;
  }
  nShown = m; vis = idx.subarray(0, m);
  // a large array: voxels hidden behind six shown neighbours are not sent to the graphics card
  culled = m > CULL_FROM && S.size >= 0.85 && cells <= 1.2e8;
  const sx = gdim[0], sxy = gdim[0] * gdim[1];
  if (culled) {
    if (occ) occ.fill(0); else occ = new Uint8Array(cells);
    for (let k = 0; k < m; k++) { const i = vis[k]; occ[(A.ijk[3 * i] - lo[0]) + (A.ijk[3 * i + 1] - lo[1]) * sx + (A.ijk[3 * i + 2] - lo[2]) * sxy] = 1; }
  }
  const c = new Float32Array(3 * m), col = new Uint8Array(4 * m), pk = new Uint8Array(4 * m), colourOf = palette();
  let d = 0;
  for (let k = 0; k < m; k++) {
    const i = vis[k];
    if (culled) {
      const x = A.ijk[3 * i] - lo[0], y = A.ijk[3 * i + 1] - lo[1], z = A.ijk[3 * i + 2] - lo[2], q = x + y * sx + z * sxy;
      if (x > 0 && y > 0 && z > 0 && x < gdim[0] - 1 && y < gdim[1] - 1 && z < gdim[2] - 1 &&
          occ[q - 1] && occ[q + 1] && occ[q - sx] && occ[q + sx] && occ[q - sxy] && occ[q + sxy]) continue;
    }
    for (let a = 0; a < 3; a++) c[3 * d + a] = D.origin[a] + (A.ijk[3 * i + a] + 0.5) * H - centre[a];
    const rgb = colourOf(i);
    col[4 * d] = rgb[0]; col[4 * d + 1] = rgb[1]; col[4 * d + 2] = rgb[2]; col[4 * d + 3] = 255;
    const id = i + 1;
    pk[4 * d] = id & 255; pk[4 * d + 1] = (id >> 8) & 255; pk[4 * d + 2] = (id >> 16) & 255; pk[4 * d + 3] = 255;
    d++;
  }
  nDraw = d;
  gl.bindBuffer(gl.ARRAY_BUFFER, centreBuf); gl.bufferData(gl.ARRAY_BUFFER, c.subarray(0, 3 * d), gl.DYNAMIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER, colourBuf); gl.bufferData(gl.ARRAY_BUFFER, col.subarray(0, 4 * d), gl.DYNAMIC_DRAW);
  gl.bindBuffer(gl.ARRAY_BUFFER, pickBuf); gl.bufferData(gl.ARRAY_BUFFER, pk.subarray(0, 4 * d), gl.DYNAMIC_DRAW);
  const st = $("stats"); st.textContent = "";
  const b = (t) => { const e = document.createElement("b"); e.textContent = t; return e; };
  st.append(b(m.toLocaleString("en")), ` of ${N.toLocaleString("en")} voxel${N === 1 ? "" : "s"} shown · leaf area `, b((area * VOL).toFixed(2) + " m²"),
            ` of ${(areaAll * VOL).toFixed(2)} m² · voxel size ${H} m`);
  if (S.hover >= 0 && !shown(S.hover)) setHover(-1);
  if (S.table) buildTable();
  legend(); pickDirty = true; draw();
}

// ---------------------------------------------------------------- lines: box, axes, scanners, highlight
const EDGES = [[0, 1], [1, 3], [3, 2], [2, 0], [4, 5], [5, 7], [7, 6], [6, 4], [0, 4], [1, 5], [2, 6], [3, 7]];
function boxLines(c, h, out) {
  const p = []; for (let k = 0; k < 8; k++) p.push([c[0] + (k & 1 ? h[0] : -h[0]), c[1] + (k & 2 ? h[1] : -h[1]), c[2] + (k & 4 ? h[2] : -h[2])]);
  for (const [a, b] of EDGES) out.push(...p[a], ...p[b]);
}
const css = n => getComputedStyle(viz).getPropertyValue(n).trim();
const rgba = (h, a) => { const c = hex(h); return [c[0] / 255, c[1] / 255, c[2] / 255, a]; };
function drawLines(arr, colour) {
  if (!arr.length) return;
  gl.bindVertexArray(lineVao); gl.bindBuffer(gl.ARRAY_BUFFER, lineBuf);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(arr), gl.DYNAMIC_DRAW);
  gl.uniform4fv(U.lcol, colour); gl.drawArrays(gl.LINES, 0, arr.length / 3);
}
const labels = $("labels");
const AXN = ["x (m)", "y (m)", "z (m)"], axes = [];
function buildTicks() {
  const step = nice(Math.max(...half) * 2 / 5);
  labels.textContent = "";
  for (let a = 0; a < 3; a++) {
    const ax = {vals: [], els: []};
    for (let v = Math.ceil(wmin[a] / step - 1e-9) * step; v <= wmax[a] + 1e-9; v += step) {
      ax.vals.push(v); const e = document.createElement("span"); e.textContent = String(Math.round(v * 1000) / 1000); labels.append(e); ax.els.push(e);
    }
    ax.title = document.createElement("span"); ax.title.className = "t"; ax.title.textContent = AXN[a]; labels.append(ax.title);
    axes.push(ax);
  }
  for (const s of D.scans) { s.el = document.createElement("span"); s.el.className = "t"; s.el.textContent = s.name; labels.append(s.el); }
}
function put(el, q) { el.style.display = q ? "" : "none"; if (q) { el.style.left = q[0] + "px"; el.style.top = q[1] + "px"; } }
function drawAxes(lines) {
  // Label the box edges that face the viewer, so that the labels stay outside the voxels:
  // x and y along the near edges of the floor (or ceiling, seen from below), z along the leftmost vertical edge.
  const az = S.az * Math.PI / 180, el = S.el * Math.PI / 180;
  const dir = [Math.cos(el) * Math.cos(az), Math.cos(el) * Math.sin(az), Math.sin(el)];
  const above = S.el >= 0, steep = Math.abs(S.el) > 55;            // floor edges, or the ceiling when looking steeply down
  const sx = dir[0] >= 0 ? 1 : -1, sy = dir[1] >= 0 ? 1 : -1, zf = (above === steep ? 1 : -1) * half[2], tl = radius * 0.02;
  let zc = null, zx = Infinity;
  for (const cx of [-1, 1]) for (const cy of [-1, 1]) { const q = project([cx * half[0], cy * half[1], 0]); if (q && q[0] < zx) { zx = q[0]; zc = [cx, cy]; } }
  const edge = [{base: [0, sy * half[1], zf], out: [0, sy, 0]}, {base: [sx * half[0], 0, zf], out: [sx, 0, 0]},
                {base: zc ? [zc[0] * half[0], zc[1] * half[1], 0] : [0, 0, 0], out: zc ? [zc[0] * 0.7, zc[1] * 0.7, 0] : [0, 0, 0]}];
  const plen = (b, o) => { const p0 = project(b), p1 = project([0, 1, 2].map(i => b[i] + o[i] * tl)); return p0 && p1 ? Math.hypot(p1[0] - p0[0], p1[1] - p0[1]) : 0; };
  for (let a = 0; a < 3; a++) {
    const ax = axes[a], e = edge[a], hide = Math.abs(dir[a]) > 0.96 || (a === 2 && !zc);
    if (a < 2) { const down = [0, 0, zf < 0 ? -1 : 1]; if (plen(e.base, down) > 1.6 * plen(e.base, e.out)) e.out = down; }   // seen end-on: tick downwards instead
    const at = (t, k) => { const p = e.base.slice(); p[a] = t; return [0, 1, 2].map(i => p[i] + e.out[i] * tl * k); };
    const off = (t, px) => {                       // a point px pixels outward from the axis, on screen
      const p0 = project(at(t, 0)), p1 = project(at(t, 1)); if (!p0 || !p1) return null;
      const dx = p1[0] - p0[0], dy = p1[1] - p0[1], l = Math.hypot(dx, dy) || 1; return [p0[0] + dx / l * px, p0[1] + dy / l * px];
    };
    const e0 = project(at(-half[a], 0)), e1 = project(at(half[a], 0)), len = e0 && e1 ? Math.hypot(e1[0] - e0[0], e1[1] - e0[1]) : 0;
    const every = len < 70 ? 0 : Math.max(1, Math.ceil(26 / (len / Math.max(ax.vals.length - 1, 1))));
    ax.vals.forEach((v, i) => {
      const on = !hide && every > 0 && i % every === 0;
      if (on) lines.push(...at(v - centre[a], 0), ...at(v - centre[a], 1));
      put(ax.els[i], on ? off(v - centre[a], 15) : null);
    });
    const o0 = project(at(0, 0)), o1 = project(at(0, 1)), sideways = o0 && o1 && Math.abs(o1[0] - o0[0]) > Math.abs(o1[1] - o0[1]);
    put(ax.title, hide || len < 70 ? null : off(0, sideways ? 56 : 38));
  }
  for (const s of D.scans) put(s.el, S.scanners ? project([0, 1, 2].map(k => s.position[k] - centre[k] + (k === 2 ? H * 2.5 : 0))) : null);
}

// ---------------------------------------------------------------- drawing
function resize() {
  const r = Math.min(window.devicePixelRatio || 1, 2), w = Math.max(1, Math.round(stage.clientWidth * r)), h = Math.max(1, Math.round(stage.clientHeight * r));
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; pickDirty = true; }
}
function scene(pick) {
  gl.useProgram(cubeProg); gl.bindVertexArray(vao);
  gl.uniformMatrix4fv(U.mvp, false, MVP); gl.uniform1f(U.size, H * S.size); gl.uniform1i(U.pick, pick ? 1 : 0);
  if (nDraw) gl.drawArraysInstanced(gl.TRIANGLES, 0, 36, nDraw);
}
let queued = false;
function draw() { if (!queued) { queued = true; requestAnimationFrame(render); } }
function render() {
  queued = false; resize();
  if (autoFit) { const d = fitDist(); if (d !== S.dist) { S.dist = d; pickDirty = true; } }
  camera();
  gl.bindFramebuffer(gl.FRAMEBUFFER, null); gl.viewport(0, 0, canvas.width, canvas.height);
  const bg = rgba(css("--surface-1"), 1);
  gl.clearColor(bg[0], bg[1], bg[2], 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  gl.enable(gl.DEPTH_TEST); gl.enable(gl.CULL_FACE);
  scene(false);
  gl.useProgram(lineProg); gl.uniformMatrix4fv(U.lmvp, false, MVP);
  const box = []; boxLines([0, 0, 0], half, box);
  drawAxes(box);
  drawLines(box, rgba(css("--axis"), 1));
  if (S.scanners) {
    const sc = [], s = Math.max(H * 1.5, S.dist * 0.012);
    for (const q of D.scans) { const p = [0, 1, 2].map(k => q.position[k] - centre[k]);
      for (let a = 0; a < 3; a++) { const u = p.slice(), v = p.slice(); u[a] -= s; v[a] += s; sc.push(...u, ...v); } }
    drawLines(sc, rgba(css("--text-primary"), 1));
  }
  if (S.hover >= 0) {
    const c = [0, 1, 2].map(a => D.origin[a] + (A.ijk[3 * S.hover + a] + 0.5) * H - centre[a]), e = H * (S.size * 0.5 + 0.04), hl = [];
    boxLines(c, [e, e, e], hl); gl.disable(gl.DEPTH_TEST); drawLines(hl, rgba(css("--text-primary"), 1)); gl.enable(gl.DEPTH_TEST);
  }
  gl.bindVertexArray(null);
}
function pickAt(x, y) {          // CSS pixels in the stage -> voxel index or -1
  const r = canvas.width / stage.clientWidth;
  if (pickW !== canvas.width || pickH !== canvas.height) {
    pickW = canvas.width; pickH = canvas.height;
    gl.bindTexture(gl.TEXTURE_2D, pickTex); gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, pickW, pickH, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.bindRenderbuffer(gl.RENDERBUFFER, pickDepth); gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, pickW, pickH);
    gl.bindFramebuffer(gl.FRAMEBUFFER, pickFbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, pickTex, 0);
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, pickDepth);
    pickDirty = true;
  }
  gl.bindFramebuffer(gl.FRAMEBUFFER, pickFbo);
  if (pickDirty) {
    camera(); gl.viewport(0, 0, pickW, pickH); gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST); gl.enable(gl.CULL_FACE); scene(true); gl.bindVertexArray(null); pickDirty = false;
  }
  const px = new Uint8Array(4);
  gl.readPixels(Math.round(x * r), pickH - 1 - Math.round(y * r), 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  return (px[0] | (px[1] << 8) | (px[2] << 16)) - 1;
}

// ---------------------------------------------------------------- tooltip, legend, table
const tip = $("tooltip");
function row(tb, k, v) { const tr = tb.insertRow(); tr.insertCell().textContent = k; tr.insertCell().textContent = v; }
function setHover(i, x, y) {
  if (i !== S.hover) { S.hover = i; draw(); }
  if (i < 0) { tip.style.display = "none"; return; }
  tip.textContent = "";
  const v = document.createElement("div"); v.className = "v";
  v.textContent = Number.isFinite(A.lad[i]) ? A.lad[i].toFixed(3) + " m² m⁻³" : "no value";
  const tb = document.createElement("table");
  row(tb, "Leaf area", Number.isFinite(A.lad[i]) ? (A.lad[i] * VOL).toFixed(4) + " m²" : "–");
  row(tb, "Quality", FLAG[A.flag[i]] || String(A.flag[i]));
  const s = A.scan[i]; row(tb, "Scan used", s >= 1 && D.scans[s - 1] ? D.scans[s - 1].name : (s >= 1 ? "scan " + s : "none"));
  row(tb, "Explored", Number.isFinite(A.explored[i]) ? A.explored[i].toFixed(0) + " %" : "–");
  row(tb, "Pulses aimed / blocked", A.nt[i].toLocaleString("en") + " / " + A.nb[i].toLocaleString("en"));
  row(tb, "Returns leaf / wood / noise", A.nil[i] + " / " + A.niw[i] + " / " + A.nin[i]);
  row(tb, "Centre x, y, z (m)", [0, 1, 2].map(a => (D.origin[a] + (A.ijk[3 * i + a] + 0.5) * H).toFixed(2)).join(", "));
  row(tb, "Index i, j, k", [0, 1, 2].map(a => A.ijk[3 * i + a]).join(", "));
  tip.append(v, tb); tip.style.display = "block";
  if (x === undefined) { const q = project([0, 1, 2].map(a => D.origin[a] + (A.ijk[3 * i + a] + 0.5) * H - centre[a])) || [20, 20]; x = q[0]; y = q[1]; }
  const w = tip.offsetWidth, h = tip.offsetHeight;
  tip.style.left = Math.max(4, Math.min(stage.clientWidth - w - 4, x + 16)) + "px";
  tip.style.top = Math.max(4, Math.min(stage.clientHeight - h - 4, y + 16)) + "px";
}
function legend() {
  const L = $("legend"); L.textContent = "";
  const ttl = t => { const e = document.createElement("div"); e.className = "ttl"; e.textContent = t; L.append(e); };
  const sw = (colour, label, n, btn, on) => {
    const e = document.createElement(btn ? "button" : "div"); e.className = "row";
    const b = document.createElement("span"); b.className = "sw"; b.style.background = colour;
    const t = document.createElement("span"); t.textContent = label;
    e.append(b, t);
    if (n !== undefined) { const c = document.createElement("span"); c.className = "n"; c.textContent = n.toLocaleString("en"); e.append(c); }
    if (btn) { e.type = "button"; e.setAttribute("aria-pressed", String(on)); e.title = "Show or hide"; e.onclick = btn; }
    L.append(e);
  };
  if (S.mode === "lad" || S.mode === "explored") {
    ttl(S.mode === "lad" ? "Leaf area density (m² m⁻³)" : "Pulses that reached the voxel (%)");
    const bar = document.createElement("div"); bar.className = "bar";
    const stops = []; for (let k = 0; k <= 12; k++) { const c = ramp(k / 12); stops.push(`rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`); }
    bar.style.background = `linear-gradient(to right, ${stops.join(",")})`;
    const ends = document.createElement("div"); ends.className = "ends";
    const a = document.createElement("span"), b = document.createElement("span");
    a.textContent = "0"; b.textContent = S.mode === "lad" ? (ladMaxData > S.cmax ? "≥ " : "") + S.cmax : "100";
    ends.append(a, b); L.append(bar, ends);
    let nn = 0; for (let k = 0; k < nShown; k++) if (!Number.isFinite(S.mode === "lad" ? A.lad[vis[k]] : A.explored[vis[k]])) nn++;
    if (nn) sw(NODATA, "no value", nn);
  } else if (S.mode === "flag") {
    ttl("Quality flag");
    const n = [0, 0, 0, 0, 0]; for (let k = 0; k < nShown; k++) n[A.flag[vis[k]]]++;
    for (let f = 0; f < 5; f++) if (flagPresent[f]) sw(FLAG_COL[f], FLAG_ICON[f] + " " + FLAG[f], n[f]);
  } else {
    ttl("Scan used");
    const n = {}; for (let k = 0; k < nShown; k++) n[A.scan[vis[k]]] = (n[A.scan[vis[k]]] || 0) + 1;
    for (const s of scansPresent) if (s <= 8) sw(scanColour(s), s >= 1 && D.scans[s - 1] ? `${s} ${D.scans[s - 1].name}` : (s ? "scan " + s : "none"), n[s] || 0,
      () => { S.scans[s] = S.scans[s] === false; rebuild(); }, S.scans[s] !== false);
    const more = scansPresent.filter(s => s > 8);
    if (more.length) {
      const on = more.some(s => S.scans[s] !== false);
      sw(scanColour(9), more.length > 1 ? `scans ${more[0]} to ${more[more.length - 1]}` : "scan " + more[0], more.reduce((t, s) => t + (n[s] || 0), 0),
        () => { for (const s of more) S.scans[s] = !on; rebuild(); }, on);
    }
  }
}
function buildTable() {
  const W = $("tablewrap"); W.textContent = "";
  const order = Array.from(vis).sort((a, b) => (Number.isFinite(A.lad[b]) ? A.lad[b] : -1) - (Number.isFinite(A.lad[a]) ? A.lad[a] : -1)).slice(0, 500);
  const t = document.createElement("table"), cap = t.createCaption();
  cap.textContent = `${order.length.toLocaleString("en")} of the ${nShown.toLocaleString("en")} voxels shown, highest leaf area density first. ` +
    (D.csv ? `All voxels are in ${D.csv}.` : "All voxels are in the voxel table (CSV file) of the run.");
  const hd = t.createTHead().insertRow();
  for (const [h, l] of [["i", 0], ["j", 0], ["k", 0], ["x (m)", 0], ["y (m)", 0], ["z (m)", 0], ["LAD (m² m⁻³)", 0], ["Leaf area (m²)", 0], ["Quality", 1], ["Scan", 1], ["Explored (%)", 0], ["Aimed", 0], ["Blocked", 0], ["Leaf returns", 0]]) {
    const th = document.createElement("th"); th.textContent = h; if (l) th.className = "l"; hd.append(th); }
  const tb = t.createTBody();
  for (const i of order) {
    const tr = tb.insertRow(), c = (v, l) => { const td = tr.insertCell(); td.textContent = v; if (l) td.className = "l"; };
    for (let a = 0; a < 3; a++) c(A.ijk[3 * i + a]);
    for (let a = 0; a < 3; a++) c((D.origin[a] + (A.ijk[3 * i + a] + 0.5) * H).toFixed(2));
    c(fmt(A.lad[i], 3)); c(Number.isFinite(A.lad[i]) ? (A.lad[i] * VOL).toFixed(4) : "–"); c(FLAG[A.flag[i]], 1);
    c(A.scan[i] >= 1 && D.scans[A.scan[i] - 1] ? D.scans[A.scan[i] - 1].name : "none", 1);
    c(Number.isFinite(A.explored[i]) ? A.explored[i].toFixed(0) : "–"); c(A.nt[i]); c(A.nb[i]); c(A.nil[i]);
  }
  W.append(t);
}

// ---------------------------------------------------------------- controls
$("title").textContent = D.title;
const flagPresent = [false, false, false, false, false], scanSet = new Set();
for (let i = 0; i < N; i++) { flagPresent[A.flag[i]] = true; scanSet.add(A.scan[i]); }
const scansPresent = Array.from(scanSet).sort((a, b) => a - b);
const ft = $("flagtoggles");
for (let f = 0; f < 5; f++) if (flagPresent[f]) {
  const l = document.createElement("label"), c = document.createElement("input");
  c.type = "checkbox"; c.checked = true; c.onchange = () => { S.flags[f] = c.checked; rebuild(); };
  l.append(c, " " + FLAG[f]); ft.append(l);
}
const sl = (id, min, max, step, val, fn) => { const e = $(id); e.min = min; e.max = max; e.step = step; e.value = val; e.oninput = () => fn(+e.value, e); return e; };
const ladStep = ladMaxData > 0 ? ladMaxData / 200 : 0.01;
const eMin = sl("ladmin", 0, ladMaxData || 1, ladStep, 0, v => { S.ladMin = Math.min(v < ladStep / 2 ? 0 : v, S.ladMax); syncLad(); rebuild(); });
const eMax = sl("ladmax", 0, ladMaxData || 1, ladStep, ladMaxData || 1, v => { S.ladMax = Math.max(v > ladMaxData - ladStep / 2 ? ladMaxData : v, S.ladMin); syncLad(); rebuild(); });
function syncLad() { eMin.value = S.ladMin; eMax.value = S.ladMax; $("ladminv").textContent = S.ladMin.toFixed(2); $("ladmaxv").textContent = S.ladMax.toFixed(2); }
const eLo = $("cliplo"), eHi = $("cliphi");
function syncClip() {
  const a = S.clipAxis; eLo.min = eHi.min = lo[a]; eLo.max = eHi.max = hi[a]; eLo.step = eHi.step = 1; eLo.value = S.clip[a][0]; eHi.value = S.clip[a][1];
  $("clipv").textContent = (D.origin[a] + S.clip[a][0] * H).toFixed(1) + " to " + (D.origin[a] + (S.clip[a][1] + 1) * H).toFixed(1) + " m";
  // cuts on the three axes add up: say which axes are cut, and offer to undo them
  const isCut = b => S.clip[b][0] > lo[b] || S.clip[b][1] < hi[b];
  for (const o of $("clipaxis").options) o.textContent = ["x", "y", "height (z)"][+o.value] + (isCut(+o.value) ? " (cut)" : "");
  $("clipreset").hidden = !(isCut(0) || isCut(1) || isCut(2));
}
$("clipreset").onclick = () => { for (let b = 0; b < 3; b++) S.clip[b] = [lo[b], hi[b]]; syncClip(); rebuild(); };
eLo.oninput = () => { const a = S.clipAxis; S.clip[a][0] = Math.min(+eLo.value, S.clip[a][1]); syncClip(); rebuild(); };
eHi.oninput = () => { const a = S.clipAxis; S.clip[a][1] = Math.max(+eHi.value, S.clip[a][0]); syncClip(); rebuild(); };
$("clipaxis").onchange = e => { S.clipAxis = +e.target.value; syncClip(); };
const cm = $("cmax"); cm.value = S.cmax; cm.step = S.cmax / 10;
cm.onchange = () => { const v = +cm.value; if (v > 0) { S.cmax = v; rebuild(); } else cm.value = S.cmax; };
$("mode").onchange = e => { S.mode = e.target.value; $("cmaxlab").style.display = S.mode === "lad" ? "" : "none"; rebuild(); };
$("hidezero").onchange = e => { S.hideZero = e.target.checked; rebuild(); };
$("size").oninput = e => {
  S.size = +e.target.value / 100;
  if (nShown > CULL_FROM && (S.size >= 0.85) !== culled) rebuild();     // small voxels let the inside show: draw all of them
  else { pickDirty = true; draw(); }
};
if (!D.scans.length) $("scanners").parentElement.style.display = "none";      // a page made from a voxel table has no scanner positions
$("scanners").onchange = e => {
  S.scanners = e.target.checked;
  S.target = [0, 0, 0]; autoFit = true;           // zoom out far enough to see them, or back in on the voxels
  pickDirty = true; draw();
};
$("ortho").onchange = e => { S.ortho = e.target.checked; pickDirty = true; draw(); };
const VIEWS = {iso: [-58, 22], top: [-90, 90], front: [-90, 0], side: [0, 0]};
function view(name) {
  [S.az, S.el] = VIEWS[name]; S.target = [0, 0, 0]; autoFit = true;
  S.ortho = name !== "iso"; $("ortho").checked = S.ortho; pickDirty = true; draw();
}
for (const b of document.querySelectorAll("[data-view]")) b.onclick = () => view(b.dataset.view);
$("tablebtn").onclick = e => { S.table = !S.table; e.target.setAttribute("aria-pressed", String(S.table)); $("tablewrap").hidden = !S.table; if (S.table) buildTable(); draw(); };
$("theme").onclick = e => { const d = !dark(); root.setAttribute("data-theme", d ? "dark" : "light"); e.target.textContent = d ? "Light" : "Dark"; rebuild(); };
$("theme").textContent = dark() ? "Light" : "Dark";
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { $("theme").textContent = dark() ? "Light" : "Dark"; rebuild(); });
$("png").onclick = () => {                     // the picture: the voxels as drawn, with the axis labels and the legend
  render();
  const k = canvas.width / Math.max(stage.clientWidth, 1), out = document.createElement("canvas");
  out.width = canvas.width; out.height = canvas.height;
  const g = out.getContext("2d"), sb = stage.getBoundingClientRect();
  const box = e => { const q = e.getBoundingClientRect(); return [q.left - sb.left, q.top - sb.top, q.width, q.height]; };
  const font = (px, w) => `${w || 400} ${px}px system-ui, -apple-system, "Segoe UI", sans-serif`;
  g.drawImage(canvas, 0, 0); g.scale(k, k); g.textBaseline = "middle";
  g.textAlign = "center"; g.font = font(11);
  for (const e of labels.children) if (e.style.display !== "none") {
    g.fillStyle = css(e.classList.contains("t") ? "--text-secondary" : "--muted"); g.fillText(e.textContent, parseFloat(e.style.left), parseFloat(e.style.top)); }
  const L = $("legend"), [lx, ly, lw, lh] = box(L);
  g.fillStyle = css("--surface-1"); g.strokeStyle = css("--border"); g.lineWidth = 1;
  g.beginPath(); g.rect(lx + 0.5, ly + 0.5, lw - 1, lh - 1); g.fill(); g.stroke();
  for (const e of L.querySelectorAll(".ttl, .bar, .ends span, .row > span")) {
    const [x, y, w, h] = box(e), row = e.closest(".row");
    g.globalAlpha = row && row.getAttribute("aria-pressed") === "false" ? 0.45 : 1;
    if (e.classList.contains("bar")) { for (let i = 0; i < w; i++) { const c = ramp(i / Math.max(w - 1, 1)); g.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`; g.fillRect(x + i, y, 1.5, h); } continue; }
    if (e.classList.contains("sw")) { g.fillStyle = e.style.background; g.fillRect(x, y, w, h); continue; }
    const right = e.classList.contains("n") || (e.parentElement.classList.contains("ends") && e === e.parentElement.lastElementChild);
    g.font = font(13, e.classList.contains("ttl") ? 600 : 400); g.textAlign = right ? "right" : "left";
    g.fillStyle = css(e.classList.contains("ttl") ? "--text-primary" : "--text-secondary");
    g.fillText(e.textContent, right ? x + w : x, y + h / 2);
  }
  out.toBlob(b => { const a = document.createElement("a"); a.href = URL.createObjectURL(b); a.download = (D.file || "voxlad") + "_view.png"; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); });
};

// ---------------------------------------------------------------- pointer, wheel, keyboard
const ptr = new Map();
let pinch = 0, moved = false, hoverAt = null, hoverQueued = false;
function doHover() { hoverQueued = false; if (hoverAt && !ptr.size) setHover(pickAt(hoverAt[0], hoverAt[1]), hoverAt[0], hoverAt[1]); }
canvas.addEventListener("contextmenu", e => e.preventDefault());
canvas.addEventListener("pointerdown", e => { canvas.setPointerCapture(e.pointerId); ptr.set(e.pointerId, [e.clientX, e.clientY, e.button, e.shiftKey]); canvas.classList.add("dragging"); moved = false;
  if (ptr.size === 2) { const p = Array.from(ptr.values()); pinch = Math.hypot(p[0][0] - p[1][0], p[0][1] - p[1][1]); } });
const up = e => { ptr.delete(e.pointerId); if (!ptr.size) canvas.classList.remove("dragging"); };
canvas.addEventListener("pointerup", up); canvas.addEventListener("pointercancel", up);
function pan(dx, dy) { const k = 2 * S.dist * Math.tan(FOV / 2) / stage.clientHeight; for (let a = 0; a < 3; a++) S.target[a] += (-dx * camRight[a] + dy * camUp[a]) * k; }
canvas.addEventListener("pointermove", e => {
  const b = stage.getBoundingClientRect(), p = ptr.get(e.pointerId);
  if (!p) { if (!ptr.size) { hoverAt = [e.clientX - b.left, e.clientY - b.top]; if (!hoverQueued) { hoverQueued = true; requestAnimationFrame(doHover); } } return; }
  const dx = e.clientX - p[0], dy = e.clientY - p[1]; p[0] = e.clientX; p[1] = e.clientY;
  if (dx || dy) moved = true;
  if (ptr.size === 2) {
    const q = Array.from(ptr.values()), d = Math.hypot(q[0][0] - q[1][0], q[0][1] - q[1][1]);
    if (pinch > 0 && d > 0) S.dist = Math.min(radius * 80, Math.max(radius * 0.05, S.dist * pinch / d));
    pinch = d; pan(dx / 2, dy / 2); autoFit = false;
  } else if (p[2] === 2 || p[2] === 1 || p[3]) { pan(dx, dy); autoFit = false; }
  else { autoFit = false; S.az -= dx * 0.35; S.el = Math.max(-90, Math.min(90, S.el + dy * 0.35)); }
  setHover(-1); pickDirty = true; draw();
});
canvas.addEventListener("pointerleave", () => { hoverAt = null; if (!ptr.size) setHover(-1); });
canvas.addEventListener("wheel", e => { e.preventDefault(); autoFit = false; S.dist = Math.min(radius * 80, Math.max(radius * 0.05, S.dist * Math.exp(e.deltaY * 0.0012))); pickDirty = true; setHover(-1); draw(); }, {passive: false});
canvas.addEventListener("dblclick", () => view("iso"));
canvas.addEventListener("keydown", e => {
  const k = e.key; let used = true; const wasFit = autoFit; autoFit = false;
  if (k === "ArrowLeft") S.az += 5; else if (k === "ArrowRight") S.az -= 5;
  else if (k === "ArrowUp") S.el = Math.min(90, S.el + 5); else if (k === "ArrowDown") S.el = Math.max(-90, S.el - 5);
  else if (k === "+" || k === "=") S.dist *= 0.9; else if (k === "-") S.dist /= 0.9;
  else if (k === "0" || k === "Home") { view("iso"); }
  else { used = false; autoFit = wasFit; }
  if (used) { e.preventDefault(); pickDirty = true; draw(); }
});
new ResizeObserver(() => { pickDirty = true; draw(); }).observe(stage);

syncLad(); syncClip(); buildTicks(); rebuild();
window.voxladViewer = {state: S, data: D, arrays: A, centre, shown: () => nShown, drawn: () => nDraw, rebuild, pickAt, hover: setHover, project};     // for scripted checks
})();
</script>
</body>
</html>
"""
