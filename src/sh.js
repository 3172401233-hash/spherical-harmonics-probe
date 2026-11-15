// Order-2 (9 coefficient) real spherical harmonics lighting core.
//
// Conventions (explicit, do not rely on pipeline defaults):
//
// World axes: right-handed, +X right, +Y up, +Z forward.
// All directions are unit vectors [x, y, z] in this world space.
//
// Coefficient order: ACN (Ambisonic Channel Number), index = l*(l+1)+m:
//   0: Y00   1: Y1-1  2: Y10  3: Y11
//   4: Y2-2  5: Y2-1  6: Y20  7: Y21  8: Y22
//
// Normalization: N3D (fully normalized, orthonormal basis) is the canonical
// internal convention. SN3D (as used by some engines' ambient slots) is
// available via toSN3D/fromSN3D. Per-band scale between them:
//   coeff_SN3D(l) = coeff_N3D(l) * sqrt(4*PI / (2l+1))

export const SH_ORDER = 2;
export const SH_COEFFICIENT_COUNT = 9;

export const SQRT_4PI = 3.5449077018110318;

// N3D basis constants.
const K0 = 0.28209479177387814; // sqrt(1/(4*pi))
const K1 = 0.4886025119029199; // sqrt(3/(4*pi))
const K2 = 1.0925484305920792; // sqrt(15/(4*pi))
const K3 = 0.31539156525352005; // sqrt(5/(16*pi))
const K4 = 0.5462742152960396; // sqrt(15/(16*pi))

// Per-band conversion factor N3D -> SN3D: sqrt(4*pi / (2l+1)).
export const SN3D_BAND_SCALE = [
  SQRT_4PI, // l = 0
  Math.sqrt((4 * Math.PI) / 3), // l = 1
  Math.sqrt((4 * Math.PI) / 5), // l = 2
];

// Band index for each ACN coefficient.
export const ACN_BAND = [0, 1, 1, 1, 2, 2, 2, 2, 2];

// Lambertian convolution coefficients A_l for irradiance reconstruction.
// E(n) = sum_lm A_l * c_lm * Y_lm(n)
export const LAMBERTIAN_BAND_A = [Math.PI, (2 * Math.PI) / 3, Math.PI / 4];

// Evaluate the 9 N3D basis functions at unit direction [x, y, z].
export function basis(direction, out = new Float64Array(SH_COEFFICIENT_COUNT)) {
  const [x, y, z] = direction;
  out[0] = K0;
  out[1] = K1 * y;
  out[2] = K1 * z;
  out[3] = K1 * x;
  out[4] = K2 * x * y;
  out[5] = K2 * y * z;
  out[6] = K3 * (3 * z * z - 1);
  out[7] = K2 * x * z;
  out[8] = K4 * (x * x - y * y);
  return out;
}

export function emptyRGBCoefficients() {
  return Array.from({ length: SH_COEFFICIENT_COUNT }, () => [0, 0, 0]);
}

// Accumulate one directional sample into radiance coefficients.
// radiance: [r, g, b], weight: solid angle of the sample (steradians).
export function accumulateSample(coefficients, direction, radiance, weight) {
  const b = basis(direction);
  for (let i = 0; i < SH_COEFFICIENT_COUNT; i++) {
    const w = b[i] * weight;
    coefficients[i][0] += radiance[0] * w;
    coefficients[i][1] += radiance[1] * w;
    coefficients[i][2] += radiance[2] * w;
  }
}

// Project an environment radiance function into SH radiance coefficients.
// radianceFn(direction) -> [r, g, b].
// sampler: iterable of { direction, weight } covering the sphere.
export function projectEnvironment(radianceFn, sampler) {
  const coefficients = emptyRGBCoefficients();
  for (const sample of sampler) {
    accumulateSample(coefficients, sample.direction, radianceFn(sample.direction), sample.weight);
  }
  return coefficients;
}

// Exact projection of a delta light: radiance = flux * delta(direction - lightDir).
export function projectDeltaLight(lightDir, flux) {
  const b = basis(lightDir);
  const coefficients = emptyRGBCoefficients();
  for (let i = 0; i < SH_COEFFICIENT_COUNT; i++) {
    coefficients[i] = [flux[0] * b[i], flux[1] * b[i], flux[2] * b[i]];
  }
  return coefficients;
}

// Convert radiance coefficients to Lambertian irradiance coefficients
// by multiplying each band by A_l. Evaluate with evaluate() to get
// irradiance E(n) in the same units the radiance was expressed in.
export function irradianceFromRadiance(coefficients) {
  return coefficients.map((c, i) => {
    const a = LAMBERTIAN_BAND_A[ACN_BAND[i]];
    return [c[0] * a, c[1] * a, c[2] * a];
  });
}

// Reconstruct the spherical signal (radiance or irradiance, depending on
// what the coefficients represent) at a unit direction.
export function evaluate(coefficients, direction) {
  const b = basis(direction);
  let r = 0;
  let g = 0;
  let bl = 0;
  for (let i = 0; i < SH_COEFFICIENT_COUNT; i++) {
    r += coefficients[i][0] * b[i];
    g += coefficients[i][1] * b[i];
    bl += coefficients[i][2] * b[i];
  }
  return [r, g, bl];
}

// Total signal energy: sum_i |c_i|^2 == integral of f^2 over the sphere
// (Parseval, exact because the N3D basis is orthonormal).
export function energy(coefficients) {
  let total = 0;
  for (const c of coefficients) {
    total += c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
  }
  return total;
}

// N3D -> SN3D and back (per-band uniform scale).
export function toSN3D(coefficients) {
  return coefficients.map((c, i) => {
    const s = SN3D_BAND_SCALE[ACN_BAND[i]];
    return [c[0] * s, c[1] * s, c[2] * s];
  });
}

export function fromSN3D(coefficients) {
  return coefficients.map((c, i) => {
    const s = 1 / SN3D_BAND_SCALE[ACN_BAND[i]];
    return [c[0] * s, c[1] * s, c[2] * s];
  });
}
