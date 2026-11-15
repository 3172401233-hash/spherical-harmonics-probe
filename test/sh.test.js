import test from 'node:test';
import assert from 'node:assert/strict';
import {
  basis,
  projectEnvironment,
  projectDeltaLight,
  irradianceFromRadiance,
  evaluate,
  energy,
  toSN3D,
  fromSN3D,
  fibonacciSphere,
  cubemapTexels,
  rotateCoefficients,
  rotationFromAxisAngle,
  SQRT_4PI,
} from '../src/index.js';

const TOL = 1e-6;

function closeTo(actual, expected, eps, message) {
  assert.ok(
    Math.abs(actual - expected) <= eps,
    `${message ?? 'value'}: expected ${expected}, got ${actual}`
  );
}

function dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function normalize(v) {
  const len = Math.hypot(v[0], v[1], v[2]);
  return [v[0] / len, v[1] / len, v[2] / len];
}

function applyMat3(m, v) {
  return [
    m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
    m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
    m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
  ];
}

test('ACN order and N3D basis values at axis directions', () => {
  // Index map: 0:Y00 1:Y1-1(y) 2:Y10(z) 3:Y11(x) 4:Y2-2(xy) 5:Y2-1(yz)
  //            6:Y20(3z^2-1) 7:Y21(xz) 8:Y22(x^2-y^2)
  const atZ = basis([0, 0, 1]);
  closeTo(atZ[0], 0.28209479177387814, TOL, 'Y00');
  closeTo(atZ[2], 0.4886025119029199, TOL, 'Y10 at +Z');
  closeTo(atZ[1], 0, TOL, 'Y1-1 at +Z');
  closeTo(atZ[3], 0, TOL, 'Y11 at +Z');
  closeTo(atZ[6], 0.63078313050504, TOL, 'Y20 at +Z');
  closeTo(atZ[8], 0, TOL, 'Y22 at +Z');

  const atX = basis([1, 0, 0]);
  closeTo(atX[3], 0.4886025119029199, TOL, 'Y11 at +X');
  closeTo(atX[6], -0.31539156525352005, TOL, 'Y20 at +X');
  closeTo(atX[8], 0.5462742152960396, TOL, 'Y22 at +X');

  const atY = basis([0, 1, 0]);
  closeTo(atY[1], 0.4886025119029199, TOL, 'Y1-1 at +Y');
  closeTo(atY[8], -0.5462742152960396, TOL, 'Y22 at +Y');

  // Orthonormality: integral of Y_i * Y_j over the sphere = delta_ij.
  const samples = [...fibonacciSphere(4096)];
  for (let i = 0; i < 9; i++) {
    for (let j = 0; j < 9; j++) {
      let sum = 0;
      for (const { direction, weight } of samples) {
        const b = basis(direction);
        sum += b[i] * b[j] * weight;
      }
      closeTo(sum, i === j ? 1 : 0, 1e-4, `orthonormality (${i},${j})`);
    }
  }
});

test('constant sky: coefficients, reconstruction brightness and energy', () => {
  const L = [2.0, 1.0, 0.5];
  const sky = () => L;
  const radiance = projectEnvironment(sky, fibonacciSphere(8192));

  // Only Y00 survives: c0 = L * sqrt(4*pi), everything else ~ 0.
  for (let ch = 0; ch < 3; ch++) {
    closeTo(radiance[0][ch], L[ch] * SQRT_4PI, 1e-6, `c0 channel ${ch}`);
    for (let i = 1; i < 9; i++) {
      closeTo(radiance[i][ch], 0, 1e-4, `c${i} channel ${ch}`);
    }
  }

  // Lambertian irradiance of a constant sky is pi * L, independent of normal.
  const irradiance = irradianceFromRadiance(radiance);
  const normals = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    normalize([1, 2, 3]),
    normalize([-3, 1, -2]),
  ];
  for (const n of normals) {
    const e = evaluate(irradiance, n);
    for (let ch = 0; ch < 3; ch++) {
      closeTo(e[ch], Math.PI * L[ch], 1e-4, `irradiance at ${n} channel ${ch}`);
    }
  }

  // Energy: integral of |f|^2 = 4*pi*|L|^2 (Parseval).
  const expectedEnergy = 4 * Math.PI * (L[0] ** 2 + L[1] ** 2 + L[2] ** 2);
  closeTo(energy(radiance), expectedEnergy, 1e-4, 'constant sky energy');
});

test('cubemap sampler: weights sum to 4*pi and match analytic projection', () => {
  let total = 0;
  for (const { weight } of cubemapTexels(64)) total += weight;
  closeTo(total, 4 * Math.PI, 1e-9, 'cubemap solid angle sum');

  const L = [1.5, 1.0, 0.25];
  const radiance = projectEnvironment(() => L, cubemapTexels(64));
  for (let ch = 0; ch < 3; ch++) {
    closeTo(radiance[0][ch], L[ch] * SQRT_4PI, 1e-4, `cubemap c0 channel ${ch}`);
    for (let i = 1; i < 9; i++) {
      closeTo(radiance[i][ch], 0, 1e-3, `cubemap c${i} channel ${ch}`);
    }
  }
});

test('single directional light: brightness peak, direction and energy', () => {
  const lightDir = normalize([1, 0.5, 0.25]);
  const flux = [3, 3, 3];
  const radiance = projectDeltaLight(lightDir, flux);
  const irradiance = irradianceFromRadiance(radiance);

  // Addition theorem: E(n) = E0 * sum_l A_l (2l+1)/(4pi) P_l(n . d).
  // At n = d:  E0 * (1/4 + 1/2 + 5/16) = 1.0625 * E0.
  // At n = -d: E0 * (1/4 - 1/2 + 5/16) = 0.0625 * E0.
  const eFront = evaluate(irradiance, lightDir)[0];
  const eBack = evaluate(irradiance, [-lightDir[0], -lightDir[1], -lightDir[2]])[0];
  closeTo(eFront, 1.0625 * flux[0], 1e-9, 'irradiance toward light');
  closeTo(eBack, 0.0625 * flux[0], 1e-9, 'irradiance away from light');
  assert.ok(eFront > eBack, 'lit side brighter than dark side');

  // The brightest reconstructed direction must align with the light.
  let bestDir = null;
  let bestValue = -Infinity;
  for (const { direction } of fibonacciSphere(4096)) {
    const v = evaluate(irradiance, direction)[0];
    if (v > bestValue) {
      bestValue = v;
      bestDir = direction;
    }
  }
  assert.ok(
    dot(bestDir, lightDir) > 0.999,
    `peak direction ${bestDir} should align with light ${lightDir}`
  );

  // Linear term of the irradiance points at the light:
  // gradient direction ~ (c3, c1, c2) for (x, y, z).
  const linear = normalize([irradiance[3][0], irradiance[1][0], irradiance[2][0]]);
  assert.ok(dot(linear, lightDir) > 0.999999, 'linear term aligns with light');

  // Energy: sum_i c_i^2 = E0^2 * sum_l (2l+1)/(4pi) = E0^2 * 9/(4pi).
  closeTo(energy(radiance), 3 * flux[0] ** 2 * (9 / (4 * Math.PI)), 1e-9, 'delta light energy');
});

test('rotated asymmetric sky: rotated coefficients match re-projection', () => {
  // Asymmetric sky with distinct content on every axis.
  const sky = (d) => {
    const v = 1 + 2 * Math.max(0, d[0]) + 0.75 * Math.max(0, -d[1]) + 0.5 * d[2] ** 2;
    return [v, 0.5 * v + 0.25 * Math.max(0, d[1]), 0.25 * v];
  };
  const rotation = rotationFromAxisAngle(normalize([0.3, 1, -0.2]), Math.PI / 3);
  const inverse = rotationFromAxisAngle(normalize([0.3, 1, -0.2]), -Math.PI / 3);

  const sampler = () => fibonacciSphere(8192);
  const original = projectEnvironment(sky, sampler());

  // Path A: rotate the coefficients.
  const rotated = rotateCoefficients(original, rotation);

  // Path B: rotate the sky itself, then project. f'(d) = f(R^-1 d).
  const rotatedSky = (d) => sky(applyMat3(inverse, d));
  const expected = projectEnvironment(rotatedSky, sampler());

  for (let i = 0; i < 9; i++) {
    for (let ch = 0; ch < 3; ch++) {
      closeTo(rotated[i][ch], expected[i][ch], 1e-3, `rotated c${i} channel ${ch}`);
    }
  }

  // Energy is invariant under rotation.
  closeTo(energy(rotated), energy(original), 1e-3, 'energy preserved by rotation');

  // Reconstruction consistency: E'(R d) == E(d) for irradiance probes.
  const irrOriginal = irradianceFromRadiance(original);
  const irrRotated = irradianceFromRadiance(rotated);
  for (const { direction } of fibonacciSphere(256)) {
    const a = evaluate(irrRotated, applyMat3(rotation, direction));
    const b = evaluate(irrOriginal, direction);
    for (let ch = 0; ch < 3; ch++) {
      closeTo(a[ch], b[ch], 1e-3, `irradiance consistency at ${direction} ch ${ch}`);
    }
  }
});

test('SN3D conversion round-trips and scales per band', () => {
  const radiance = projectEnvironment((d) => [1 + d[0], 1 + d[1], 1 + d[2]], fibonacciSphere(2048));
  const sn3d = toSN3D(radiance);
  const restored = fromSN3D(sn3d);
  for (let i = 0; i < 9; i++) {
    for (let ch = 0; ch < 3; ch++) {
      closeTo(restored[i][ch], radiance[i][ch], 1e-12, `round trip c${i} ch${ch}`);
    }
  }
  // Band scales: sqrt(4pi/(2l+1)).
  closeTo(sn3d[0][0] / radiance[0][0], Math.sqrt(4 * Math.PI), 1e-12, 'l0 scale');
  closeTo(sn3d[3][0] / radiance[3][0], Math.sqrt((4 * Math.PI) / 3), 1e-12, 'l1 scale');
});
