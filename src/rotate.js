// SH coefficient rotation for order-2 probes.
//
// rotateCoefficients applies an active rotation R (3x3, row-major, world
// space) to the environment: the rotated signal satisfies
//   f'(d) = f(R^-1 * d)
// Use it when the skybox or the whole scene is rotated by R, or to move a
// probe between a rotated frame and world space (use the inverse rotation).
//
// The 9x9 rotation matrix is built by numerically projecting rotated basis
// functions back onto the basis. Band-limited rotation stays within the
// band, so with a dense Fibonacci quadrature this is accurate to ~1e-6.

import { basis, SH_COEFFICIENT_COUNT } from './sh.js';
import { fibonacciSphere } from './sampling.js';

function transpose3(m) {
  return [
    [m[0][0], m[1][0], m[2][0]],
    [m[0][1], m[1][1], m[2][1]],
    [m[0][2], m[1][2], m[2][2]],
  ];
}

function applyMat3(m, v) {
  return [
    m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
    m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
    m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
  ];
}

// Build the 9x9 SH rotation matrix M such that c' = M * c rotates the
// signal by R (active rotation, orthonormal 3x3 with det +1).
export function shRotationMatrix(rotation, quadratureSamples = 4096) {
  const inv = transpose3(rotation); // rotation must be orthonormal
  const M = Array.from({ length: SH_COEFFICIENT_COUNT }, () =>
    new Float64Array(SH_COEFFICIENT_COUNT)
  );
  const bAt = new Float64Array(SH_COEFFICIENT_COUNT);
  const bRot = new Float64Array(SH_COEFFICIENT_COUNT);
  for (const { direction, weight } of fibonacciSphere(quadratureSamples)) {
    basis(direction, bAt);
    basis(applyMat3(inv, direction), bRot);
    for (let i = 0; i < SH_COEFFICIENT_COUNT; i++) {
      const w = bAt[i] * weight;
      for (let j = 0; j < SH_COEFFICIENT_COUNT; j++) {
        M[i][j] += w * bRot[j];
      }
    }
  }
  return M;
}

// Rotate RGB SH coefficients by a world-space rotation matrix (3x3).
export function rotateCoefficients(coefficients, rotation, quadratureSamples = 4096) {
  const M = shRotationMatrix(rotation, quadratureSamples);
  const out = Array.from({ length: SH_COEFFICIENT_COUNT }, () => [0, 0, 0]);
  for (let i = 0; i < SH_COEFFICIENT_COUNT; i++) {
    for (let j = 0; j < SH_COEFFICIENT_COUNT; j++) {
      const m = M[i][j];
      out[i][0] += m * coefficients[j][0];
      out[i][1] += m * coefficients[j][1];
      out[i][2] += m * coefficients[j][2];
    }
  }
  return out;
}

// Compose a rotation matrix from axis + angle (right-hand rule, radians).
export function rotationFromAxisAngle(axis, angle) {
  const [x, y, z] = axis;
  const len = Math.hypot(x, y, z);
  const nx = x / len;
  const ny = y / len;
  const nz = z / len;
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  const t = 1 - c;
  return [
    [t * nx * nx + c, t * nx * ny - s * nz, t * nx * nz + s * ny],
    [t * nx * ny + s * nz, t * ny * ny + c, t * ny * nz - s * nx],
    [t * nx * nz - s * ny, t * ny * nz + s * nx, t * nz * nz + c],
  ];
}
