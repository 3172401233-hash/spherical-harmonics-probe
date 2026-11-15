// Direction samplers with solid-angle weights for SH projection.

// Near-uniform sphere sampling via the Fibonacci lattice.
// Every sample carries weight 4*pi/count.
export function* fibonacciSphere(count) {
  const goldenAngle = Math.PI * (3 - Math.sqrt(5));
  const weight = (4 * Math.PI) / count;
  for (let i = 0; i < count; i++) {
    const y = 1 - (2 * (i + 0.5)) / count;
    const radius = Math.sqrt(Math.max(0, 1 - y * y));
    const phi = goldenAngle * i;
    yield {
      direction: [radius * Math.cos(phi), y, radius * Math.sin(phi)],
      weight,
    };
  }
}

// Cubemap face axes in world space (right-handed, +Y up, +Z forward).
// Faces in order: +X, -X, +Y, -Y, +Z, -Z.
// For each face, (right, up, forward) map texel (s, t) in [-1, 1]^2 to
// direction normalize(forward + s*right + t*up).
const CUBE_FACES = [
  { right: [0, 0, -1], up: [0, 1, 0], forward: [1, 0, 0] },
  { right: [0, 0, 1], up: [0, 1, 0], forward: [-1, 0, 0] },
  { right: [1, 0, 0], up: [0, 0, -1], forward: [0, 1, 0] },
  { right: [1, 0, 0], up: [0, 0, 1], forward: [0, -1, 0] },
  { right: [1, 0, 0], up: [0, 1, 0], forward: [0, 0, 1] },
  { right: [-1, 0, 0], up: [0, 1, 0], forward: [0, 0, -1] },
];

// Iterate all texels of a cubemap of the given face size, yielding the
// world-space direction of each texel center and its solid angle.
// Weights are normalized so they sum to exactly 4*pi.
export function* cubemapTexels(size) {
  const raw = [];
  let totalWeight = 0;
  for (const face of CUBE_FACES) {
    for (let v = 0; v < size; v++) {
      for (let u = 0; u < size; u++) {
        const s = (2 * (u + 0.5)) / size - 1;
        const t = (2 * (v + 0.5)) / size - 1;
        const dir = [
          face.forward[0] + s * face.right[0] + t * face.up[0],
          face.forward[1] + s * face.right[1] + t * face.up[1],
          face.forward[2] + s * face.right[2] + t * face.up[2],
        ];
        const len = Math.hypot(dir[0], dir[1], dir[2]);
        dir[0] /= len;
        dir[1] /= len;
        dir[2] /= len;
        // Differential solid angle: dW = ds*dt / (1 + s^2 + t^2)^(3/2),
        // with ds = dt = 2/size.
        const norm = 1 + s * s + t * t;
        const weight = (4 / (size * size)) / (norm * Math.sqrt(norm));
        totalWeight += weight;
        raw.push({ direction: dir, weight });
      }
    }
  }
  const scale = (4 * Math.PI) / totalWeight;
  for (const sample of raw) {
    yield { direction: sample.direction, weight: sample.weight * scale };
  }
}
