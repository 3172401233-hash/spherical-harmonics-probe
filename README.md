# spherical-harmonics-probe

JavaScript 二阶球谐（SH, order-2 / 9 系数）光照探针库：从天空盒方向采样生成 SH 系数，
做低频辐照度重建，支持探针与场景整体旋转，并提供能量校验。

## 约定（显式声明，不依赖渲染管线默认值）

- **世界坐标轴**：右手系，`+X` 右、`+Y` 上、`+Z` 前。所有方向为单位向量 `[x, y, z]`。
- **系数顺序**：ACN（Ambisonic Channel Number），`index = l*(l+1)+m`：

  | index | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
  |-------|---|---|---|---|---|---|---|---|---|
  | basis | Y00 | Y1-1 | Y10 | Y11 | Y2-2 | Y2-1 | Y20 | Y21 | Y22 |

- **归一化**：内部一律使用 N3D（完全归一化、正交归一基）。SN3D 通过
  `toSN3D` / `fromSN3D` 转换，每带缩放因子为 `sqrt(4π / (2l+1))`。
- **基函数**（N3D，直接以世界方向分量表示）：

  ```
  Y00  = 0.2820947918
  Y1-1 = 0.4886025119 * y        Y10 = 0.4886025119 * z        Y11 = 0.4886025119 * x
  Y2-2 = 1.0925484306 * x*y      Y2-1 = 1.0925484306 * y*z     Y20 = 0.3153915653 * (3z² - 1)
  Y21  = 1.0925484306 * x*z      Y22 = 0.5462742153 * (x² - y²)
  ```

## 使用

```js
import {
  projectEnvironment, irradianceFromRadiance, evaluate,
  cubemapTexels, fibonacciSphere, rotateCoefficients, rotationFromAxisAngle,
} from './src/index.js';

// 1. 投影：天空盒方向采样 + 立体角权重（cubemapTexels 已按纹素立体角加权，
//    权重归一化到总和 4π；也可用 fibonacciSphere(N) 做函数式天空）。
const radianceSH = projectEnvironment(sampleSkybox, cubemapTexels(64));

// 2. Lambertian 卷积：每带乘 A_l = [π, 2π/3, π/4]，得到辐照度系数。
const irradianceSH = irradianceFromRadiance(radianceSH);

// 3. 场景/探针整体旋转（主动旋转，f'(d) = f(R⁻¹ d)）。
const R = rotationFromAxisAngle([0, 1, 0], Math.PI / 4);
const rotatedSH = rotateCoefficients(irradianceSH, R);

// 4. 重建：E(n) = Σ c_i * Y_i(n)。
const e = evaluate(rotatedSH, normal);
```

## 材质侧消费（角色 / 环境材质）

CPU 侧把 `irradianceFromRadiance` 输出的 9 个 RGB 系数按 **ACN 顺序**上传为
`vec3 uSH[9]`（uniform 或探针纹理，逐像素插值时逐系数线性插值即可，SH 是线性的）。
着色器用与 CPU 完全相同的 N3D 基函数求值：

```glsl
uniform vec3 uSH[9]; // 辐照度系数，N3D，ACN 顺序，已乘 A_l

vec3 shIrradiance(vec3 n) { // n: 世界空间法线，右手系 +Y 上 +Z 前
  return uSH[0] * 0.2820947918
       + uSH[1] * (0.4886025119 * n.y)
       + uSH[2] * (0.4886025119 * n.z)
       + uSH[3] * (0.4886025119 * n.x)
       + uSH[4] * (1.0925484306 * n.x * n.y)
       + uSH[5] * (1.0925484306 * n.y * n.z)
       + uSH[6] * (0.3153915653 * (3.0 * n.z * n.z - 1.0))
       + uSH[7] * (1.0925484306 * n.x * n.z)
       + uSH[8] * (0.5462742153 * (n.x * n.x - n.y * n.y));
}

// 漫反射项：albedo / π * shIrradiance(worldNormal)
```

注意事项（对应"切昼夜后转动机位亮度偏掉"的典型根因）：

- **顺序**：上传 / 下载数组必须同为 ACN；Y1-1/Y10/Y11 分别对应 y/z/x，不是 x/y/z。
- **归一化**：若引擎环境光槽位是 SN3D（如部分引擎的 ambient SH），用 `toSN3D`
  转换后再上传，不要混用两套基常量。
- **A_l 卷积**：材质直接重建漫反射时，系数必须已乘 `[π, 2π/3, π/4]`；
  若上传的是 radiance 系数而着色器里又乘一次，会导致整体亮度随视角偏移。
- **旋转**：天空盒 / 场景旋转时用 `rotateCoefficients`（主动旋转 `f'(d)=f(R⁻¹d)`）
  同步旋转探针系数；法线在世界空间求值，坐标系不一致会表现为"转动机位亮度漂"。

## 测试

```
node --test
```

覆盖：ACN/N3D 基函数取值与正交归一性、常量天空（系数、πL 辐照度、4π|L|² 能量）、
cubemap 立体角权重（总和 4π）、单侧 delta 光源（峰值方向、1.0625/0.0625 前后亮度、
9E0²/4π 能量）、旋转后的非对称天空（旋转系数 vs 重投影一致、能量不变、
重建一致性 `E'(Rd)=E(d)`）、SN3D 往返转换。
