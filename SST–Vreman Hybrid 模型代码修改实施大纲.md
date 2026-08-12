# SST–Vreman Hybrid 模型代码修改实施大纲

## 1. 修改目标

在现有 **SST / SST-DDES / SST-IDDES** CFD 代码基础上，实现第一版 **SST–Vreman Hybrid RANS/LES 模型**。

总体思想为：

- RANS 区域由标准 \(k-\omega\) SST 模型负责；
- LES 区域由 Vreman SGS 模型直接提供亚格子涡黏度；
- 使用基于 SST-DDES 屏蔽函数和网格分辨率判据构造连续混合权重 \(\gamma\)；
- 第一版暂时不修改 SST 的 \(k\)、\(\omega\) 输运方程；
- 仅将进入动量方程、能量方程和黏性通量计算的湍流黏度由原来的 \(\mu_t^{SST}\) 修改为混合黏度 \(\mu_t^H\)。

最终模型：

\[
\nu_t^H
=
(1-\gamma)\nu_t^{SST}
+
\gamma\nu_t^{Vreman}.
\]

必须满足两个严格极限：

\[
\gamma=0
\Rightarrow
\nu_t^H=\nu_t^{SST},
\]

\[
\gamma=1
\Rightarrow
\nu_t^H=\nu_t^{Vreman}.
\]

---

## 2. 修改原则

Codex 修改代码时遵循以下原则。

1. **尽可能少修改现有 SST 模块。**
2. 不删除现有 SST、DES、DDES 或 IDDES 实现。
3. 新模型作为独立选项加入，例如：
   - `SST`
   - `SST_DDES`
   - `SST_IDDES`
   - `SST_VREMAN`
4. 第一阶段不改变标准 SST 的 \(k\) 方程和 \(\omega\) 方程。
5. SST 模型仍在全计算域计算一个背景 RANS 黏度：
   \[
   \nu_t^{SST}.
   \]
6. Vreman 单独计算：
   \[
   \nu_t^V.
   \]
7. 新建混合权重：
   \[
   \gamma.
   \]
8. 动量方程实际使用：
   \[
   \nu_t^H.
   \]
9. 所有新增量必须能够输出用于后处理和调试。
10. 所有除法、平方根和双曲函数计算必须加入数值保护。

---

# 3. 第一阶段：梳理现有代码结构

在修改任何模型公式之前，先检查并记录以下代码位置。

### 3.1 找出 SST 模型实现位置

确认：

- \(k\) 输运方程；
- \(\omega\) 输运方程；
- SST 的 \(F_1\)、\(F_2\)；
- SST 湍流黏度计算；
- production limiter；
- turbulent diffusion；
- wall distance；
- \(S_{ij}\) 或 strain-rate magnitude；
- \(\Omega_{ij}\) 或 vorticity magnitude。

重点找到当前代码中：

\[
\nu_t^{SST}
=
\frac{a_1k}
{\max(a_1\omega,F_2S)}
\]

或其等价实现。

记录当前变量名，例如可能是：

```cpp
mut
nut
eddyViscosity
turbulentViscosity
mu_t
```

不要在尚未确定其用途时直接覆盖该变量。

---

## 3.2 找出湍流黏度最终进入哪些方程

搜索所有使用：

```cpp
mu_t
nut
eddyViscosity
```

的位置。

至少确认：

- 动量方程黏性应力；
- 能量方程；
- 热通量；
- 有效黏度；
- turbulence production；
- 壁面函数；
- SST 自身扩散项。

必须区分：

### 背景 SST 黏度

定义为：

```cpp
nu_t_sst
```

用于 SST 模型自身。

### 实际混合黏度

定义为：

```cpp
nu_t_hybrid
```

用于主控制方程。

**第一版不能粗暴地把代码中所有 `nu_t_sst` 全部替换为 `nu_t_hybrid`。**

否则 SST 的输运方程本身也会被悄悄改变。

---

# 4. 第二阶段：新增 Vreman SGS 模块

建议新建独立函数，例如：

```cpp
double ComputeVremanViscosity(...);
```

或者相应的类/模块。

---

## 4.1 计算速度梯度张量

定义：

\[
\alpha_{ij}
=
\frac{\partial u_j}{\partial x_i}.
\]

程序中必须得到完整的 \(3\times3\) 速度梯度：

\[
\boldsymbol{\alpha}
=
\begin{bmatrix}
\partial u/\partial x &
\partial v/\partial x &
\partial w/\partial x
\\
\partial u/\partial y &
\partial v/\partial y &
\partial w/\partial y
\\
\partial u/\partial z &
\partial v/\partial z &
\partial w/\partial z
\end{bmatrix}.
\]

优先复用求解器已有的 velocity-gradient reconstruction。

不要为了 Vreman 单独使用低阶梯度算法，否则 SGS 模型可能与主空间离散不一致。

---

## 4.2 获取三个方向网格尺度

定义：

\[
\Delta_1,\qquad
\Delta_2,\qquad
\Delta_3.
\]

对于结构六面体网格，优先采用局部三个方向物理尺寸：

\[
\Delta_1=h_1,\qquad
\Delta_2=h_2,\qquad
\Delta_3=h_3.
\]

不要第一版直接全部写成：

\[
V^{1/3}.
\]

如果代码采用曲线坐标或非正交网格，需要检查当前程序如何定义 DES/LES 网格尺度，并尽可能利用已有几何信息。

新增输出：

```cpp
delta1
delta2
delta3
```

用于验证。

---

# 5. 第三阶段：实现 Vreman 张量

计算：

\[
\beta_{ij}
=
\sum_{m=1}^{3}
\Delta_m^2
\alpha_{mi}\alpha_{mj}.
\]

建议 Codex 显式实现，不要一开始进行过度优化。

伪代码：

```cpp
for (int i = 0; i < 3; ++i)
{
    for (int j = 0; j < 3; ++j)
    {
        beta[i][j] = 0.0;

        for (int m = 0; m < 3; ++m)
        {
            beta[i][j] +=
                delta[m] * delta[m]
                * alpha[m][i]
                * alpha[m][j];
        }
    }
}
```

然后计算：

\[
A_\alpha
=
\alpha_{ij}\alpha_{ij}
=
\sum_{i,j}\alpha_{ij}^2.
\]

再计算：

\[
B_\beta
=
\beta_{11}\beta_{22}-\beta_{12}^2
+
\beta_{11}\beta_{33}-\beta_{13}^2
+
\beta_{22}\beta_{33}-\beta_{23}^2.
\]

注意 C++ 数组从 0 开始，因此对应：

```cpp
Bbeta =
    beta[0][0] * beta[1][1]
    - beta[0][1] * beta[0][1]
    + beta[0][0] * beta[2][2]
    - beta[0][2] * beta[0][2]
    + beta[1][1] * beta[2][2]
    - beta[1][2] * beta[1][2];
```

---

# 6. 第四阶段：Vreman 黏度计算

定义：

\[
\nu_t^V
=
C_V
\sqrt{
\frac{B_\beta}
{A_\alpha}
}.
\]

第一版：

\[
C_V=0.07.
\]

写成输入参数，不要硬编码：

```text
VREMAN_CONSTANT = 0.07
```

建议配置文件能够修改。

---

## 6.1 数值保护

必须处理：

\[
A_\alpha\rightarrow0
\]

以及由于浮点误差出现：

\[
B_\beta<0.
\]

实现类似：

```cpp
const double epsA = 1.0e-30;

Aalpha = max(Aalpha, epsA);
Bbeta  = max(Bbeta, 0.0);

if (Bbeta <= epsB)
{
    nu_t_vreman = 0.0;
}
else
{
    nu_t_vreman =
        Cv * sqrt(Bbeta / Aalpha);
}
```

不要使用人为较大的：

```cpp
nu_t_vreman_min
```

因为 Vreman 的重要性质之一就是允许：

\[
\nu_t^V=0.
\]

---

# 7. 第五阶段：保留 SST 背景黏度

标准 SST 继续计算：

\[
\boxed{
\nu_t^{SST}
=
\frac{a_1k}
{\max(a_1\omega,F_2S)}
}
\]

将其保存成独立变量：

```cpp
nu_t_sst
```

或者：

```cpp
mu_t_sst = rho * nu_t_sst;
```

不能直接被 Vreman 覆盖。

---

# 8. 第六阶段：构造 SST RANS 长度尺度

定义：

\[
\boxed{
l_R
=
\frac{\sqrt{k}}
{C_\mu\omega}
}
\]

第一版：

\[
C_\mu=0.09.
\]

加入保护：

\[
k\ge k_{\min},
\]

\[
\omega\ge\omega_{\min}.
\]

例如：

```cpp
double kSafe     = max(k, kMin);
double omegaSafe = max(omega, omegaMin);

l_rans =
    sqrt(kSafe) /
    (Cmu * omegaSafe);
```

输出：

```cpp
l_rans
```

---

# 9. 第七阶段：实现 SST-DDES 型 shielding

第一版建议采用：

\[
\boxed{
f_d
=
1-\tanh[(20r_d)^3]
}
\]

其中：

\[
r_d
=
\frac{\nu+\nu_t^{SST}}
{
\kappa^2d_w^2
G
},
\]

速度梯度尺度可根据现有 SST-DDES 实现采用：

\[
G=
\sqrt{
\frac12
(S^2+\Omega^2)
}.
\]

如果当前代码已有 SST-DDES：

**优先直接复用现有 \(r_d\) 和 \(f_d\) 实现，不重新写一套。**

这是最重要的代码复用点之一。

必须检查：

```cpp
fd ≈ 0
```

是否对应附着边界层，

而：

```cpp
fd ≈ 1
```

是否对应解除 RANS shielding 的区域。

输出：

```cpp
r_d
f_d
```

---

# 10. 第八阶段：定义 LES 网格判据

第一版定义一个标量混合尺度：

\[
\Delta_H.
\]

最简单版本先采用：

\[
\boxed{
\Delta_H=h_{\max}
}
\]

其中：

\[
h_{\max}
=
\max(h_1,h_2,h_3).
\]

然后：

\[
\boxed{
l_G
=
C_{DES}\Delta_H
}
\]

`C_DES` 暂时沿用当前 SST-DDES 的值。

不要把 Vreman 自身的 \(\Delta_i\) 和这里的 \(\Delta_H\) 混为同一个变量。

二者功能不同：

- \(\Delta_i\)：Vreman SGS 张量计算；
- \(\Delta_H\)：判断是否具备 LES 分辨率。

建议变量命名：

```cpp
delta_vreman[3]
delta_hybrid
l_grid
```

---

# 11. 第九阶段：构造网格分辨率函数

第一版定义：

\[
\boxed{
f_{res}
=
\max
\left(
0,
1-\frac{l_G}{l_R}
\right)
}
\]

再限制：

\[
0\le f_{res}\le1.
\]

代码：

```cpp
double ratio =
    l_grid / max(l_rans, epsLength);

f_res =
    max(0.0, 1.0 - ratio);

f_res =
    min(1.0, f_res);
```

物理检查：

### 网格太粗

\[
l_G\ge l_R
\]

应得到：

\[
f_{res}=0.
\]

### 网格很细

\[
l_G\ll l_R
\]

应得到：

\[
f_{res}\rightarrow1.
\]

输出：

```cpp
f_res
```

---

# 12. 第十阶段：构造总混合权重

定义：

\[
\boxed{
\gamma=f_d f_{res}
}
\]

同时限制：

\[
0\le\gamma\le1.
\]

代码：

```cpp
gamma = f_d * f_res;
gamma = max(0.0, min(1.0, gamma));
```

验证三个极限。

### 附着边界层

\[
f_d\approx0
\]

因此：

\[
\gamma\approx0.
\]

### 分离区 + 粗网格

\[
f_d\approx1,
\qquad
f_{res}\approx0
\]

因此：

\[
\gamma\approx0.
\]

### 分离区 + LES 网格

\[
f_d\approx1,
\qquad
f_{res}\approx1
\]

因此：

\[
\gamma\approx1.
\]

---

# 13. 第十一阶段：构造实际混合涡黏度

定义：

\[
\boxed{
\nu_t^H
=
(1-\gamma)\nu_t^{SST}
+
\gamma\nu_t^{Vreman}
}
\]

代码：

```cpp
nu_t_hybrid =
    (1.0 - gamma) * nu_t_sst
    + gamma * nu_t_vreman;
```

然后：

\[
\boxed{
\mu_t^H
=
\rho\nu_t^H
}
\]

注意不要加入：

```cpp
max(nu_t_sst, nu_t_vreman)
```

也不要加入：

```cpp
min(...)
```

第一版严格测试线性 blending。

---

# 14. 第十二阶段：明确哪些方程使用哪一个黏度

这是整个代码修改最重要的部分之一。

## 主流动控制方程使用

\[
\boxed{
\nu_t^H
}
\]

包括：

- Reynolds/SGS stress；
- 动量黏性通量；
- 总有效黏度；
- LES/RANS 对主流场的实际作用。

---

## SST 自身第一版继续使用

\[
\boxed{
\nu_t^{SST}
}
\]

包括：

- \(k\) production；
- \(\omega\) production；
- SST diffusion；
- SST cross diffusion；
- SST limiter；
- shielding sensor 中需要的背景 RANS 黏度。

即第一版明确采用：

```text
Background SST:
    nu_t_sst

Flow equations:
    nu_t_hybrid
```

不要让 Codex 自动全局替换变量。

---

# 15. 能量方程和湍流热通量

如果代码求解可压缩 N-S，需要检查能量方程。

原先如果：

\[
k_t
=
\frac{c_p\mu_t}
{Pr_t},
\]

则主流能量方程建议使用：

\[
\boxed{
\mu_t=\mu_t^H
}
\]

因此：

\[
k_t^H
=
\frac{c_p\mu_t^H}
{Pr_t}.
\]

第一版暂时保持原湍流 Prandtl 数。

不要同时引入新的 SGS Prandtl 模型，否则无法判断 Vreman 的单独作用。

---

# 16. 第十三阶段：新增诊断变量

至少必须能够输出：

\[
\nu_t^{SST}
\]

\[
\nu_t^{Vreman}
\]

\[
\nu_t^H
\]

\[
\gamma
\]

\[
f_d
\]

\[
f_{res}
\]

\[
l_R
\]

\[
l_G
\]

\[
A_\alpha
\]

\[
B_\beta
\]

以及：

\[
\frac{\nu_t^{Vreman}}
{\nu_t^{SST}+\varepsilon}.
\]

建议 Tecplot 输出名称：

```text
Nut_SST
Nut_Vreman
Nut_Hybrid
Hybrid_Gamma
DDES_fd
GridResolutionFactor
L_RANS
L_GRID
Vreman_Aalpha
Vreman_Bbeta
Nut_Vreman_to_SST
```

这些量对判断模型是否正确工作至关重要。

---

# 17. 第十四阶段：增加运行时开关

建议配置文件加入：

```text
TURBULENCE_MODEL = SST_VREMAN_HYBRID

VREMAN_CONSTANT = 0.07

HYBRID_USE_SHIELDING = 1
HYBRID_USE_RESOLUTION_SENSOR = 1

HYBRID_OUTPUT_DIAGNOSTICS = 1
```

不要把新模型直接覆盖现有 SST-IDDES。

必须保证：

```text
TURBULENCE_MODEL = SST
```

时结果与修改前完全一致。

---

# 18. 第十五阶段：单元测试

在跑复杂算例之前，先写小型测试验证 Vreman 数学实现。

## Test A：均匀速度

\[
\nabla\mathbf u=0
\]

要求：

\[
A_\alpha=0
\]

最终：

\[
\boxed{
\nu_t^V=0
}
\]

不得产生 NaN。

---

## Test B：简单速度梯度场

人为输入一个固定：

\[
\alpha_{ij}.
\]

独立使用 Python/手算计算：

\[
\beta_{ij},
\quad
B_\beta,
\quad
\nu_t^V.
\]

要求 C++ 结果与理论值一致。

---

## Test C：检查坐标索引

特别检查：

\[
\alpha_{ij}
=
\partial u_j/\partial x_i
\]

是否因为代码中的矩阵存储习惯被写成了：

\[
\partial u_i/\partial x_j.
\]

这是 Vreman 实现中非常容易出现的错误。

---

# 19. 第十六阶段：退化性测试

## 测试 1：强制 \(\gamma=0\)

设置：

```cpp
gamma = 0.0;
```

要求：

\[
\nu_t^H=\nu_t^{SST}.
\]

完整 CFD 结果必须与标准 SST 在数值误差范围内一致。

如果不一致，说明修改影响了 SST 方程其他部分。

---

## 测试 2：强制 \(\gamma=1\)

设置：

```cpp
gamma = 1.0;
```

要求：

\[
\nu_t^H=\nu_t^V.
\]

用于验证纯 Vreman 分支。

---

## 测试 3：关闭 Vreman

临时设置：

\[
\nu_t^V=\nu_t^{SST}.
\]

则无论：

\[
\gamma
\]

是多少，都应该得到：

\[
\nu_t^H=\nu_t^{SST}.
\]

这是验证 blending 实现的很好办法。

---

# 20. 第十七阶段：建议验证算例顺序

不要直接从多段翼型开始。

### 算例 1：Channel Flow

目的：

- 检查 RANS shielding；
- 检查近壁 \(\gamma\)；
- 确认 SST 边界层未被破坏。

重点看：

\[
\gamma(y),
\]

\[
\nu_t^{SST}(y),
\]

\[
\nu_t^V(y).
\]

---

### 算例 2：Backward-Facing Step

这是第一项真正重要的 hybrid 验证。

重点区域：

1. 入口附着边界层；
2. 台阶分离点；
3. 自由剪切层；
4. 回流区；
5. 再附着区。

期望：

入口边界层：

\[
\gamma\approx0.
\]

分离剪切层：

\[
\gamma\uparrow.
\]

LES 分辨区域：

\[
\gamma\rightarrow1.
\]

再附着近壁区域：

\[
\gamma
\]

应该重新减小。

重点比较：

- 再附着长度；
- \(C_f\)；
- 平均速度；
- Reynolds stress；
- resolved turbulent kinetic energy。

---

### 算例 3：周期均匀湍流或衰减湍流

目的：

验证：

\[
\boxed{\text{纯Vreman LES行为}}
\]

避免 SST shielding 干扰。

---

### 算例 4：RA16SC1 / 30P30N 等高升力构型

只有前面验证完成后再进入。

重点检查：

- slat cusp shear layer；
- slat cavity；
- main-element wake；
- flap separation；
- RANS/LES 接口；
- 声学相关高频耗散。

---

# 21. 第一版暂时不要做的修改

Codex 第一轮修改时明确禁止以下操作：

1. 不修改 SST \(k\) 方程耗散项；
2. 不修改 \(\omega\) 方程；
3. 不加入 \(f_e\)；
4. 不加入完整 IDDES WMLES 分支；
5. 不重新标定 \(C_V\)；
6. 不修改湍流 Prandtl 数；
7. 不修改 WENO/MUSCL 数值格式；
8. 不修改时间推进；
9. 不引入 dynamic Vreman；
10. 不将 Vreman 黏度强制设置一个正的下限；
11. 不把所有 SST 黏度变量全局替换成 hybrid 黏度；
12. 不同时优化代码性能。

第一阶段目标只有一个：

\[
\boxed{
\text{确认SST-RANS/Vreman-LES混合思想本身是否工作}
}
\]

---

# 22. 第二阶段可能的模型升级

只有 Version 1 验证完成后，再考虑。

## 22.1 修改 \(k\) production

可以尝试：

\[
\boxed{
P_k^H
=
(1-\gamma)P_k^{SST}
}
\]

使 LES 区：

\[
\gamma\rightarrow1
\]

时：

\[
P_k^{SST}\rightarrow0.
\]

这样 SST 的背景 \(k\) 不会在 LES 区持续增长。

---

## 22.2 修改 \(\omega\) production

必须与 \(k\) production 修改配套研究，不能单独随意乘：

\[
1-\gamma.
\]

需要检查 SST production ratio 和近壁渐近性质。

---

## 22.3 引入接口保护函数

如果发现：

\[
\nu_t^{SST}
\rightarrow
\nu_t^V
\]

下降过快，而 resolved turbulence 尚未建立，可以增加类似 IDDES elevating function 的：

\[
f_e.
\]

但这一项必须根据实际的：

\[
\nu_t^H,
\quad
k_{resolved},
\quad
\gamma
\]

结果再设计。

---

# 23. Codex 实现顺序

严格按照以下顺序提交修改：

1. **定位 SST 湍流黏度及其所有用途。**
2. **新增独立 Vreman 模块。**
3. **完成 Vreman 单元测试。**
4. **新增 `nu_t_sst` 和 `nu_t_vreman`，暂不改变主求解器。**
5. **复用/实现 SST-DDES 的 \(f_d\)。**
6. **实现 \(l_R\)、\(l_G\)、\(f_{res}\)。**
7. **实现 \(\gamma\)。**
8. **实现 `nu_t_hybrid`。**
9. **仅将主流动方程改用 `nu_t_hybrid`。**
10. **保持 SST 输运方程使用 `nu_t_sst`。**
11. **增加诊断输出。**
12. **进行 \(\gamma=0\) 退化测试。**
13. **进行 \(\gamma=1\) 测试。**
14. **运行简单 CFD 验证。**
15. **确认无误后再进行性能优化。**

每一步完成后都应保持代码可以编译运行。

不要一次性修改整个湍流模型。

---

# 24. Codex 每次修改后必须报告

每一次代码修改，都要求 Codex 给出：

### 修改文件

例如：

```text
src/turbulence/SST.cpp
src/turbulence/Vreman.cpp
src/flux/ViscousFlux.cpp
```

### 修改函数

例如：

```text
ComputeSSTViscosity()
ComputeVremanViscosity()
ComputeHybridViscosity()
```

### 新增变量

例如：

```text
nu_t_sst
nu_t_vreman
nu_t_hybrid
f_d
f_res
gamma
```

### 公式对应关系

明确指出每一段代码对应：

\[
\nu_t^V,
\quad
f_d,
\quad
f_{res},
\quad
\gamma,
\quad
\nu_t^H
\]

中的哪个公式。

### 风险

特别说明：

- 是否改变了标准 SST 行为；
- 是否改变了 \(k,\omega\) 方程；
- 是否改变了 wall boundary condition；
- 是否改变了原 SST-DDES/IDDES。

---

# 25. 最终验收标准

第一版 SST–Vreman Hybrid 完成后至少满足：

### 数学正确性

\[
B_\beta\ge0
\]

数值上得到保护；

\[
\nu_t^V\ge0;
\]

\[
0\le\gamma\le1;
\]

\[
\nu_t^H\ge0.
\]

### RANS 极限

\[
\gamma=0
\]

时严格恢复标准 SST。

### LES 极限

\[
\gamma=1
\]

时动量方程严格使用 Vreman SGS。

### 边界层保护

正常附着边界层：

\[
f_d\approx0
\]

从而：

\[
\gamma\approx0.
\]

### 网格保护

粗网格区域：

\[
f_{res}\approx0
\]

不得仅因流动分离就直接切入 Vreman LES。

### 可诊断性

必须能够输出：

\[
f_d,\quad
f_{res},\quad
\gamma,
\quad
\nu_t^{SST},
\quad
\nu_t^V,
\quad
\nu_t^H.
\]

---

# 26. 第一版模型最终公式汇总

Codex 应以以下公式为第一版实现目标：

\[
\boxed{
\nu_t^{SST}
=
\frac{a_1k}
{\max(a_1\omega,F_2S)}
}
\]

\[
\boxed{
\alpha_{ij}
=
\frac{\partial u_j}{\partial x_i}
}
\]

\[
\boxed{
\beta_{ij}
=
\sum_{m=1}^{3}
\Delta_m^2
\alpha_{mi}\alpha_{mj}
}
\]

\[
\boxed{
B_\beta
=
\beta_{11}\beta_{22}-\beta_{12}^2
+
\beta_{11}\beta_{33}-\beta_{13}^2
+
\beta_{22}\beta_{33}-\beta_{23}^2
}
\]

\[
\boxed{
\nu_t^V
=
C_V
\sqrt{
\frac{B_\beta}
{\alpha_{ij}\alpha_{ij}}
}
}
\]

\[
\boxed{
l_R
=
\frac{\sqrt{k}}
{C_\mu\omega}
}
\]

\[
\boxed{
l_G
=
C_{DES}\Delta_H
}
\]

\[
\boxed{
f_{res}
=
\max
\left(
0,
1-\frac{l_G}{l_R}
\right)
}
\]

\[
\boxed{
f_d
=
1-\tanh[(20r_d)^3]
}
\]

\[
\boxed{
\gamma
=
f_df_{res}
}
\]

最终：

\[
\boxed{
\nu_t^H
=
(1-\gamma)\nu_t^{SST}
+
\gamma\nu_t^V
}
\]

第一版模型的核心原则概括为：

\[
\boxed{
\text{SST决定RANS应力}
\quad+\quad
\text{Vreman决定LES亚格子应力}
\quad+\quad
\text{DDES型shielding决定是否离开RANS}
\quad+\quad
\text{网格分辨率判据决定是否允许LES}
}
\]

实现第一版时，应优先追求**模型逻辑清晰、标准 SST 可严格恢复、各个中间量可诊断**，而不是立即追求完整 IDDES 功能或复杂的经验修正。