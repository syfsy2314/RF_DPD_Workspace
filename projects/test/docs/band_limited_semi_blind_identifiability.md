# 带限观测下结构化半盲辨识的可辨识性分析

> 研究方向理论骨架（第一版）
> 核心命题：在带限观测 + 采集路径 IIR 滤波器盲未知下，联合辨识功率放大器（PA）系数与滤波器相位色散，恢复全带信号。

---

## 1. 研究定位

**一句话**：带限观测下，输入已知、PA 基函数已知、采集路径 IIR 滤波器（有理结构）系数盲未知的级联系统辨识，用于恢复全带信号。

**三条腿**：

- **带限观测**：反馈/采集路径只保留 $[-B,B]$，带外被截断；
- **结构化半盲**：输入 $x$ 已知、PA 基函数 $\psi_i$ 已知，只有 PA 系数 $\theta$ 与滤波器 $H$ 未知（不是纯盲）；
- **级联（Hammerstein 型）**：非线性 PA 在前、线性滤波器 $H$ 在后。

**区别于纯盲辨识**：现有盲 Hammerstein 辨识（Bai 2002、Yu 2014）靠"过采样 + 循环平稳"多拿信息，且输入未知；本问题输入已知，走"结构化半盲"路线，不依赖过采样。

---

## 2. 问题定义

$$
y=\mathcal{B}\big[H(z)\,\text{PA}(x)\big],\qquad \text{PA}(x)=\sum_i \theta_i\,\psi_i(x)
$$

- $x$：已知基带输入；
- $\psi_i$：已知 PA 基函数（记忆多项式）；
- $\theta_i$：PA 系数（盲未知）；
- $H(z)=B(z)/A(z)$：采集路径 IIR 滤波器（盲未知，仅知其为有理函数，有相位色散）；
- $\mathcal{B}$：带限观测算子（理想低通，只留 $[-B,B]$）。

**目标**：恢复全带 $\hat{s}=\sum_i\hat\theta_i\,\psi_i(x)\approx s$，服务于 DPD。

---

## 3. 核心机制：循环谱给出相位差

循环谱关系：

$$
S_y^\alpha(f)=H(f)\,H^*(f-\alpha)\,S_x^\alpha(f)
$$

输入循环谱 $S_x^\alpha(f)$ 已知（调制方式已知）时：

$$
\Delta(f,\alpha)=\angle H(f)-\angle H(f-\alpha)=\angle\!\left(\frac{S_y^\alpha(f)}{S_x^\alpha(f)}\right)
$$

**循环谱给出的是相位差，不是相位本身**。通过相位差累加（相位差链）重建完整相位：

$$
\angle H(f)=\angle H(f_0)+\sum_m \Delta(f_m,\alpha_m)
$$

- 普通功率谱只含幅度 $|H|^2$，相位盲；
- 循环谱含 $H(f)H^*(f-\alpha)$，相位差被显式编码 → 非最小相位（含色散）可辨识。

---

## 4. 可辨识性定理

**设定**：$y=\mathcal{B}[H(z)\,\text{PA}(x)]$，$H(z)$ 为稳定因果有理（IIR）滤波器。

**假设 A1–A4**：

- **A1（结构）**：$H$ 的零点允许非最小相位（允许相位色散）。
- **A2（循环平稳）**：输入 $x$ 的循环频率 $\alpha=k/T_s$（$k\in\mathbb{Z}$），最小正循环频率 $\alpha_{\min}=1/T_s$。
- **A3（带限）**：观测只保留 $|f|\le B$，且 $B<B_s$。
- **A4（充分激励）**：对所有满足 $|f|\le B$ 且 $|f-\alpha|\le B$ 的 $(\alpha,f)$，循环谱 $S_x^\alpha(f)\ne 0$。

**结论**：在 A1–A4 下，$H$ 的相位 $\angle H(f)$：

- (a) $|f|\le B-\alpha_{\min}$：可唯一辨识（至常数相位）；
- (b) $B-\alpha_{\min}<|f|\le B$：不可辨识（边缘盲区，宽度 $\alpha_{\min}$）；
- (c) $|f|>B$：不可辨识（带外截断）。

**证明思路（相位差链）**：

1. 循环谱给相位差 $\Delta(f,\alpha)$，在 $|f|\le B$ 且 $|f-\alpha|\le B$ 时可测；
2. 取参考点 $f_0$，构造链 $f_0\to f_0+\alpha_{\min}\to\cdots\to f$，每步两端落带内；
3. 累加得 $\angle H(f)=\angle H(f_0)+\sum_m\Delta(f_m,\alpha_m)$，唯一（差常数 $\angle H(f_0)$）；
4. 边缘 $|f|\in(B-\alpha_{\min},B]$：任何相位差步踩到带外，链断裂，不可辨识。

**三个含义**：

1. 可辨识相位带宽 $=B-\alpha_{\min}$；
2. 常数相位模糊 $=\angle H(f_0)$，等价于 $H\to e^{j\varphi}H$（复尺度模糊的相位部分）；
3. 符号率越低（$\alpha_{\min}$ 越小），盲区越窄——低符号率信号反而更易辨识相位。

---

## 5. OFDM 循环谱的实际支撑

OFDM 循环平稳有两个来源，其循环谱支撑不是"满菱形"，而是"离散 α 上的切片"。

**① CP 诱导（主来源）**：自相关以符号周期 $T_s$ 为周期，循环频率 $\alpha=m/T_s$（$m\in\mathbb{Z}$）。

$$
S_x^\alpha(f)\ne 0 \iff \alpha=\frac{m}{T_s}\ \text{且}\ |f|\le B_s,\ |f-\alpha|\le B_s
$$

支撑 = 离散 α 上的菱形切片，强度随 $|m|$ 增大衰减。

**② 导频诱导（辅助来源）**：导频子载波数据跨符号相关，循环频率 $\alpha=p\,\Delta k\,\Delta f$（$\Delta k$ 为导频子载波间隔，$\Delta f$ 为子载波间隔），支撑集中在导频附近，功率较强。

**实际可用循环频率集合**：

$$
\mathcal{C}=\left\{\frac{m}{T_s}\right\}_{CP}\ \cup\ \left\{p\,\Delta k\,\Delta f\right\}_{导频}
$$

因 $1/T_s<\Delta f$，CP 提供更密网格，是相位差链主力；导频为补充锚点。

---

## 6. 定理的 OFDM 实际化

把 A2/A4 换成 OFDM 实际结构后：

1. **可辨识相位带宽**：$B-\alpha_{\min}$，$\alpha_{\min}=1/T_s$（CP 网格最密）；
2. **高阶 CP 循环谱衰减**：$m$ 越大 $S_x^{m/T_s}(f)$ 越弱，相位差信噪比越低 → 链在高阶 $m$ 上"踩空"；
3. **导频补点**：导频诱导循环频率可修补 CP 网格弱区，把可辨识带宽往回拉。

**可写结论**：

> OFDM 循环谱支撑由 CP（$\alpha=m/T_s$，密网格但高阶衰减）与导频（$\alpha=p\Delta k\Delta f$，疏网格但功率强）共同决定；带限 $B$ 下 $H$ 相位的可辨识带宽近似 $B-1/T_s$，实际盲区宽度由高阶 CP 衰减与导频补点共同决定。

---

## 7. 与现有工作的差异（方法级）

### 7.1 三篇关键文献的方法内核

**Bai 2002（纯盲 Hammerstein–Wiener）**

- 输入未知（纯盲），仅从输出测量恢复内部变量；
- 关键手段：输出端**快速采样（过采样）→ 循环平稳 → 盲辨识**；
- 非线性结构未知、不假设白噪声。

**Yu 2014（过采样确定性 Hammerstein）**

- 输入未知（纯盲）；
- **过采样技术 → 确定性盲辨识**；
- 先盲辨识线性动态部分，再估计非线性函数。

**兰荣 2020（循环平稳 + 低采样率 DPD）**

- 输入 $x$ 已知（半盲，DPD 场景）；
- PA = 记忆多项式 $y(n)=\sum_{k,q}b_{kq}x(n-q)|x(n-q)|^{k-1}$，辨识参数 $b_{kq}$（**只 θ**）；
- 观测 = PA 输出直接欠采样（**无滤波器 H 这一层**）；
- 利用循环平稳证明低采样率可辨识。

### 7.2 方法级差异表

| 维度 | Bai 2002 | Yu 2014 | 兰荣 2020 | **本文** |
|---|---|---|---|---|
| 输入 $x$ | 未知（纯盲） | 未知（纯盲） | 已知（半盲） | **已知（半盲）** |
| 辨识对象 | θ + 线性块 | 线性块 + θ | 只 θ | **θ + H 双盲** |
| 滤波器 H | 盲辨识（过采样） | 盲辨识（过采样） | 无 | **盲未知 IIR（带限）** |
| 采样方向 | 过采样 | 过采样 | 低采样率 | **带限（少信息）** |
| 信息机制 | 过采样→循环平稳 | 过采样 | 循环平稳（低采样存活） | **结构先验 + 残存循环平稳** |

### 7.3 增量总结

兰荣已证明"低采样率 + 半盲 + **只辨识 θ**"可行；本文挑战"带限 + 半盲 + **θ 与 H 双盲**"——多辨识一个盲未知 H，且带限砍掉了信息，是三篇中最难的组合，也是空白所在。

**增量**：$H$ 盲未知 + 相位色散显式恢复 + 带限约束，三者同时成立。

---

## 8. 待解难点与假设（需验证/放松）

1. **尺度模糊（固有）**：$\theta\to\alpha\theta$、$H\to H/\alpha$ 观测不变。终点指标为 DPD 的 NMSE/ACPR 时，差尺度的 $s$ 或已够用，需明确归一化约定。
2. **带外不可辨识**：带限截断 $Y\approx 0$，带外恢复只能靠有理延拓，随色散量↑/带宽↓/噪声↑病态加剧。
3. **相位缠绕（phase wrapping）**：真实测量有 $2\pi$ 缠绕，定理为理想无噪声版，实际需相位展开。
4. **A4 实际验证**：OFDM 循环谱在哪些 $(\alpha,f)$ 上非零，需推导 CP + 导频的解析表达式坐实（高阶 CP 衰减、导频补点的定量形式）。

---

## 9. 参考文献（待精读）

1. Bai, E.-W., "A blind approach to the Hammerstein–Wiener model identification," *Automatica*, 2002.
2. Yu, C., Zhang, C., Xie, L., "A New Deterministic Identification Approach to Hammerstein Systems," *IEEE TSP*, 2014.
3. Tong, L., Xu, G., Hassibi, B., Kailath, T., "Blind channel identification based on second-order statistics," *IEEE Trans. IT*.
4. Moulines, E., et al., "Subspace methods for the blind identification of multichannel FIR filters," *IEEE Trans. SP*, 1995.
5. Gardner, W. A., *Cyclostationarity in Communications and Signal Processing*, 1994.
6. Dogancay, K., "Blind compensation of nonlinear distortion for bandlimited signals," *IEEE TCAS-I*, 2005.
7. Ma, Y., et al., "Wideband Digital Predistortion Using Spectral Extrapolation of Band-Limited Feedback Signal," *IEEE TCAS-I*, 2014.
8. Wu, X., et al., "A Low-Sampling-Rate Digital Predistortion Method Based on Inverse Filter Signal Recovery," *IEEE Trans. Broadcasting*, 2025.
9. 兰荣 等，"基于循环平稳特性的低采样频率宽带数字预失真研究," *电子与信息学报*, 2020.
10. Slepian, D., "Prolate spheroidal wave functions," *Bell Syst. Tech. J.*
