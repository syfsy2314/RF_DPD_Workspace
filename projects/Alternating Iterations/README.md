# Alternating Iterations —— 面向非平坦线性相位反馈链路的结构化幅度自校准

## 1. 项目定位

学位论文第四章 4.2 节的算法与实验载体。解决的问题是：

> 第三章的 LD-PLS 辨识要求模型侧带限覆盖算子与实际反馈链路匹配；当模拟反馈滤波器**带内幅度非平坦**（且无法预先精确测量）时，传统覆盖（仅用已知数字滤波器 h2）产生结构化回归矩阵失配，PA 参数辨识与全带重建出现偏差。

**核心思路**：利用第 0 轮传统覆盖得到的粗略 PA 全带输出，从带限反馈观测中估计一个**实系数、偶对称、零相位**的结构化幅度算子 `A(a)`，修正模型侧覆盖为 `h_model = conv(hA, h2)`，再用 LD-PLS 重新辨识。与第三章的关系：**不是替代，是滤波器覆盖前提的推广**（匹配覆盖 → 失配覆盖 + 自校准）。

观测模型（理论文档式(1)）：

```
y_bf = H_L · A_★ · Ψ · θ_★ + η
```

## 2. 目录结构

```
Alternating Iterations/
├── main.m                          主实验入口：通带纹波扫描 + 幅度自校准对照
├── estimate_amplitude_operator.m   核心函数：实偶对称零相位幅度算子估计
├── run_amplitude_alternating_iterations.m  交替迭代驱动函数
├── tests/
│   └── test_amplitude_alternating_iterations_mg.m   MG 对照的单次迭代测试
├── results/
│   ├── result.md                   纹波扫描结果汇总表（最终结果）
│   ├── NMSE_vs_ripple_BW80.png     纹波扫描 NMSE 曲线（BW=80）
│   └── plot_NMSE.m                 结果绘图脚本
├── doc/
│   ├── 幅度自校准理论_统一符号与严格论证 (2).md   理论主文档（命题 1–5）
│   ├── amplitude_self_calibration_theory.md      理论早期版本
│   ├── 晏嘉文 2026年8月5号.docx                   进展记录
│   └── main_legacy.m               归档：早期线性响应复增益对照脚本
└── README.md
```

## 3. 运行入口

`main.m` 是唯一主入口，做**通带纹波扫描**：逐点扫描 H1 通带纹波 `Rp_list`，每个纹波点运行三组对照——

| 组 | 覆盖算子 | 含义 |
|---|---|---|
| **Ref** | `conv(h1, h2)`（真实复合滤波器） | 性能上限（完美已知反馈链路） |
| **Initial** | 仅 `h2`（传统覆盖） | 未校准基线 = "迭代 0 次" |
| **Iteration** | 一次自校准修正后的 `hModel` | 方法贡献 |

输出：逐点打印、汇总表、NMSE-纹波三曲线图。数据路径：

```matlab
D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat
D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat
```

## 4. H1 的构造：固定阶数 + 通带纹波独立可调

扫描实验的关键设计是把 H1 分解为「平坦低通 × 余弦波纹」两部分，使**阶数与波纹解耦**（避免等波纹 FIR 扫纹波时阶数跟着变、以及频率采样法逼近误差淹没波纹的问题）：

- **平坦低通 `h_flat`**：`firpm` 一次设计、阶数固定，通带纹波 0.01 dB（远小于最小扫描值 0.1 dB），只负责截止与过渡带；
- **通带波纹 `h_ripple`**：余弦幅度波纹 `1 + ε·cos(2πf·τ/Fs)` 的**精确时域表示**是三抽头零相位 FIR `[ε/2, 1, ε/2]`（左右各隔 τ 个样本），零逼近误差；
- **合成**：`h1 = conv(h_flat, h_ripple)`，对应频域幅度相乘。

波纹幅度与 Rp 一一对应：`ε = (10^(Rp/20) − 1) / (10^(Rp/20) + 1)`；`τ = round(nRipple·Fs/(BW+5))` 控制通带内波纹周期数。两者均为线性相位 FIR，群延迟为常数，会被 `align_circular_signal` 吸收，故实验测的是纯幅度纹波。

## 5. 各函数说明

### `estimate_amplitude_operator.m`（核心）

```matlab
function [hModel, hA, a, fitNmseDb] = ...
    estimate_amplitude_operator(yCoarse, yFeedback, h2, R, eta)
```

由粗略 PA 全带输出估计实偶对称零相位幅度算子。流程对应理论文档式(10)–(13)：

| 步骤 | 内容 |
|---|---|
| 构造回归矩阵 Z | 先 `filter(h2,1,yCoarse)`，再按 `z_r(n)=z0(n−r)+z0(n+r)` 生成 R+1 列 |
| 目标对齐 | 反馈信号对 Z 的第 1 列做**一次**公共时延/常相位对齐（禁止逐列对齐） |
| 有效区间 | 首尾各裁 `len(h2)−1+R` 个保护样本 |
| 岭回归 | 实参数岭回归，`λ = eta·trace(gram)/n` |
| 尺度归一化 | 强制中心频率增益 `A(0)=1`，消除与 θ 的双线性尺度歧义 |

### `run_amplitude_alternating_iterations.m`（交替迭代驱动）

```matlab
function [yResults, gainResults] = run_amplitude_alternating_iterations( ...
    iterationCount, yInitial, gainInitial, yFeedback, h2, ...
    x, K, M, X, S, lamda, operatorR, operatorEta)
```

外层交替：`θ̂^(i) → 估计幅度算子 → 更新 hModel → 反馈对齐 → LD-PLS 重辨识 θ̂^(i+1)`。输出矩阵第 1 列为初始结果，第 i+1 列为第 i 次迭代结果。

### `tests/test_amplitude_alternating_iterations_mg.m`

以 MG（`get_MG`）代替 LD-PLS 的单次交替测试，验证自校准框架对辨识器的无关性。⚠️ 数据路径仍指向旧机器（`C:\Users\admin\Desktop\测试\第四次\`），运行前需改为 `D:\RF_DPD_Workspace\datasets\...`。

### `doc/main_legacy.m`

归档的早期脚本（线性响应估计复增益 + MG/OPLS 对照），已由 `main.m` 取代，仅留作历史备查，不作为入口。

## 6. 关键参数

| 参数 | 当前值 | 含义 |
|---|---|---|
| `BW`, `Fs` | 80, 1000 | 采集带宽、采样率 |
| `Rp_list` | `0.1:0.2:2` | H1 通带纹波扫描范围（dB） |
| `nRipple` | 5 | 通带内波纹周期数 |
| `N1` | firpmord 自动（约 815） | 平坦低通固定阶数 |
| `operator_R` / `operator_eta` | 250 / 1e-4 | 幅度算子单边阶数 / 岭正则 |
| `iteration_count` | 1 | 外层交替次数 |
| `S`, `lamda` | 0.3, 1 | LD-PLS 潜变量比例与岭正则 |
| `K`, `M` | 7, 5 | MP 非线性阶数、记忆深度 |

## 7. 实验结果

`results/result.md` 记录最终纹波扫描结果（BW=80，标称纹波 0.1–1.9 dB，步进 0.2 dB）。随反馈链路通带纹波从 0.1 dB 增大到 1.9 dB，传统覆盖（Initial）全带恢复 NMSE 由 -38.32 dB 退化至 -32.45 dB，而一次幅度自校准后仅由 -38.35 dB 降至 -36.18 dB，与性能上限 Ref（约 -38.45 dB）的差距由约 6.1 dB 收窄至约 2.3 dB，验证了幅度自校准的有效性。NMSE-纹波曲线见 `results/NMSE_vs_ripple_BW80.png`。

## 8. 理论文档指引

`doc/幅度自校准理论_统一符号与严格论证 (2).md` 为准，包含：

- 命题 1：传统覆盖估计误差的精确分解（回归矩阵失配项）
- 命题 2：约束岭回归闭式解（式(13)）
- 命题 3：首次算子估计扰动界（式(16)，σ_min 条件）
- 命题 4：给定算子时 PA 参数误差分解（式(20)(21)）
- 命题 5：交替更新只保证反馈域目标不增（非全带 NMSE 单调）
- §14 可声称 / 不可声称清单；§15 创新点表述

## 9. 已知边界与后续方向

幅度自校准主线（通带纹波扫描 + 一次自校准验证）已完成；以下为已知边界与后续方向：

- **A(ω) ≥ 0 未验证**：当前无正值约束，估计算子局部为负时是"实偶零相位残余算子"而非严格幅度（理论 §5.1 的两选一待落地）
- **实测未接入**：`datasets/测试/2026.6.10/test2/{50–90}MHz/` 已有实测记录，尚未进入本框架
- **非线性相位**：本方法严格适用于"幅度非平坦、相位近似线性"；非线性相位是后续方向
