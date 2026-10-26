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
├── main.m                          主实验脚本（含 2.1.1 复增益对照）
├── amplitude_operator/             自校准核心
│   ├── amplitude_alternating_iterations.m   纹波扫描主实验（推荐入口）
│   ├── run_amplitude_alternating_iterations.m 交替迭代驱动函数
│   ├── estimate_amplitude_operator.m        幅度算子估计核心函数
│   ├── test/
│   │   └── test_amplitude_alternating_iterations_mg.m   MG 对照的单次迭代测试
│   └── result/
│       ├── result.md               纹波扫描结果汇总表（手工轮次）
│       └── plot_NMSE.m             结果绘图脚本
└── doc/
    ├── 幅度自校准理论_统一符号与严格论证 (2).md   理论主文档（命题 1–5）
    ├── amplitude_self_calibration_theory.md      理论早期版本
    └── 晏嘉文 2026年8月5号.docx                   进展记录
```

## 3. 各函数说明

### `estimate_amplitude_operator.m`（核心）

```matlab
function [hModel, hA, a, fitNmseDb] = ...
    estimate_amplitude_operator(yCoarse, yFeedback, h2, R, eta)
```

由粗略 PA 全带输出估计幅度算子。流程对应理论文档式(10)–(13)：

| 步骤 | 内容 |
|---|---|
| 构造回归矩阵 Z | 先 `filter(h2,1,yCoarse)`，再按 `z_r(n)=z0(n-r)+z0(n+r)` 生成 R+1 列（实偶对称零相位算子的时域表示） |
| 目标对齐 | 反馈信号对 Z 的第 1 列做**一次**公共时延/常相位对齐（禁止逐列对齐，避免破坏算子结构） |
| 有效区间 | 去除 FIR 瞬态与循环移位环绕样本：首尾各裁 `len(h2)-1+R` 个保护样本 |
| 岭回归 | 实参数岭回归（Gram 强制 Hermitian），正则系数 `λ = eta·trace(gram)/n`，对应理论式(12) |
| 尺度归一化 | 强制中心频率增益 `A(0)=1`（式(9)），消除与 θ 的双线性尺度歧义 |

**输出**：`a`（归一化算子系数）、`hA`（长度 2R+1 因果偶对称 FIR）、`hModel = conv(hA, h2)`（修正后的模型侧覆盖）、`fitNmseDb`（算子拟合残差，诊断用）。

### `run_amplitude_alternating_iterations.m`（交替迭代驱动）

```matlab
function [yResults, gainResults] = run_amplitude_alternating_iterations( ...
    iterationCount, yInitial, gainInitial, yFeedback, h2, ...
    x, K, M, X, S, lamda, operatorR, operatorEta)
```

外层交替：`θ̂^(i) → 估计幅度算子 → 更新 hModel → 反馈对齐 → LD-PLS 重辨识 θ̂^(i+1)`。
输出矩阵**第 1 列 = 初始恢复（迭代 0 次）**，第 i+1 列 = 第 i 次迭代结果，历史不被覆盖。

### `amplitude_alternating_iterations.m`（主实验，推荐入口）

纹波自动扫描版：H1 通带纹波 `Rp_list = 0.01:0.5:6` 逐点运行三组对照——

| 组 | 覆盖算子 | 含义 |
|---|---|---|
| **Ref** | `conv(h1, h2)`（真实复合滤波器） | 性能上限（完美已知反馈链路） |
| **Initial** | 仅 `h2`（传统覆盖） | 未校准基线 = "迭代 0 次" |
| **Iteration** | 一次自校准修正后的 `hModel` | 方法贡献 |

输出逐点打印、汇总表、NMSE-纹波三曲线图。**注意 Ref 组内含复合滤波器群时延对齐处理**（将 y_bf 对齐到 `filter(h_Ref,1,x)` 的线性参考），保证天花板不受时间错位污染。

### `main.m`

早期主脚本：H1 纹波固定（当前 5.8 dB）的单次实验 + 2.1.1 复增益方法（线性响应估计 `w_lin`）对照。纹波扫描功能已迁入 `amplitude_operator/amplitude_alternating_iterations.m`。

### `test/test_amplitude_alternating_iterations_mg.m`

以 MG（`get_MG`）代替 LD-PLS 作为辨识器的单次交替测试，用于验证自校准框架对辨识器的无关性。⚠️ 数据路径仍指向旧机器（`C:\Users\admin\Desktop\测试\第四次\`），运行前需改为 `D:\RF_DPD_Workspace\datasets\...`。

### `result/plot_NMSE.m`

读取手工轮次得到的纹波扫描数据并绘图（Ref / Initial / IterationFinal 三线）。

## 4. 关键参数

| 参数 | 当前值 | 含义 |
|---|---|---|
| `R` / `operator_R` | 250 | 幅度算子单边阶数（FIR 长 501，独立实参数 251） |
| `eta` / `operator_eta` | 1e-4 | 算子岭正则（无量纲，按 Gram 迹归一） |
| `iteration_count` | 1 | 外层交替次数（理论建议一次 + 验证早停） |
| `S`, `lamda` | 0.3, 1 | LD-PLS 潜变量比例与岭正则 |
| `BW` | 70 | 采集带宽（MHz） |

## 5. 主要实验结果（见 `result/result.md`）

纹波扫描（BW=70，LD-PLS 辨识）：

- **Ref**：-39 dB 附近，与纹波无关（天花板）
- **Initial**：随纹波单调恶化（-39 → -27 dB，6 dB 纹波时）
- **一次自校准**：回升至 -37 dB 左右，挽回大部分失配，但受首次估计扰动界（理论命题 3：初始误差经 `H_L·A_★` 进入算子回归残差）限制，不达天花板

## 6. 理论文档指引

`doc/幅度自校准理论_统一符号与严格论证 (2).md` 为准，包含：

- 命题 1：传统覆盖估计误差的精确分解（回归矩阵失配项）
- 命题 2：约束岭回归闭式解（式(13)）
- 命题 3：首次算子估计扰动界（式(16)，σ_min 条件）
- 命题 4：给定算子时 PA 参数误差分解（式(20)(21)）
- 命题 5：交替更新只保证反馈域目标不增（非全带 NMSE 单调）
- §14 可声称 / 不可声称清单；§15 创新点表述

## 7. 已知边界与待办

- **A(ω) ≥ 0 未验证**：当前无正值约束，估计算子局部为负时是"实偶零相位残余算子"而非严格幅度（理论 §5.1 的两选一待落地）
- **实测未接入**：`datasets/测试/2026.6.10/test2/{50–90}MHz/` 已有实测记录，尚未进入本框架
- **非线性相位**：本方法严格适用于"幅度非平坦、相位近似线性"；非线性相位是后续方向