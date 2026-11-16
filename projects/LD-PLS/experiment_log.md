# LD-PLS 实验记录

## 函数接口

固定潜变量数函数现名为 `get_LD_PLS`；仓库中的调用点同步更新。该函数仍以比例 `S` 决定 `O = round(S * K * (M + 1))`。

新增 `get_adaptive_LD_PLS(K, M, X, x, h, y_bf, lambda, T_delta, kappa_G_max)`，返回 `[y_LD_PLS, gain_rms, O]`。它沿用原始 LD-PLS 方向，不使用全带宽能量惩罚。`O` 由停止条件确定，最大值为非线性基函数列数；`T_delta=0` 关闭边际收益条件，条件数超限仍会停止。当前候选潜变量被拒绝时不会计入 `O`，第一个潜变量始终保留。

## 自适应停止准则（2026-11-16）

数据：`data/x.mat` 和 `data/y0.mat`。模拟 H1+H2 链路，带宽为 80、50、30 MHz；前 60% 用于拟合，随后 20% 选择参数，末尾 20% 仅用于测试。模型参数为 K=7、M=5、岭系数 10，最多提取 21 个潜变量。比较脚本见 `tests/verify_stop_logic_comparison.m` 和 `../Adaptive latent variable selection/tests/verify_adaptive_stop_sweep.m`。

停止时拒绝当前候选潜变量，至少保留一个。分别比较 `小增益 AND 高条件数` 与 `小增益 OR 高条件数`。原阈值 `T_delta=9e-4`、`kappa_G_max=1e3` 下，OR 停得过早。分别在验证集上选择阈值后，两种逻辑的全带宽能量惩罚模型在三个带宽上的测试 NMSE 基本相同。选择较简单的 OR 逻辑，并设 `T_delta=0` 关闭增益条件，此时只由条件数触发停止。

| 方法 | `kappa_G_max` | 80 MHz：A / NMSE | 50 MHz：A / NMSE | 30 MHz：A / NMSE |
| --- | ---: | ---: | ---: | ---: |
| 原始 LD-PLS | `1e4` | 10 / -38.98 dB | 7 / -36.75 dB | 8 / -35.60 dB |
| 全带宽能量惩罚，验证集选择固定 μ | `10^4.5` | 12 / -39.94 dB | 9 / -38.92 dB | 8 / -36.27 dB |

全带宽能量惩罚模型对应的 μ 为 0.001、0.1、0.1。上述阈值由这三组带宽的验证段选择，不能视为其他数据或带宽上的通用常数。`Adaptive latent variable selection/main.m` 使用原始 LD-PLS，因此设 `kappa_G_max=1e4`；其外部桌面数据文件在当前环境不存在，尚未直接运行该入口的 42 MHz 实验。

复现实验：从仓库根目录运行 `matlab -batch "startup; verify_adaptive_stop_sweep"`，或运行 `matlab -batch "startup; verify_stop_logic_comparison"`。

## 提取方向的自适应全带宽惩罚（2026-11-16）

在原方向惩罚公式中，用第 k 个潜变量的未惩罚方向 `w0_k` 计算 `r_k = E_full(w0_k) / E_band(w0_k)`，再设 `mu_k = alpha * r_k^p`。比较 `p=0.5` 与 `p=1`；`alpha` 从 0 和 `logspace(-8,2,21)` 中按各带宽验证 NMSE 选择。两者均使用 `kappa_G_max=10^4.5`、`T_delta=0`；数据划分、H1+H2、K、M 和岭系数同上。`alpha=0` 与原始 LD-PLS 在所有潜变量数的验证曲线上一致。脚本见 `tests/verify_adaptive_direction_penalty.m`。

| 带宽 | 固定 μ：A / 测试 NMSE | `p=0.5`：α / A / 测试 NMSE | `p=1`：α / A / 测试 NMSE |
| --- | ---: | ---: | ---: |
| 80 MHz | 0.001 / 12 / -39.94 dB | 0.000316 / 12 / -39.94 dB | 0.0001 / 12 / -39.94 dB |
| 50 MHz | 0.1 / 9 / -38.92 dB | 0.003162 / 9 / -38.99 dB | 0.001 / 8 / -39.12 dB |
| 30 MHz | 0.1 / 8 / -36.27 dB | 0.316228 / 8 / -36.28 dB | 0.003162 / 7 / -36.04 dB |

平方根变化在三个带宽上与固定 μ 基本持平，50 MHz 有约 0.07 dB 改善；线性变化在 50 MHz 改善约 0.20 dB，但 30 MHz 下降约 0.23 dB。若强制跨带宽共用一个 α 和 p，验证集选择 `alpha=0.001, p=1`，测试 NMSE 为 -39.92、-39.12、-35.69 dB，30 MHz 的下降更明显。当前自适应公式没有稳定超越固定 μ，因此保留为实验脚本，尚不替换正式方法。后两种变化方式是在观察同一测试段后继续探索的，确认收益需要新的独立数据。

复现：`matlab -batch "startup; verify_adaptive_direction_penalty"`。
