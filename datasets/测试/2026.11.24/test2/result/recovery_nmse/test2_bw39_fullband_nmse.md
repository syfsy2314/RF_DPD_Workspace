# test2 · ±39 MHz 全带恢复 NMSE

三种方法共用 38,400 点复数 IQ、同一带限观测和同一 MP 基函数矩阵（K=7、M=5）。采样率为 1 GS/s，H1 的通带与阻带边缘为 44/49 MHz，H2 为 39/42 MHz。LD-PLS 的岭系数固定为 λ=10。

| 方法 | 最佳配置 | 全带 NMSE（dB） | 恢复 IQ |
|---|---|---:|---|
| DOMP | 8 个基函数项 | **−37.230763501** | [y_test2_bw39_domp_O8.mat](../../DOMP/y_test2_bw39_domp_O8.mat) |
| LD-PLS，有全带惩罚 | μ=0.1，12 个潜变量 | −36.638703425 | [y_test2_bw39_ldpls_mu0p1_O12.mat](../../LD-PLS/y_test2_bw39_ldpls_mu0p1_O12.mat) |
| LD-PLS，无全带惩罚 | 9 个潜变量 | −35.869965567 | [y_test2_bw39_ldpls_mu0_O9.mat](../../LD-PLS/y_test2_bw39_ldpls_mu0_O9.mat) |

DOMP 扫描 1–42 项。两种 LD-PLS 扫描 1–41 个潜变量，有惩罚版沿用 14 个正 μ 值的网格。全部样本用于拟合和全带评分。与 40 MHz 的样本内最优值相比，DOMP、有惩罚 LD-PLS 和无惩罚 LD-PLS 的 NMSE 分别升高 0.82、1.20 和 1.16 dB。 [40 MHz 明细](test2_bw40_fullband_nmse.md)
