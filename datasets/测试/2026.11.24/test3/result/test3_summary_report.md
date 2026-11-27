# test3 实验总报告

各方法目录按处理阶段存放数据。recovery 保存恢复 MAT，transmit 保存 DPD 发送 DAT，capture 保存 PA 回采 BIN。各方法原有的恢复记录归档在 recovery/records。

[LD-PLS](../LD-PLS/test3_ldpls_summary.md) · [MG](../MG/test3_mg_summary.md) · [MR](../MR/test3_mr_summary.md) · [SBL](../SBL/test3_sbl_summary.md) · [全带 MP](../Fullband/test3_fullband_mp_summary.md)

## 1. 全带恢复 NMSE

合并 LD-PLS、MG、MR、SBL 的恢复记录。模型使用 test2 的配对 IQ，参照 test2/y0.mat，用全部 38,400 个样点计算全带恢复 NMSE。以下 21 组均从保存的 MAT 文件重新读取变量 y，并与原记录核对。 [全带参考输出](../../test2/y0.mat)

采样率为 1 GS/s。表中的 BW 是 H2 的单边通带边缘，H2 的阻带边缘为 BW+3 MHz。H1 的通带和阻带边缘分别为 BW+5、BW+10 MHz。各方法使用 K=7、M=5。LD-PLS 取各带宽扫描得到的全带惩罚最优参数，λ=10。MG 设 s=3，MR 设 a1=1e-4、a2=0，SBL 迭代 3 次。

### 全带恢复 NMSE (dB)

| H2 单边带宽 (MHz) | LD-PLS（有全带惩罚） | MG | MR | SBL |
|---:|---:|---:|---:|---:|
| 40 | **[−37.837](../LD-PLS/recovery/y_test2_bw40_ldpls_mu0p003_O11.mat)** | — | — | — |
| 50 | **[−37.863](../LD-PLS/recovery/y_test2_bw50_ldpls_mu0p01_O8.mat)** | [−36.740](../MG/recovery/y_test2_bw50_mg_s3.mat) | [−35.565](../MR/recovery/y_test2_bw50_mr_a1_0p0001_a2_0.mat) | [−34.692](../SBL/recovery/y_test2_bw50_sbl_iter3.mat) |
| 60 | **[−37.965](../LD-PLS/recovery/y_test2_bw60_ldpls_mu0p1_O21.mat)** | [−37.477](../MG/recovery/y_test2_bw60_mg_s3.mat) | [−37.497](../MR/recovery/y_test2_bw60_mr_a1_0p0001_a2_0.mat) | [−35.926](../SBL/recovery/y_test2_bw60_sbl_iter3.mat) |
| 70 | **[−38.408](../LD-PLS/recovery/y_test2_bw70_ldpls_mu1_O14.mat)** | [−37.761](../MG/recovery/y_test2_bw70_mg_s3.mat) | [−37.773](../MR/recovery/y_test2_bw70_mr_a1_0p0001_a2_0.mat) | [−36.196](../SBL/recovery/y_test2_bw70_sbl_iter3.mat) |
| 80 | **[−38.612](../LD-PLS/recovery/y_test2_bw80_ldpls_mu0p003_O40.mat)** | [−37.802](../MG/recovery/y_test2_bw80_mg_s3.mat) | [−38.539](../MR/recovery/y_test2_bw80_mr_a1_0p0001_a2_0.mat) | [−36.575](../SBL/recovery/y_test2_bw80_sbl_iter3.mat) |
| 91 | **[−38.643](../LD-PLS/recovery/y_test2_bw91_ldpls_mu0p01_O40.mat)** | [−37.905](../MG/recovery/y_test2_bw91_mg_s3.mat) | [−38.540](../MR/recovery/y_test2_bw91_mr_a1_0p0001_a2_0.mat) | [−36.796](../SBL/recovery/y_test2_bw91_sbl_iter3.mat) |

40 MHz 只有 LD-PLS 恢复记录。五个共同带宽下，LD-PLS 的 NMSE 均最低，其中 80、91 MHz 与 MR 分别相差 0.073、0.103 dB。表中数值链接到对应的恢复信号。

### LD-PLS 扫描最优参数

| H2 单边带宽 (MHz) | μ | 潜变量数 |
|---:|---:|---:|
| 40 | 0.003 | 11 |
| 50 | 0.01 | 8 |
| 60 | 0.1 | 21 |
| 70 | 1 | 14 |
| 80 | 0.003 | 40 |
| 91 | 0.01 | 40 |

## 2. 实测 DPD 结果

复算 test3 的 20 组 PA 输出回采，涵盖 LD-PLS（有全带惩罚）、MG（s=3）、MR（a1=1e-4、a2=0）及 SBL（3 次迭代）。全带 MP 基准取 test2 的历史 DPD 回采。 [全带 MP 参考回采](../../test2/MP/y0_dpd.bin)

按 get_x.m 的口径将每组回采与原始输入 x 做循环时延和相位对齐，用全部 38,400 个样点计算 EVM 和全带 DPD NMSE。两者共用 RMS 归一化误差。ACPR 采样率为 1 GS/s，主信道和两个邻道的测量带宽均为 100 MHz，左右邻道中心偏移为 −100 MHz 和 +100 MHz。平均 ACPR 在线性功率比上平均后转为 dBc。 [原始输入](../../../../data_38400_v0.6_test1.dat)

带宽列对应恢复模型的 H2 单边通带边缘。各方法的 y_dpd_90.bin 按对应的 91 MHz 发送文件归组。40 MHz 的 LD-PLS 有恢复信号和发送文件，尚无 DPD 回采。

### 实测结果

| H2 单边带宽 (MHz) | 方法 | EVM (%) | 左 ACPR (dBc) | 右 ACPR (dBc) | 平均 ACPR (dBc) | 全带 DPD NMSE (dB) | 回采文件 |
|---:|---|---:|---:|---:|---:|---:|---|
| 50 | LD-PLS（有全带惩罚） | 1.393 | −47.73 | −48.22 | −47.97 | −37.121 | [y_dpd_50.bin](../LD-PLS/capture/y_dpd_50.bin) |
| 50 | MG | 1.721 | −42.91 | −46.67 | −44.40 | −35.286 | [y_dpd_50.bin](../MG/capture/y_dpd_50.bin) |
| 50 | MR | 1.629 | −44.33 | −45.83 | −45.01 | −35.762 | [y_dpd_50.bin](../MR/capture/y_dpd_50.bin) |
| 50 | SBL | 2.182 | −38.64 | −40.92 | −39.63 | −33.224 | [y_dpd_50.bin](../SBL/capture/y_dpd_50.bin) |
| 60 | LD-PLS（有全带惩罚） | 1.319 | −47.15 | −49.09 | −48.02 | −37.594 | [y_dpd_60.bin](../LD-PLS/capture/y_dpd_60.bin) |
| 60 | MG | 1.555 | −45.70 | −47.39 | −46.46 | −36.168 | [y_dpd_60.bin](../MG/capture/y_dpd_60.bin) |
| 60 | MR | 1.505 | −45.69 | −48.25 | −46.78 | −36.447 | [y_dpd_60.bin](../MR/capture/y_dpd_60.bin) |
| 60 | SBL | 1.879 | −40.81 | −43.05 | −41.79 | −34.521 | [y_dpd_60.bin](../SBL/capture/y_dpd_60.bin) |
| 70 | LD-PLS（有全带惩罚） | 1.294 | −48.68 | −50.42 | −49.47 | −37.760 | [y_dpd_70.bin](../LD-PLS/capture/y_dpd_70.bin) |
| 70 | MG | 1.525 | −45.97 | −48.28 | −46.97 | −36.333 | [y_dpd_70.bin](../MG/capture/y_dpd_70.bin) |
| 70 | MR | 1.442 | −47.86 | −49.01 | −48.40 | −36.821 | [y_dpd_70.bin](../MR/capture/y_dpd_70.bin) |
| 70 | SBL | 1.795 | −42.22 | −43.79 | −42.94 | −34.920 | [y_dpd_70.bin](../SBL/capture/y_dpd_70.bin) |
| 80 | LD-PLS（有全带惩罚） | 1.293 | −47.98 | −51.12 | −49.27 | −37.768 | [y_dpd_80.bin](../LD-PLS/capture/y_dpd_80.bin) |
| 80 | MG | 1.533 | −45.38 | −48.93 | −46.80 | −36.287 | [y_dpd_80.bin](../MG/capture/y_dpd_80.bin) |
| 80 | MR | 1.402 | −47.65 | −49.93 | −48.64 | −37.064 | [y_dpd_80.bin](../MR/capture/y_dpd_80.bin) |
| 80 | SBL | 1.725 | −43.22 | −45.19 | −44.09 | −35.265 | [y_dpd_80.bin](../SBL/capture/y_dpd_80.bin) |
| 91 | LD-PLS（有全带惩罚） | 1.279 | −48.44 | −51.11 | −49.57 | −37.861 | [y_dpd_90.bin](../LD-PLS/capture/y_dpd_90.bin) |
| 91 | MG | 1.620 | −46.23 | −49.14 | −47.44 | −35.809 | [y_dpd_90.bin](../MG/capture/y_dpd_90.bin) |
| 91 | MR | 1.374 | −48.04 | −50.69 | −49.17 | −37.240 | [y_dpd_90.bin](../MR/capture/y_dpd_90.bin) |
| 91 | SBL | 1.754 | −43.70 | −45.41 | −44.47 | −35.120 | [y_dpd_90.bin](../SBL/capture/y_dpd_90.bin) |
| 全带 | MP（test2 基准） | 1.181 | −49.89 | −49.30 | −49.59 | −38.553 | [y0_dpd.bin](../../test2/MP/y0_dpd.bin) |

### 结果解读

- 同带宽比较，有全带惩罚的 LD-PLS 在五个带宽下的 EVM、平均 ACPR 和全带 DPD NMSE 均最低，MR 在这三项指标上均居第二。
- LD-PLS 的 EVM 和 NMSE 从 50 到 91 MHz 逐点改善。平均 ACPR 整体改善，70 到 80 MHz 回退 0.19 dB，91 MHz 达到 −49.57 dBc。
- MR 的 EVM、平均 ACPR 和 NMSE 均逐点改善。MG 在 70 MHz 的 EVM 和 NMSE 最低，SBL 在 80 MHz 最低，两者的平均 ACPR 均在 91 MHz 最低。
- 91 MHz 的 LD-PLS 为 test3 四种方法综合最优，EVM 为 1.279%，平均 ACPR 为 −49.57 dBc，NMSE 为 −37.861 dB。
- test2 全带 MP 基准的 NMSE 为 −38.553 dB，比 test3 的 91 MHz LD-PLS 低 0.693 dB。两组平均 ACPR 仅相差 0.014 dB。该基准来自历史回采，跨次测量可能影响对照。
