# test3 全带恢复 NMSE 汇总

合并 LD-PLS、MG、MR、SBL 的恢复记录。模型使用 test2 的配对 IQ，参照 test2/y0.mat，用全部 38,400 个样点计算全带恢复 NMSE。以下 21 组均从保存的 MAT 文件重新读取变量 y，并与原记录核对。 [全带参考输出](../../../test2/y0.mat)

采样率为 1 GS/s。表中的 BW 是 H2 的单边通带边缘，H2 的阻带边缘为 BW+3 MHz。H1 的通带和阻带边缘分别为 BW+5、BW+10 MHz。各方法使用 K=7、M=5。LD-PLS 取各带宽扫描得到的全带惩罚最优参数，λ=10。MG 设 s=3，MR 设 a1=1e-4、a2=0，SBL 迭代 3 次。

## 全带恢复 NMSE (dB)

| H2 单边带宽 (MHz) | LD-PLS（有全带惩罚） | MG | MR | SBL |
|---:|---:|---:|---:|---:|
| 40 | **[−37.837](../../LD-PLS/recovery/y_test2_bw40_ldpls_mu0p003_O11.mat)** | — | — | — |
| 50 | **[−37.863](../../LD-PLS/recovery/y_test2_bw50_ldpls_mu0p01_O8.mat)** | [−36.740](../../MG/recovery/y_test2_bw50_mg_s3.mat) | [−35.565](../../MR/recovery/y_test2_bw50_mr_a1_0p0001_a2_0.mat) | [−34.692](../../SBL/recovery/y_test2_bw50_sbl_iter3.mat) |
| 60 | **[−37.965](../../LD-PLS/recovery/y_test2_bw60_ldpls_mu0p1_O21.mat)** | [−37.477](../../MG/recovery/y_test2_bw60_mg_s3.mat) | [−37.497](../../MR/recovery/y_test2_bw60_mr_a1_0p0001_a2_0.mat) | [−35.926](../../SBL/recovery/y_test2_bw60_sbl_iter3.mat) |
| 70 | **[−38.408](../../LD-PLS/recovery/y_test2_bw70_ldpls_mu1_O14.mat)** | [−37.761](../../MG/recovery/y_test2_bw70_mg_s3.mat) | [−37.773](../../MR/recovery/y_test2_bw70_mr_a1_0p0001_a2_0.mat) | [−36.196](../../SBL/recovery/y_test2_bw70_sbl_iter3.mat) |
| 80 | **[−38.612](../../LD-PLS/recovery/y_test2_bw80_ldpls_mu0p003_O40.mat)** | [−37.802](../../MG/recovery/y_test2_bw80_mg_s3.mat) | [−38.539](../../MR/recovery/y_test2_bw80_mr_a1_0p0001_a2_0.mat) | [−36.575](../../SBL/recovery/y_test2_bw80_sbl_iter3.mat) |
| 91 | **[−38.643](../../LD-PLS/recovery/y_test2_bw91_ldpls_mu0p01_O40.mat)** | [−37.905](../../MG/recovery/y_test2_bw91_mg_s3.mat) | [−38.540](../../MR/recovery/y_test2_bw91_mr_a1_0p0001_a2_0.mat) | [−36.796](../../SBL/recovery/y_test2_bw91_sbl_iter3.mat) |

40 MHz 只有 LD-PLS 恢复记录。五个共同带宽下，LD-PLS 的 NMSE 均最低，其中 80、91 MHz 与 MR 分别相差 0.073、0.103 dB。表中数值链接到对应的恢复信号。

## LD-PLS 扫描最优参数

| H2 单边带宽 (MHz) | μ | 潜变量数 |
|---:|---:|---:|
| 40 | 0.003 | 11 |
| 50 | 0.01 | 8 |
| 60 | 0.1 | 21 |
| 70 | 1 | 14 |
| 80 | 0.003 | 40 |
| 91 | 0.01 | 40 |

恢复 NMSE 的参考是 PA 原始全带输出 y0。实测 DPD 的 NMSE 以原始输入 x 为参考，结果见 [test3 实测 DPD 报告](test3_dpd_report.md)。
