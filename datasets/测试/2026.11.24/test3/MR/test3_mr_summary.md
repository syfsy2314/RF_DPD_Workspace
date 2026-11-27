# test3 MR 实验总结

## 1. 全带恢复 NMSE

以 test2 的全带 PA 输出 y0 为参考，用全部 38,400 个样点评分。采样率为 1 GS/s，K=7、M=5。设 a1=1e-4、a2=0。 [参考输出](../../test2/y0.mat)

表中 BW 为 H2 的单边通带边缘。H1 的通带和阻带边缘为 BW+5、BW+10 MHz，H2 的阻带边缘为 BW+3 MHz。

| BW (MHz) | 恢复 NMSE (dB) | 恢复信号 | DPD 发送文件 |
|---:|---:|---|---|
| 50 | −35.565250 | [MAT](recovery/y_test2_bw50_mr_a1_0p0001_a2_0.mat) | [DAT](transmit/z_test2_bw50_mr_a1_0p0001_a2_0_mp75.dat) |
| 60 | −37.496754 | [MAT](recovery/y_test2_bw60_mr_a1_0p0001_a2_0.mat) | [DAT](transmit/z_test2_bw60_mr_a1_0p0001_a2_0_mp75.dat) |
| 70 | −37.773169 | [MAT](recovery/y_test2_bw70_mr_a1_0p0001_a2_0.mat) | [DAT](transmit/z_test2_bw70_mr_a1_0p0001_a2_0_mp75.dat) |
| 80 | −38.538737 | [MAT](recovery/y_test2_bw80_mr_a1_0p0001_a2_0.mat) | [DAT](transmit/z_test2_bw80_mr_a1_0p0001_a2_0_mp75.dat) |
| 91 | −38.539834 | [MAT](recovery/y_test2_bw91_mr_a1_0p0001_a2_0.mat) | [DAT](transmit/z_test2_bw91_mr_a1_0p0001_a2_0_mp75.dat) |

## 2. 实测 DPD 结果

按 get_x.m 的流程，将回采与原始输入 x 对齐，再用全部 38,400 个样点评分。主信道和邻道测量带宽均为 100 MHz，左右邻道中心偏移 ±100 MHz，平均 ACPR 在线性功率比上计算。

| BW (MHz) | EVM (%) | 左 ACPR (dBc) | 右 ACPR (dBc) | 平均 ACPR (dBc) | DPD 全带 NMSE (dB) | 回采 |
|---:|---:|---:|---:|---:|---:|---|
| 50 | 1.629 | −44.33 | −45.83 | −45.01 | −35.762 | [BIN](capture/y_dpd_50.bin) |
| 60 | 1.505 | −45.69 | −48.25 | −46.78 | −36.447 | [BIN](capture/y_dpd_60.bin) |
| 70 | 1.442 | −47.86 | −49.01 | −48.40 | −36.821 | [BIN](capture/y_dpd_70.bin) |
| 80 | 1.402 | −47.65 | −49.93 | −48.64 | −37.064 | [BIN](capture/y_dpd_80.bin) |
| 91 | 1.374 | −48.04 | −50.69 | −49.17 | −37.240 | [BIN](capture/y_dpd_90.bin) |

91 MHz 的 EVM 和 DPD NMSE 最低，分别为 1.374% 和 −37.240 dB。最低平均 ACPR 出现在 91 MHz，为 −49.17 dBc。

91 MHz 回采文件保留原名 y_dpd_90.bin，按对应的 91 MHz 发送文件归组。

[test3 跨方法总报告](../result/test3_summary_report.md)
