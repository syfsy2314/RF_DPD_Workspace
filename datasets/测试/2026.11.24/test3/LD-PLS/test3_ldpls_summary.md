# test3 LD-PLS 实验总结

## 1. 全带恢复 NMSE

以 test2 的全带 PA 输出 y0 为参考，用全部 38,400 个样点评分。采样率为 1 GS/s，K=7、M=5。使用全带能量惩罚，固定 λ=10，各带宽采用该组扫描的最优 μ 和潜变量数。 [参考输出](../../test2/y0.mat)

表中 BW 为 H2 的单边通带边缘。H1 的通带和阻带边缘为 BW+5、BW+10 MHz，H2 的阻带边缘为 BW+3 MHz。

| BW (MHz) | μ | 潜变量数 | 恢复 NMSE (dB) | 恢复信号 | DPD 发送文件 |
|---:|---:|---:|---:|---|---|
| 40 | 0.003 | 11 | −37.836630 | [MAT](recovery/y_test2_bw40_ldpls_mu0p003_O11.mat) | [DAT](transmit/z_test2_bw40_ldpls_mu0p003_O11_mp75.dat) |
| 50 | 0.01 | 8 | −37.863193 | [MAT](recovery/y_test2_bw50_ldpls_mu0p01_O8.mat) | [DAT](transmit/z_test2_bw50_ldpls_mu0p01_O8_mp75.dat) |
| 60 | 0.1 | 21 | −37.964730 | [MAT](recovery/y_test2_bw60_ldpls_mu0p1_O21.mat) | [DAT](transmit/z_test2_bw60_ldpls_mu0p1_O21_mp75.dat) |
| 70 | 1 | 14 | −38.407574 | [MAT](recovery/y_test2_bw70_ldpls_mu1_O14.mat) | [DAT](transmit/z_test2_bw70_ldpls_mu1_O14_mp75.dat) |
| 80 | 0.003 | 40 | −38.611905 | [MAT](recovery/y_test2_bw80_ldpls_mu0p003_O40.mat) | [DAT](transmit/z_test2_bw80_ldpls_mu0p003_O40_mp75.dat) |
| 91 | 0.01 | 40 | −38.642872 | [MAT](recovery/y_test2_bw91_ldpls_mu0p01_O40.mat) | [DAT](transmit/z_test2_bw91_ldpls_mu0p01_O40_mp75.dat) |

## 2. 实测 DPD 结果

按 get_x.m 的流程，将回采与原始输入 x 对齐，再用全部 38,400 个样点评分。主信道和邻道测量带宽均为 100 MHz，左右邻道中心偏移 ±100 MHz，平均 ACPR 在线性功率比上计算。

| BW (MHz) | EVM (%) | 左 ACPR (dBc) | 右 ACPR (dBc) | 平均 ACPR (dBc) | DPD 全带 NMSE (dB) | 回采 |
|---:|---:|---:|---:|---:|---:|---|
| 50 | 1.393 | −47.73 | −48.22 | −47.97 | −37.121 | [BIN](capture/y_dpd_50.bin) |
| 60 | 1.319 | −47.15 | −49.09 | −48.02 | −37.594 | [BIN](capture/y_dpd_60.bin) |
| 70 | 1.294 | −48.68 | −50.42 | −49.47 | −37.760 | [BIN](capture/y_dpd_70.bin) |
| 80 | 1.293 | −47.98 | −51.12 | −49.27 | −37.768 | [BIN](capture/y_dpd_80.bin) |
| 91 | 1.279 | −48.44 | −51.11 | −49.57 | −37.861 | [BIN](capture/y_dpd_90.bin) |

91 MHz 的 EVM 和 DPD NMSE 最低，分别为 1.279% 和 −37.861 dB。最低平均 ACPR 出现在 91 MHz，为 −49.57 dBc。

91 MHz 回采文件保留原名 y_dpd_90.bin，按对应的 91 MHz 发送文件归组。

40 MHz 已完成恢复和发送文件导出，暂无 DPD 回采。

[test3 跨方法总报告](../result/test3_summary_report.md)
