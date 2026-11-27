# test3 SBL 实验总结

## 1. 全带恢复 NMSE

以 test2 的全带 PA 输出 y0 为参考，用全部 38,400 个样点评分。采样率为 1 GS/s，K=7、M=5。迭代 3 次。 [参考输出](../../test2/y0.mat)

表中 BW 为 H2 的单边通带边缘。H1 的通带和阻带边缘为 BW+5、BW+10 MHz，H2 的阻带边缘为 BW+3 MHz。

| BW (MHz) | 恢复 NMSE (dB) | 恢复信号 | DPD 发送文件 |
|---:|---:|---|---|
| 50 | −34.692122 | [MAT](recovery/y_test2_bw50_sbl_iter3.mat) | [DAT](transmit/z_test2_bw50_sbl_iter3_mp75.dat) |
| 60 | −35.925724 | [MAT](recovery/y_test2_bw60_sbl_iter3.mat) | [DAT](transmit/z_test2_bw60_sbl_iter3_mp75.dat) |
| 70 | −36.195538 | [MAT](recovery/y_test2_bw70_sbl_iter3.mat) | [DAT](transmit/z_test2_bw70_sbl_iter3_mp75.dat) |
| 80 | −36.575387 | [MAT](recovery/y_test2_bw80_sbl_iter3.mat) | [DAT](transmit/z_test2_bw80_sbl_iter3_mp75.dat) |
| 91 | −36.795526 | [MAT](recovery/y_test2_bw91_sbl_iter3.mat) | [DAT](transmit/z_test2_bw91_sbl_iter3_mp75.dat) |

## 2. 实测 DPD 结果

按 get_x.m 的流程，将回采与原始输入 x 对齐，再用全部 38,400 个样点评分。主信道和邻道测量带宽均为 100 MHz，左右邻道中心偏移 ±100 MHz，平均 ACPR 在线性功率比上计算。

| BW (MHz) | EVM (%) | 左 ACPR (dBc) | 右 ACPR (dBc) | 平均 ACPR (dBc) | DPD 全带 NMSE (dB) | 回采 |
|---:|---:|---:|---:|---:|---:|---|
| 50 | 2.182 | −38.64 | −40.92 | −39.63 | −33.224 | [BIN](capture/y_dpd_50.bin) |
| 60 | 1.879 | −40.81 | −43.05 | −41.79 | −34.521 | [BIN](capture/y_dpd_60.bin) |
| 70 | 1.795 | −42.22 | −43.79 | −42.94 | −34.920 | [BIN](capture/y_dpd_70.bin) |
| 80 | 1.725 | −43.22 | −45.19 | −44.09 | −35.265 | [BIN](capture/y_dpd_80.bin) |
| 91 | 1.754 | −43.70 | −45.41 | −44.47 | −35.120 | [BIN](capture/y_dpd_90.bin) |

80 MHz 的 EVM 和 DPD NMSE 最低，分别为 1.725% 和 −35.265 dB。最低平均 ACPR 出现在 91 MHz，为 −44.47 dBc。

91 MHz 回采文件保留原名 y_dpd_90.bin，按对应的 91 MHz 发送文件归组。

[test3 跨方法总报告](../result/test3_summary_report.md)
