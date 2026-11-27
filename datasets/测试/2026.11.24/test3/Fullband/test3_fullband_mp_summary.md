# 全带 MP 基准

## 1. 全带恢复 NMSE

全带 MP 直接用对齐后的 test2 原始全带 PA 回采 y0 和输入 x 生成预失真信号，恢复 NMSE 不适用。发送文件用 get_z_mp 的 K=7、M=5 生成，与 test2/MP/z0.dat 逐样点一致。 [原始全带回采](../../test2/MP/y0.bin) · [发送文件](transmit/z_test2_fullband_direct_mp75.dat)

## 2. 实测 DPD 结果

对比基准采用 test2 的全带 MP 历史回采，按 test3 的统一口径对齐原始输入 x，并用全部 38,400 个样点评分。

| EVM (%) | 左 ACPR (dBc) | 右 ACPR (dBc) | 平均 ACPR (dBc) | DPD 全带 NMSE (dB) | 基准回采 |
|---:|---:|---:|---:|---:|---|
| 1.181 | −49.89 | −49.30 | −49.59 | −38.553 | [BIN](../../test2/MP/y0_dpd.bin) |

test3 的全带 MP 回采归档在本目录。 [test3 回采](capture/y_dpd.bin) · [test3 跨方法总报告](../result/test3_summary_report.md)
