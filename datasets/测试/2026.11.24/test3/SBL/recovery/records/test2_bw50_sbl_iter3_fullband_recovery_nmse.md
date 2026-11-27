# test2 SBL 全带恢复 NMSE

沿用 test2 配对 IQ 和 main.m 中的 SBL 设置，K=7，M=5，迭代 3 次。H1 通带和阻带边缘为 55、60 MHz，H2 为 50、53 MHz。

用 test2 的 y0 对全部 38400 个样点计算全带 NMSE。从 MAT 文件重新加载变量 y 后复算。

| H2 单边带宽 (MHz) | 迭代次数 | 全带 NMSE (dB) | 读回复核 NMSE (dB) | 恢复信号 |
|---:|---:|---:|---:|---|
| 50 | 3 | -34.692121971 | -34.692121971 | [y_test2_bw50_sbl_iter3.mat](../y_test2_bw50_sbl_iter3.mat) |
