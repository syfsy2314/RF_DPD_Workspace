# test2 MG 全带恢复 NMSE

沿用 test2 配对 IQ，K=7，M=5，s=3。H1 通带和阻带边缘为 96、101 MHz，H2 为 91、94 MHz。

用 test2 的 y0 对全部 38400 个样点计算全带 NMSE。对 MAT 文件中的变量 y 重新加载并复算。

| H2 单边带宽 (MHz) | s | 全带 NMSE (dB) | 读回复核 NMSE (dB) | 恢复信号 |
|---:|---:|---:|---:|---|
| 91 | 3 | -37.904567797 | -37.904567797 | [y_test2_bw91_mg_s3.mat](../y_test2_bw91_mg_s3.mat) |
