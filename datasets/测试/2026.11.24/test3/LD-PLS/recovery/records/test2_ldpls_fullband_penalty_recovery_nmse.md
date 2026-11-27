# test2 有全带惩罚 LD-PLS 恢复 NMSE 复核

数据来源为 test2 的配对 IQ。输入从仓库 datasets 下的 data_38400_v0.6_test1.dat 读取并做峰值归一化，参考为 test2/y0.mat。归档位置为 test3/LD-PLS。

每个带宽使用此前扫描得到的最佳 μ 和潜变量数重新拟合。每个 MAT 文件保存变量 y，为 38,400 × 1 的复数 IQ。NMSE 用 MATLAB get_NMSE 对全部 38,400 个样点计算，复核值由重新加载的 y 得出。

H2 单边通带边缘为 BW，阻带边缘为 BW+3 MHz。H1 的通带、阻带边缘分别为 BW+5、BW+10 MHz。K=7，M=5，λ=10。

| H2 单边带宽 (MHz) | μ | 潜变量数 | 原扫描 NMSE (dB) | 读回复核 NMSE (dB) | 恢复信号 |
|---:|---:|---:|---:|---:|---|
| 40 | 0.003 | 11 | -37.836630373 | -37.836630373 | [y_test2_bw40_ldpls_mu0p003_O11.mat](../y_test2_bw40_ldpls_mu0p003_O11.mat) |
| 50 | 0.01 | 8 | -37.863192723 | -37.863192723 | [y_test2_bw50_ldpls_mu0p01_O8.mat](../y_test2_bw50_ldpls_mu0p01_O8.mat) |
| 60 | 0.1 | 21 | -37.964729741 | -37.964729741 | [y_test2_bw60_ldpls_mu0p1_O21.mat](../y_test2_bw60_ldpls_mu0p1_O21.mat) |
| 70 | 1 | 14 | -38.407573731 | -38.407573731 | [y_test2_bw70_ldpls_mu1_O14.mat](../y_test2_bw70_ldpls_mu1_O14.mat) |
| 80 | 0.003 | 40 | -38.611904550 | -38.611904550 | [y_test2_bw80_ldpls_mu0p003_O40.mat](../y_test2_bw80_ldpls_mu0p003_O40.mat) |
| 91 | 0.01 | 40 | -38.642871693 | -38.642871693 | [y_test2_bw91_ldpls_mu0p01_O40.mat](../y_test2_bw91_ldpls_mu0p01_O40.mat) |

复核结果均与原扫描一致，最大绝对差为 3.54e-10 dB。
