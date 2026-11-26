# test2 · 椭圆 IIR H1 幅相校准全带恢复

| 恢复配置 | H2 单边通带 | 全带 NMSE（dB） | 恢复 IQ |
|---|---:|---:|---|
| 开关 1 直接覆盖，模型侧仅 H2 | 80 MHz | −19.466444396 | [y_test2_switch1_direct_bw80.mat](../../Phase_Amplitude_Calibration/y_test2_switch1_direct_bw80.mat) |
| 开关 2 幅相补偿 | 80 MHz | −20.207420106 | [y_test2_switch2_amp_phase_bw80.mat](../../Phase_Amplitude_Calibration/y_test2_switch2_amp_phase_bw80.mat) |
| 开关 2 幅相补偿 | 44 MHz | **−37.863783831** | [y_test2_switch2_amp_phase_bw44_mu0p01_O11_all.mat](../../Phase_Amplitude_Calibration/y_test2_switch2_amp_phase_bw44_mu0p01_O11_all.mat) |

使用 test2 的实测 y0 和配对原始 IQ，采样率为 1 GS/s。H1 是 12 阶椭圆低通 IIR，通带和阻带边缘为 85/90 MHz，通带波纹为 1 dB，阻带衰减为 80 dB。开关 1 的模拟反馈为 H2(H1(y0))。开关 2 用 x 和模拟的 H1(x) 估计幅相，再将校准 FIR 与 H2 卷积。

三组恢复信号均已与 y0 对齐，并在全部 38,400 个样点上计算 NMSE。两组 80 MHz 采用原项目固定参数 μ=1、13 个潜变量、λ=1。44 MHz 组使用完整记录做幅相估计并拟合模型，固定 λ=1，扫描 1–41 个潜变量及 μ 网格 0、0.0001、0.0003、0.001、0.003、0.006、0.01、0.015、0.02、0.03、0.06、0.1、0.3、1、3。每个 MAT 文件中的变量 y 均为 38,400 × 1 的复数 IQ。

## 44 MHz 参数对照

| 配置 | μ | 潜变量数 | 全带 NMSE（dB） |
|---|---:|---:|---:|
| 原固定参数 | 1 | 13 | −36.579490200 |
| 无惩罚扫描最优 | 0 | 12 | −37.600448114 |
| **联合扫描最优** | **0.01** | **11** | **−37.863783831** |

44 MHz 的拟合、选参和评分共用同一条记录。
