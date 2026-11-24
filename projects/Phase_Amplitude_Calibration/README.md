# 开关 2 相位与幅度同步校准

本项目延续 Phase_Calibration/main.m 的 PA 数据、椭圆 IIR H1、互谱频点选择、模型侧覆盖、LD-PLS 参数和评分区间。开关 2 使用同一条已知输入 x，得到 x_f=H1(x)，从同一次观测提取残余相位和归一化幅度；切回开关 1 后，用模型侧的 H1 估计与 H2 卷积恢复 PA 输出。

H1 的 85/90 MHz 通带/阻带边缘和 80 dB 阻带衰减不变，通带波纹设为 1 dB。H2 比较 80 MHz 与 44 MHz 通带边缘，与 Phase_Calibration/main.m 的当前参数一致。

## 对照

同一真实反馈 H2(H1(y0)) 下，比较仅 H2、开关 2 相位校准、开关 2 相位与幅度同步校准、已知 H1 的离线参考。所有组使用同一组 K=7、M=5、S=0.3、lamda=1 的 LD-PLS 参数，并在样本 5001–33400 评分。相位和幅度仅在原实验的可靠频段估计，在该频段边缘 5 MHz 内逐渐回到单位响应。未观测频段不声称校准准确。

从仓库根目录运行：

    matlab -batch "startup; run('projects/Phase_Amplitude_Calibration/main.m')"

在 MATLAB 桌面中直接打开幅相估计图和 PSD figure（无需重跑实验）：

    startup
    run('projects/Phase_Amplitude_Calibration/results/plot_Calibration.m')
    run('projects/Phase_Amplitude_Calibration/results/plot_PSD.m')

主实验结束时也会保留这两个可见 figure。幅相图覆盖 H1 的 ±85 MHz 通带，开关 2 观测只在可靠的约 ±49 MHz 内绘制，另显示模型侧 FIR 在完整通带的响应。PSD 图只展示 H2 44 MHz 下的真实 PA 输出与两种恢复结果，横轴为完整采样频段 ±500 MHz。结果写入 results/results.mat，并导出 calibration_response.png、nmse_comparison.png 与 psd_comparison.png；数值与限制见 results/report.md。PSD 使用相同的 Welch 参数和评分样本，每路信号分别归一化到单位 RMS，以比较频谱形状。tests/diag_model_edge.m 量化模型侧 FIR 在 H2 44 MHz 通带边缘的误差。仿真中的真实 y0 与 H1 参数仅用于生成开关 1 反馈、离线参考和评分，不用于开关 2 估计。
