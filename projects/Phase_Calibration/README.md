# 未知反馈滤波器相位校准

本项目研究：能否由 PA 输入与带限反馈估计 H1 残余相位，再通过 H2 覆盖和 LD-PLS 恢复 PA 全带输出。

## 实验入口

从仓库根目录运行 `startup` 后，选择所需实验：

| 入口 | 用途 | 运行命令 |
|---|---|---|
| [main.m](main.m) | 保留的七组相位范围与 H2 带宽对照 | `run('projects/Phase_Calibration/main.m')` |
| [verify_linear_backoff_phase.m](tests/verify_linear_backoff_phase.m) | 用 x_f=H1(x) 校准相位，检查去除 PA 相位后的恢复性能 | `r = verify_linear_backoff_phase();` |
| [verify_input_feedback_phase.m](tests/verify_input_feedback_phase.m) | PA BLA 残余相位及 H1 相位估计误差归因 | `r = verify_input_feedback_phase();` |

图均为独立、居中的 MATLAB figure，不保存 PNG 或 FIG。数值及报告写入对应的 `results/` 子目录。

## 相位测量与测试手段

### 测量对象

本项目测量两个信号之间随频率变化的相位差。信号定义为 `y=PA(x)`、`y_f=H1(y)`、`x_f=H1(x)`。

| 信号对 | 测量内容 | 用途与限制 |
|---|---|---|
| x,y_f | PA BLA 与 H1 的级联相位 | 实际输入与反馈估计，可能混入 PA 自身残余相位 |
| y,y_f | H1 相位 | 使用真实 PA 输出的离线参考 |
| x,x_f | 理想线性 PA 条件下的 H1 相位 | 模拟校准记录；真实低功率 PA 仍可能有线性记忆相位 |

`verify_input_feedback_phase` 用 Welch 互谱作相位误差归因。`main` 和线性区实验中的 y,y_f 离线参考，则先拟合 2048 抽头完整复数 FIR（岭参数 1e-6），再由其频率响应提取相位。真实 H1 的 `freqz` 响应仅用于最终离线核对。

### 对齐与互谱估计

先用 `align_circular_signal` 去除整数/分数采样时延与常相位差，再排除首尾各 5000 点，避免启动和循环移位边缘影响。剩余样本为 5001–33400。

对输入 x 和对齐后的反馈 u（可取 y_f 或 x_f），估计：

\[
\widehat G_{ux}(f)=\frac{S_{ux}(f)}{S_{xx}(f)},\qquad
\phi(f)=\operatorname{unwrap}(\angle\widehat G_{ux}(f)).
\]

互谱采用工程中已核对方向的 `cpsd(u,x)`，输入功率谱用 `pwelch(x)`。固定参数为 4096 点周期 Hann 窗、2048 点重叠、8192 点 FFT、中心化双边频率坐标；互谱和功率谱各做 5 频点移动平均。先平滑复数互谱，再取相位，不直接平均相位角，也不逐点计算 `FFT(u)/FFT(x)`。

### 如何限制测量频段

从同一频率网格中选取目标频段，同时排除激励不足或相干度低的频点。例如测量 ±45 MHz：

```matlab
% u 已与 xAligned 对齐；Fs 和 f 均以 MHz 为单位
Fs = 1000;
B = 45;
idx = 5001:(numel(xAligned) - 5000);
options = {hann(4096, 'periodic'), 2048, 8192, Fs, 'centered'};
[Sux, f] = cpsd(u(idx), xAligned(idx), options{:});
Sux = movmean(Sux, 5);
Sxx = movmean(pwelch(xAligned(idx), options{:}), 5);
Suu = movmean(pwelch(u(idx), options{:}), 5);
coherence = min(1, abs(Sux).^2 ./ (Sxx .* Suu));
mask = abs(f) <= B ...
    & 10 * log10(Sxx / max(Sxx)) >= -30 ...
    & coherence >= 0.8;
G = Sux ./ Sxx;
phi = unwrap(angle(G(mask)));
omega = 2 * pi * f(mask) / Fs;
A = [ones(nnz(mask), 1), omega];
phiResidual = phi - A * (A \ phi);
```

`abs(f)<=B` 指定测量范围，另外两项固定为输入功率 ≥ 峰值 −30 dB、相干度 ≥ 0.8。相干度为 `|Sux|^2/(Sxx*Suu)`。当前实测输入可靠覆盖约 ±48.83 MHz；把 B 改成 80 不会增加带外激励，不能据此声称已测准 ±80 MHz。

相位展开应沿连续频率区间进行。若可靠频点出现分离频段，应分别展开并明确各段参考，不能假设跨缺口的相位绕转已知。当前实验的可靠频段为连续区间。

### 残余相位与误差指标

展开相位还包含常相位和整体时延。上述最小二乘拟合去除 `c0+c1*omega`，得到残余非线性相位 `phiResidual`。不同方法必须使用相同频点、展开规则和拟合规则；改变拟合频段会改变最佳直线，也会改变残余相位，因此不能直接混比。

离线相位误差为：

\[
\mathrm{RMS}_{\phi}=\sqrt{\frac{1}{N_f}\sum_f
\bigl(\widehat\phi_{H1,r}(f)-\phi_{H1,r}(f)\bigr)^2}.
\]

结果同时报告 rad 和度（乘以 180/π）、实际可靠频率范围及覆盖率。初始验证使用 0.10 rad 门限；它只判断已观测频段，不能保证高精度全带恢复。PA 残余相位和互谱偏差的离线分解也必须使用同一去直线算子，两项可能增强或抵消。

### 相位测量频段与 H2 通带的区别

相位测量范围由频点筛选决定；H2 是实际反馈与模型侧覆盖所用的滤波器。两者分别设置，H2 改为 45 MHz 不会自动把相位估计范围改成 ±45 MHz。当前主程序先在输入可靠的约 ±49 MHz 内估计相位，再比较 H2 ±45 MHz 和 ±80 MHz。

相位曲线通过 4001 抽头 FIR 进入模型侧，并与 H2 卷积。局部相位在约 43.83–48.83 MHz 的边缘 5 MHz 内逐渐归零，所以“估计到 ±49 MHz”不等于在整个 ±49 MHz 内都完整修正。宽频离线参考在 ±80 MHz 内保持修正，80–85 MHz 归零。这是 H2 50 MHz 和 45 MHz 实验差异的关键。

恢复测试固定同一条真实 PA 反馈、LD-PLS 参数（K=7、M=5、S=0.3、lamda=1）和样本 5001–33400。恢复后与真实 y 对齐，再用 `get_NMSE` 评分；该指标分别把两个信号归一化到单位 RMS，因此不直接度量绝对增益误差。

### 估计相位如何用于修正

以 x,x_f 校准为例，由 `S_xfx/S_xx` 取得 H1 相位，去掉常相位和线性时延后得到 `phiResidual`。实际 x,y_f 校准也执行相同计算，但结果可能混入 PA BLA 相位。

程序在频率网格上构造单位幅度响应 `Hphase(f)=exp(1j*w(f)*phiResidual(f))`。其中 w 在可靠频段内部为 1，边缘 5 MHz 内下降至 0，频段外为 0。逆 FFT 后取中心 4001 抽头，得到带有已知实现时延的因果 FIR `hPhase`；因此实际 FIR 幅度只近似为 1。再构造模型侧算子：

```matlab
hPhase = create_phase_only_fir(fPhase, phiResidual, Fs, 8192, 4001, 5);
hModel = conv(hPhase, h2(:));
yFeedback = filter(h2, 1, yfAligned);
reference = filter(hModel, 1, xAligned);
[feedbackAligned, ~, ~] = align_circular_signal(yFeedback, reference);
X = geth_even_1d(xAligned, K, M);
yRecovered = get_OPLS(K, M, X, xAligned, hModel, feedbackAligned, S, lamda);
```

修正发生在模型侧：真实测量为 `yFeedback=H2*H1*y`；设 PA 模型输出为 `yModel=X*theta`，LD-PLS 用 `H2*Hphase*X` 作为反馈域基函数，拟合与之对齐的 `yFeedback`。拟合出 theta 后，使用未经过反馈算子的 `X*theta` 得到全带输出。程序中的 `get_OPLS` 还执行线性分量剥离和 PLS 正则拟合。

这里相位响应使用正号 `exp(+1j*phiResidual)`，因为模型要重现 H1 对 PA 输出施加的相移。常相位与实现时延通过对齐处理，增益由拟合处理。本轮没有构造 H1 的幅度逆算子，也没有把测量反馈直接逆滤波；恢复性能取决于这个相位模型是否准确，以及 PA 基函数是否足够表达真实输出。

## 当前结论与设置

- 实测记录 38,400 点，采样率 1 GHz；H1 为 85/90 MHz、0.1 dB 波纹的 13 阶椭圆 IIR。
- x 的可靠频段约为 ±48.83 MHz；H2 当前比较 ±45 MHz 和 ±80 MHz 通带。
- H2 45 MHz 下，x,y_f 相位校准的恢复 NMSE 为 −31.28 dB；理想线性区 x,x_f 校准为 −35.61 dB，接近 y,y_f 离线参考 −35.63 dB。
- 恢复输出辅助迭代仅有小幅收益，已停止；其结果和 H2 50 MHz 历史对照保留在[理论与实验摘要](results/phase_experiment_summary.md)。

## 保留的辅助函数

- [create_phase_only_fir.m](create_phase_only_fir.m)：将局部残余相位转换为纯相位 FIR。
- [estimate_full_fir_operator.m](estimate_full_fir_operator.m)：由已知输入输出拟合完整 FIR，仅用于离线参考。

`y=PA(x)`，`y_f=H1(y)`，`x_f=H1(x)`。真实 y 与 H1 参数只用于生成实验观测、离线参考和评分，不进入 x,y_f 或 x,x_f 的相位估计。所有 MHz 频段写作正负边缘，避免混淆双边总带宽。
