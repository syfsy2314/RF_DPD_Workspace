%泄漏地板对照实验：不同阻带衰减下的条件数-带宽曲线
%预期：若 50 MHz 以下 κ 掉头由阻带泄漏造成，则阻带越深，拐点越往低带宽移动/消失
%%
clear;
close all;
fs = 1000;
%% 导入数据
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\x_v0.6_v1.mat;
load D:\RF_DPD_Workspace\datasets\测试\2026.6.10\test1\y0_v0.6_v1.mat;
%% 扫描设置
BW_list = 100:-5:30;                 % 采集带宽扫描（MHz），从大到小
nBW = numel(BW_list);
Dstop_list = [1e-4, 1e-6, 1e-8];     % 阻带衰减对照（振幅指标）
Dpass2 = 0.0057501127785;            % H2 通带波纹（与 t.m 一致）
Dpass1 = (10^(0.5/20) - 1) / (10^(0.5/20) + 1);
Dstop1 = 0.0001;
dens = 20;
K = 7;
M = 5;
kappa_n = zeros(nBW, numel(Dstop_list));
kappa_x = zeros(nBW, numel(Dstop_list));
%% 三种阻带指标 × 带宽扫描
for iD = 1:numel(Dstop_list)
    Dstop2 = Dstop_list(iD);
    for iBW = 1:nBW
        BW = BW_list(iBW);
        % OMUX 等效 FIR（与 t.m 一致）
        [N1, Fo, Ao, W] = firpmord([BW+5, BW+10]/(fs/2), [1 0], [Dpass1, Dstop1]);
        H1_FIR = dsp.FIRFilter('Numerator', firpm(N1, Fo, Ao, W, {dens}));
        % 覆盖 FIR：阻带指标取 Dstop2（ firpmord 自动给出所需阶数）
        [N2, Fo, Ao, W] = firpmord([BW, BW+3]/(fs/2), [1 0], [Dpass2, Dstop2]);
        H2_FIR = dsp.FIRFilter('Numerator', firpm(N2, Fo, Ao, W, {dens}));
        h2 = H2_FIR.Numerator;
        % 观测链路
        y_f = H1_FIR(y0);
        [y_f, x_al, ~] = align_circular_signal(y_f, x);
        % 带限基函数与条件数
        H_x = geth_even_BLMP(x_al, K, M, h2);
        Hn = H_x ./ vecnorm(H_x);
        s1 = svd(H_x, 'econ');
        s2 = svd(Hn,  'econ');
        kappa_x(iBW, iD) = s1(1) / s1(end);
        kappa_n(iBW, iD) = s2(1) / s2(end);
        fprintf('Dstop=%.0e, BW=%3d MHz: cond(H_n) = %.3g\n', ...
            Dstop2, BW, kappa_n(iBW, iD));
    end
end
%% 绘图：κ(H_n) 三条曲线对比
figure;
plot(BW_list, log10(kappa_n), 'o-', 'LineWidth', 1.2);
xlabel('采集带宽 (MHz)');
ylabel('log_{10}(条件数)');
set(gca, 'XDir', 'reverse');
legend(arrayfun(@(d) sprintf('Dstop=%.0e', d), Dstop_list, 'UniformOutput', false), ...
    'Location', 'best');
title('阻带衰减对条件数-带宽曲线的影响（列归一化基）');
grid on;
%% 绘图：κ(H_x) 三条曲线对比
figure;
plot(BW_list, log10(kappa_x), 'o-', 'LineWidth', 1.2);
xlabel('采集带宽 (MHz)');
ylabel('log_{10}(条件数)');
set(gca, 'XDir', 'reverse');
legend(arrayfun(@(d) sprintf('Dstop=%.0e', d), Dstop_list, 'UniformOutput', false), ...
    'Location', 'best');
title('阻带衰减对条件数-带宽曲线的影响（原始基）');
grid on;
