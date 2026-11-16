function results = verify_energy_penalty_across_conditions()
%VERIFY_ENERGY_PENALTY_ACROSS_CONDITIONS 固定 mu=0/1 的多条件对照。
% 三组采集数据、五组 H1+H2 带宽、两种训练长度共 30 个条件。
% 每个条件独立拟合两种方向，潜变量数由停止准则自适应决定。
% 同时检查共享 kappa 阈值，以及无惩罚版先前使用的 kappa=1e4。
% 测试段恒为最后 20%，且不参与参数或潜变量数选择。

%% 数据和固定参数
projectDir = fileparts(fileparts(mfilename('fullpath')));
workspaceDir = fileparts(fileparts(projectDir));
datasetDir = fullfile(workspaceDir, 'datasets', '测试');
dataNames = ["第四次" "第五次" "第七次"];
bandwidths = [80 65 50 40 30];  % MHz
trainFractions = [0.3 0.6];

fs = 1000;  % MHz
K = 7;
M = 5;
lambda = 10;
T_delta = 0;
kappa_G_max = 10^4.5;
kappa_G_original = 1e4;
maxO = round(0.5 * K * (M + 1));
muValues = [0 1];

nData = numel(dataNames);
nBandwidths = numel(bandwidths);
nFractions = numel(trainFractions);
nCases = nData * nBandwidths * nFractions;
dataName = strings(nCases, 1);
bandwidthMHz = zeros(nCases, 1);
trainFraction = zeros(nCases, 1);
orderNoPenalty = zeros(nCases, 1);
orderPenalty = zeros(nCases, 1);
orderOriginalThreshold = zeros(nCases, 1);
trainNoPenalty = zeros(nCases, 1);
trainPenalty = zeros(nCases, 1);
trainOriginalThreshold = zeros(nCases, 1);
testNoPenalty = zeros(nCases, 1);
testPenalty = zeros(nCases, 1);
testOriginalThreshold = zeros(nCases, 1);
deltaTest = zeros(nCases, 1);
deltaOriginalThreshold = zeros(nCases, 1);
deltaGrid = zeros(nData, nBandwidths, nFractions);
deltaOriginalGrid = zeros(nData, nBandwidths, nFractions);
caseIndex = 0;

for id = 1:nData
    %% 每组采集独立构造全带宽基函数
    input = load(fullfile(datasetDir, dataNames(id), 'x.mat'), 'x');
    output = load(fullfile(datasetDir, dataNames(id), 'y0.mat'), 'y0');
    x = input.x(:);
    y0 = output.y0(:);
    X = geth_even_1d(x, K, M);
    test = (floor(0.8 * numel(x)) + 1):numel(x);

    for ib = 1:nBandwidths
        %% H1 + H2：与 main.m 相同的测量链路
        BW = bandwidths(ib);
        h1 = design_filter([BW + 5 BW + 10], 0.5, 1e-4, fs);
        h2 = design_filter([BW BW + 3], 0.1, 1e-6, fs);
        y_f = filter(h1, 1, y0);
        [y_f, ~, ~] = align_circular_signal(y_f, x);
        y_bf = filter(h2, 1, y_f);
        H_x = filter(h2, 1, X);

        for it = 1:nFractions
            %% 同一训练段分别拟合 mu=0 和 mu=1
            train = 1:floor(trainFractions(it) * numel(x));
            hByMu = zeros(size(X, 2), numel(muValues));
            orderByMu = zeros(1, numel(muValues));
            trainNmse = zeros(1, numel(muValues));
            testNmse = zeros(1, numel(muValues));
            for im = 1:numel(muValues)
                [hByMu(:, im), orderByMu(im)] = get_energy_regularized_LD_PLS( ...
                    X(train, :), H_x(train, :), y_bf(train), lambda, ...
                    muValues(im), T_delta, kappa_G_max, maxO);
                trainNmse(im) = get_NMSE(X(train, :) * hByMu(:, im), y0(train));
                testNmse(im) = get_NMSE(X(test, :) * hByMu(:, im), y0(test));
            end
            % 第二条无惩罚基线沿用原 LD-PLS 的条件数阈值，排除阈值偏置。
            [hOriginal, originalO] = get_energy_regularized_LD_PLS( ...
                X(train, :), H_x(train, :), y_bf(train), lambda, ...
                0, T_delta, kappa_G_original, maxO);
            originalTrainNmse = get_NMSE(X(train, :) * hOriginal, y0(train));
            originalTestNmse = get_NMSE(X(test, :) * hOriginal, y0(test));

            caseIndex = caseIndex + 1;
            dataName(caseIndex) = dataNames(id);
            bandwidthMHz(caseIndex) = BW;
            trainFraction(caseIndex) = trainFractions(it);
            orderNoPenalty(caseIndex) = orderByMu(1);
            orderPenalty(caseIndex) = orderByMu(2);
            orderOriginalThreshold(caseIndex) = originalO;
            trainNoPenalty(caseIndex) = trainNmse(1);
            trainPenalty(caseIndex) = trainNmse(2);
            trainOriginalThreshold(caseIndex) = originalTrainNmse;
            testNoPenalty(caseIndex) = testNmse(1);
            testPenalty(caseIndex) = testNmse(2);
            testOriginalThreshold(caseIndex) = originalTestNmse;
            % NMSE 越低越好；正值表示添加全带宽惩罚改善了测试结果。
            deltaTest(caseIndex) = testNmse(1) - testNmse(2);
            deltaOriginalThreshold(caseIndex) = originalTestNmse - testNmse(2);
            deltaGrid(id, ib, it) = deltaTest(caseIndex);
            deltaOriginalGrid(id, ib, it) = deltaOriginalThreshold(caseIndex);
        end
    end
end

%% 逐条件结果及汇总
cases = table(dataName, bandwidthMHz, trainFraction, ...
    orderNoPenalty, orderPenalty, orderOriginalThreshold, ...
    trainNoPenalty, trainPenalty, trainOriginalThreshold, ...
    testNoPenalty, testPenalty, testOriginalThreshold, ...
    deltaTest, deltaOriginalThreshold);
disp(cases);
fprintf('共享 kappa=%g：mu=1 优于 mu=0：%d/%d 组；平均改善 %.2f dB，中位数 %.2f dB。\n', ...
    kappa_G_max, ...
    sum(deltaTest > 0), nCases, mean(deltaTest), median(deltaTest));
fprintf('原方法 kappa=%g：mu=1 优于 mu=0：%d/%d 组；平均改善 %.2f dB，中位数 %.2f dB。\n', ...
    kappa_G_original, sum(deltaOriginalThreshold > 0), nCases, ...
    mean(deltaOriginalThreshold), median(deltaOriginalThreshold));
for id = 1:nData
    selected = dataName == dataNames(id);
    fprintf('%s：共享阈值改善 %d/%d，原阈值改善 %d/%d；最差 %.2f / %.2f dB。\n', ...
        dataNames(id), sum(deltaTest(selected) > 0), sum(selected), ...
        sum(deltaOriginalThreshold(selected) > 0), sum(selected), ...
        min(deltaTest(selected)), min(deltaOriginalThreshold(selected)));
end

%% 增益热图：两行对应共享阈值与原方法阈值
figure('Name', '全带宽能量惩罚的测试 NMSE 改善');
colorLimit = max(abs([deltaTest; deltaOriginalThreshold]));
for it = 1:nFractions
    subplot(2, nFractions, it);
    imagesc(deltaGrid(:, :, it));
    clim([-colorLimit colorLimit]);
    colorScale = colorbar;
    colorScale.Label.String = '改善 (dB)';
    xticks(1:nBandwidths);
    xticklabels(string(bandwidths));
    yticks(1:nData);
    yticklabels(dataNames);
    xlabel('H2 通带带宽 (MHz)');
    ylabel('采集数据');
    title(sprintf('共享阈值，训练 %.0f%%', 100 * trainFractions(it)));
    subplot(2, nFractions, nFractions + it);
    imagesc(deltaOriginalGrid(:, :, it));
    clim([-colorLimit colorLimit]);
    colorScale = colorbar;
    colorScale.Label.String = '改善 (dB)';
    xticks(1:nBandwidths);
    xticklabels(string(bandwidths));
    yticks(1:nData);
    yticklabels(dataNames);
    xlabel('H2 通带带宽 (MHz)');
    ylabel('采集数据');
    title(sprintf('原方法阈值，训练 %.0f%%', 100 * trainFractions(it)));
end
sgtitle('测试 NMSE 差值：无惩罚 - 有惩罚；正值表示有惩罚更好');

results.cases = cases;
results.deltaGrid = deltaGrid;
results.deltaOriginalGrid = deltaOriginalGrid;
results.muValues = muValues;
results.T_delta = T_delta;
results.kappa_G_max = kappa_G_max;
results.kappa_G_original = kappa_G_original;
end

function h = design_filter(edges, rippleDb, stopRipple, fs)
%DESIGN_FILTER 由 MHz 单位的通/阻带边缘设计 FIR 低通滤波器。
passRipple = (10^(rippleDb / 20) - 1) / (10^(rippleDb / 20) + 1);
[order, frequency, amplitude, weights] = firpmord( ...
    edges / (fs / 2), [1 0], [passRipple stopRipple]);
h = firpm(order, frequency, amplitude, weights, {20});
end
