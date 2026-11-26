function [I, Q, complex_signal] = import_iq_binary(filename, iq_pairs_count)
% IMPORT_IQ_BINARY 导入二进制格式的交错IQ数据
%   [I, Q, complex_signal] = import_iq_binary(filename, iq_pairs_count)
%
% 输入参数:
%   filename - 二进制文件名
%   iq_pairs_count - 预期的IQ对数
%
% 输出参数:
%   I - I信号数组(int16)
%   Q - Q信号数组(int16)
%   complex_signal - 复数信号(double)

    % 计算预期字节数 (每组IQ=4字节)
    expected_bytes = iq_pairs_count * 4;

    % 获取文件信息
    file_info = dir(filename);
    if isempty(file_info)
        error('文件不存在: %s', filename);
    end

    % 检查文件大小是否匹配
    if file_info.bytes ~= expected_bytes
        warning('文件大小(%d字节)与预期(%d字节)不匹配', ...
                file_info.bytes, expected_bytes);
        iq_pairs_count = floor(file_info.bytes / 4);
    end

    % 打开文件
    fid = fopen(filename, 'rb');
    if fid == -1
        error('无法打开文件: %s', filename);
    end

    % 读取所有数据(int16格式)
    raw_data = fread(fid, iq_pairs_count * 2, 'int16');
    fclose(fid);

    % 分离I和Q通道 (交错存储: I0,Q0,I1,Q1,...)
    I = raw_data(1:2:end);  % 奇数位是I
    Q = raw_data(2:2:end);  % 偶数位是Q

    % 转换为复数信号
    %complex_signal = I + 1i*Q;
    complex_signal = double(I) + 1i*double(Q);

    % 验证数据长度
    if length(raw_data) ~= iq_pairs_count * 2
        warning('读取的数据长度(%d)与预期(%d)不符', ...
                length(raw_data), iq_pairs_count * 2);
    end
end