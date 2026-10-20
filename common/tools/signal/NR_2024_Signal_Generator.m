function [signal,Fs] = NR_2024_Signal_Generator(bandwidth,FFT,ScSpacing)
%% Function handles
papr_red = @Clipping_andFiltering;

windowing = @txWindowOFDM;
%% OFDM waveform parameters
Carrier_BW = bandwidth;                               %Channel BW
OSF = 16;                                            %Oversampling factor
FFT_size = FFT*OSF;                                %Size of IFFT/FFT
% ScSpacing = 60*1e3;                                 %Subcarrier spacing
BW_utilization = .95;                               %Ocuppied channel BW (%)
NumSc = ceil(Carrier_BW*BW_utilization/ScSpacing)+1;%Number of active subcarriers (3168 according to 5G NR numerology)
if(mod(NumSc,2))
    NumSc=NumSc-1;
end
Fs = FFT_size*ScSpacing;                            %Sample rate of the generated signal
nOFDMSymbols = 7*2;                                 %Number of OFDM symbols to be generated
M_QAM = 64;                                         %Subcarrier modulation M-QAM

%% Symbol generation
symbol_index = randi([0 M_QAM-1], NumSc,nOFDMSymbols);  %Indexes of the M-QAM symbols to be generated
qamSymbols = qammod(symbol_index, M_QAM...
    ,'UnitAveragePower',true);                      %Perform M-QAM modulation with unitary average power.

%% Zero-padding insertion and symbol mapping.
ifft_input = [qamSymbols(1:NumSc/2,:); zeros(FFT_size-NumSc,nOFDMSymbols); qamSymbols(end-NumSc/2+1:end,:)];

%% OFDM modulation
OFDM_symbols = ifft(ifft_input);                    %OFDM modulated signal

% o = inf_ofdm_cp(:);
% [ofdm_sig_ps,o,~] = align_circular_signal(ofdm_sig_ps,o);
% a = reshape(ofdm_sig_ps(1:end),[],3);%塑形
% b = a(160+1:end,:);%裁剪
% c = fft(b,[],1);%逆变换
% d = c(:);%拉直
% scatterplot(d);%看图
%% PAPR reduction
%Per symbol PAPR reduction based on iterative clipping and filtering, that
%is, time-domain hard clipping plus frequency domain filtering of
%out-of-band emissions.
Sc_Index = [1:NumSc/2 FFT_size-NumSc/2+1:FFT_size].';   %We select the indexes of the active subcarriers.
PAPRlimit = 7;
disp_results = 0;
clipped_OFDM_signal = zeros(size(OFDM_symbols));
for k = 1:nOFDMSymbols
    signal = papr_red(OFDM_symbols(:,k),Sc_Index,PAPRlimit);
    clipped_OFDM_signal(:,k) = signal;
end
%% CP insertion
CP_duration = 144*OSF;
CP_OFDM = [clipped_OFDM_signal(end-CP_duration+1:end,:); clipped_OFDM_signal];
% o = inf_ofdm_cp(:);
% [ofdm_sig_ps,o,~] = align_circular_signal(ofdm_sig_ps,o);


%% Windowing
%Symbol-based raised-cosine windowing for spectral containment enhancement.
%this provides a smoother transition between consecutive OFDM symbols in
%time domain.
Nt = 16*OSF;
for i = 1:nOFDMSymbols
    
    [out_w, postfix] = windowing(CP_OFDM(:,i),CP_duration,Nt,'RC'); %加窗
    
    if i==1 % first generated symbol has no postfix from previous symbol!
        postfix_prev = zeros(size(postfix));
    end
    
    out_w(1:Nt) = out_w(1:Nt)+postfix_prev;
    
    if i==nOFDMSymbols % for the last symbol to be generated, append postfix directly
        signal = [signal; out_w; postfix];
    else
        signal = [signal; out_w];
    end
    
    postfix_prev = postfix;
    
end
% o = signal(10241:end-80);
% a = reshape(o,[],14);%塑形
% b = a(720+1:end,:);%裁剪
% c = fft(b,[],1);%逆变换
% d = c(:);%拉直
% scatterplot(d);%看图
end %函数结束
%% Functions

function [out,PAPR_out,Niter] = Clipping_andFiltering(in,subc,PAPRlimit)

% Peak-to-Average Power Ratio Reduction (PAPRR) for OFDM with Iterative
% Clipping and Filtering (ICF).
% INPUTS:
% - in = input time-domain OFDM symbol WITHOUT a cyclic prefix
% - subc = vector that contains the indices of the ACTIVE subcarriers
% - PAPRlimit = design target for the PAPR
%
% OUTPUTS:
% - out = output time-domain OFDM signal
% - PAPR_out = output PAPR after PAPR reduction
% - Niter = Number of iterations to reach the target PAPR

% PAPR MEASUREMENT
[Ppeak,~] = max(abs(in).^2);
Pmean = mean(abs(in).^2);
PAPR_orig = 10*log10(Ppeak/Pmean);

out = in;
PAPR_out = PAPR_orig;
Niter = 0;
while PAPR_out>PAPRlimit
    
    % ENVELOPE CLIPPING:
    [phi,A] = cart2pol(real(out),imag(out)); % cartesian-to-polar conversion
    clip_level_dB = 10*log10(mean(abs(A).^2)) + PAPRlimit - 0.1; % set clip level 0.1 dB below target --> will converge faster!
    clip_level = 10^(clip_level_dB/20);
    A(A>clip_level) = clip_level;
    in_clip = A.*exp(1i*phi);
    
    % FILTERING:
    IN_clip = fft(in_clip);
    OUT = zeros(size(IN_clip));
    OUT(subc) = IN_clip(subc); % Filter, i.e., pick only the active subcarriers and null the others
    
    out = ifft(OUT);
    
    Peak_out = max(abs(out).^2);
    P_out = mean(abs(out).^2);
    PAPR_out = 10*log10(Peak_out/P_out);
    
    Niter = Niter+1;
    
    
end

end

function [out, postfix] = txWindowOFDM(in,Ncp,Nt,varargin)

% Apply windowing to the prefix of an OFDM symbol, and calculate a
% windowed postfix. The postfix should be summed to the prefix of the
% following OFDM symbol.
%
% INPUTS:
% - IN   = Single OFDM symbol in time domain with cyclic prefix included.
% - NCP  = Cyclic prefix length
% - NT   = Window taper region length
% - WIN  = Window type. Currently only raised cosine ('RC'; default) and
%          rectangular ('rect') are supported.
%
% OUTPUTS:
% - OUT  = output OFDM symbol with windowed prefix. Same length as IN.
% - POSTFIX = windowed postfix
%

% Set default values for the inputs
win = 'RC';

if nargin>=4
    win = varargin{1};
end

in = in(:);

switch win
    case 'RC'
        rc = 1/2*(1+cos(pi*[-Nt+1/2:Nt-1/2]/Nt)); % raised cosine pulse
        win_pre = rc(1:Nt);  % prefix taper
        win_post = rc(Nt+1:2*Nt); % postfix taper
    case 'rect'
        win_pre = ones(1,Nt);
        win_post = zeros(1,Nt);
    otherwise
        rc = 1/2*(1+cos(pi*[-Nt+1/2:Nt-1/2]/Nt)); % raised cosine pulse
        win_pre = rc(1:Nt);  % prefix taper
        win_post = rc(Nt+1:2*Nt); % postfix taper
end
out = [win_pre(:).*in(1:Nt); in(Nt+1:end)];
postfix = win_post(:).*in(Ncp+1:Ncp+Nt);
end


