function [gBB,xsupply,PARAM] = IMS2017_generate_signal(GainBB)

%   function [gBB,xsupply,PARAM] = IMS2017_generate_signal(GainBB)
%   This function takes as input:
%   - GainBB - Gain value ranging from 0 to 1 to fix the amplitude of
%   dual-band signal. The absolute value of the dual-band signal that will
%   be sent to the remoteUPCLab cannot exceed 1.
%
%   This function outputs:
%   - gBB - is a dual band signal composed by 2 OFDM signals (of 5MHz and 10MHz bandwidth) at -40MHz and 40MHz, relative to the center frequency of the PA.
%   - xsupply - is the "slow envelope" supply signal used to dynamically supply the PA. The participants MUST NOT MODIFY the way the xsupply function is calculated, 
%       i.e., they are not allowed to modify the 揝upply_Envelope_Generation.m?function. If the participants need to recalculate the supply signal from a different dual-band 
%       IQ signal than the original gBB (for example, after applying CFR techniques to the original gBB signal), they must use the aforementioned function to generate the new 
%       xsupply signal to be sent to the remoteUPCLab. Finally, the participants may optionally consider time-alignments between the supply envelope and the dual-band IQ signal sent to the remoteUPCLab. 
%   - PARAM - This is an struct containing relevant information, for example:
%      * OFDM signals bandwidth 
%      * OFDM signals separation 
%      * Length of the data vectors
%      * The sampling Frequency
%
%   Originally created by - Filipe Barradas; Date - February 2016
%   Modified by - Gabriel Montoro and Pere L. Gilabert; Date - February 2017

%% Fixed Settings
%%%%% Dual-band signal settings
PARAM.MOD.f_ofdm1 =-40e6; f_ofdm1=PARAM.MOD.f_ofdm1;
PARAM.MOD.f_ofdm2 = 40e6; f_ofdm2=PARAM.MOD.f_ofdm2;
PARAM.MOD.bw_lte1 =10e6; bw_lte1=PARAM.MOD.bw_lte1;
PARAM.MOD.bw_lte2 =10e6; bw_lte2=PARAM.MOD.bw_lte2;
GuanyBB_ofdm1=GainBB;
GuanyBB_ofdm2=GainBB;

%%%%% Clock and Length settings
PARAM.L.LBB=20000; N_samples=PARAM.L.LBB;
PARAM.CLK.LOFreq=2.14e9;
PARAM.CLK.CLKBB=50e6*4; fs=PARAM.CLK.CLKBB;
PARAM.MOD.BandPassBW=min(PARAM.CLK.CLKBB/2,abs(PARAM.MOD.f_ofdm1)+abs(PARAM.MOD.f_ofdm2)+1.5*(max(PARAM.MOD.bw_lte1,PARAM.MOD.bw_lte2))); 
%%%%% Signal capture settings
PARAM.DEMOD.LSampling=PARAM.L.LBB;
PARAM.DEMOD.SamplingRate=PARAM.CLK.CLKBB;
PARAM.DEMOD.BaseBandBW=2*(PARAM.MOD.BandPassBW/2);

T=N_samples/fs;
Nc_lte = 1;
on_slot=[1 1];
seed_ofdm = 'shuffle';
Ns = size(on_slot,1);

%% OFDM 1
signal_type = 'OFDM';

signal_opts.end_time = T;
signal_opts.sampling_frequency = fs;%采样率
signal_opts.frequency_separation = 0;%频率间隔
signal_opts.number_of_carriers = Nc_lte;%子载波数 1个
signal_opts.bandwidth = bw_lte1;
signal_opts.qam_points = 128; % 128qam
signal_opts.ofdm_carriers = 512;%% nfft点数
signal_opts.ofdm_cyclic_prefix = 512/32; %循环前缀 16
signal_opts.ofdm_guard_carriers = 512/64; %保护子载波

rng(seed_ofdm);
x_ofdm = zeros(T*fs,Ns);
for k1 = 1:Ns
    x_ofdm_slot = get_complex_envelope(signal_type,signal_opts);
    x_ofdm(:,k1) = x_ofdm_slot/sqrt(mean(abs(x_ofdm_slot).^2))*on_slot(k1,2);
end
x_ofdm1 = reshape(x_ofdm,[],1);
%% OFDM 2
signal_type = 'OFDM';

signal_opts.end_time = T;
signal_opts.sampling_frequency = fs;
signal_opts.frequency_separation = 0;
signal_opts.number_of_carriers = Nc_lte;
signal_opts.bandwidth = bw_lte2;
signal_opts.qam_points = 128;
signal_opts.ofdm_carriers = 512;
signal_opts.ofdm_cyclic_prefix = 512/32;
signal_opts.ofdm_guard_carriers = 512/64;

rng(seed_ofdm);
x_ofdm = zeros(T*fs,Ns);
for k1 = 1:Ns
    x_ofdm_slot = get_complex_envelope(signal_type,signal_opts);
    x_ofdm(:,k1) = x_ofdm_slot/sqrt(mean(abs(x_ofdm_slot).^2))*on_slot(k1,2);
end
x_ofdm2 = reshape(x_ofdm,[],1);

%% DUAL-BAND
t=[0:1:N_samples-1]';

[x_ofdm1]=dpd_normalization(x_ofdm1,'RMS');
[x_ofdm2]=dpd_normalization(x_ofdm2,'RMS');

Scaling_Factor1=6;
Scaling_Factor2=6;

x_ofdm1=dpd_FFTFilter(fs,0,bw_lte1/2,x_ofdm1,'LP'); 
x_ofdm2=dpd_FFTFilter(fs,0,bw_lte2/2,x_ofdm2,'LP'); 

x_ofdm1=x_ofdm1/Scaling_Factor1*GuanyBB_ofdm1;
x_ofdm2=x_ofdm2/Scaling_Factor2*GuanyBB_ofdm2;

omega1=2*pi*f_ofdm1/fs;
omega2=2*pi*f_ofdm2/fs;

x = x_ofdm1.*exp(1j*omega1*t)+x_ofdm2.*exp(1j*omega2*t);
gBB = x-mean(x);
%% SLOW ENVELOPE GENERATION FOR DYNAMIC SUPPLY
 [xsupply]=Supply_Envelope_Generation(gBB,PARAM);

end

function [x,varargout] = get_complex_envelope(mod,opts)
%   A function which returns a random information complex envelope for
%   modulation of a carrier signal.
%
%   The current possible modulations are:
%           -> N carrier GSM (N random information sequences) - 'GSM'
%           -> N carrier OFDM (N random information sequences) - 'OFDM'
%                   -> PSK, QPSK, 16-QAM, 64-QAM - 2,4,16,64
%
%   Created by - Filipe Barradas
%   Date - January 2013

%% Switch signal generator
switch lower(mod)
    case 'gsm'
        te = opts.end_time;
        fs = opts.sampling_frequency;
        f_sep = opts.frequency_separation;
        c_bw = opts.carrier_bandwidth;
        N = opts.number_of_carriers;
        train = opts.train;
        varargout = cell(1,N);
        [x,varargout{:}] = get_gsm(te,fs,c_bw,f_sep,N,train);
    case 'ofdm'
        te = opts.end_time;
        fs = opts.sampling_frequency;
        f_sep = opts.frequency_separation;
        N = opts.number_of_carriers;
        bw = opts.bandwidth;
        qam = opts.qam_points;
        ofdm_c = opts.ofdm_carriers;
        cp = opts.ofdm_cyclic_prefix;
        ofdm_g = opts.ofdm_guard_carriers;
        x = get_ofdm(te,fs,f_sep,N,bw,qam,ofdm_c,cp,ofdm_g);
end

end

function [x] = get_ofdm(end_time,sampling_frequency,f_sep,N,bandwidth,sub_mod,ofdm_c,cp,ofdm_g,seed)
%% Control
%Specify the sampling frequency
fs = sampling_frequency; %Hz

%Specify the final time
te = end_time; %s

%Specify the bandwidt
bw = bandwidth*(ofdm_c+1)/(ofdm_c+1-ofdm_g*2);

%% Generate information sequences
inf_seq = randi([0 1],ceil(te*bw/(ofdm_c-2*ofdm_g))*(ofdm_c-2*ofdm_g)*log2(sub_mod),N);

%% Generate modulator
%MSK modulation
%%mod = modem.qammod('M',sub_mod,'InputType','bit');
% sub_mod是调制阶数
%% Modulate
%inf_qam = modulate(mod,inf_seq); 
inf_qam = qammod(inf_seq,sub_mod,'InputType','bit');

%inf_qam = qammod(inf_seq,sub_mod);
%% 有做修改 

inf_qam = inf_qam./max(max(abs(inf_qam)));

inf_ofdm = zeros([ofdm_c+1,floor(length(inf_qam)/(ofdm_c-2*ofdm_g)),N]);
inf_ofdm([ofdm_g+1:ceil((end-1)/2),ceil((end-1)/2)+2:end-ofdm_g],:,:) = reshape(inf_qam(1:floor(length(inf_qam)/(ofdm_c-2*ofdm_g))*(ofdm_c-2*ofdm_g),:),[ofdm_c-2*ofdm_g,floor(length(inf_qam)/(ofdm_c-2*ofdm_g)),N]);
inf_ofdm = ifftshift(inf_ofdm,1);
inf_ofdm = [inf_ofdm(1:ceil((end-1)/2),:,:);zeros([(ofdm_c+1)*(round(fs/bw)-1),size(inf_ofdm,2),N]);inf_ofdm(ceil((end-1)/2)+1:end,:,:)];

inf_ofdm = ifft(inf_ofdm,[],1)*ofdm_c;%% 在这里解调是可以看到的

inf_ofdm_cp = [inf_ofdm(end-cp*round(fs/(2*bw))+1:end,:,:);inf_ofdm];%加循环序列 这一步之后就看不到了 循环前缀160个点 加了3个循环前缀

ofdm_sig = reshape(inf_ofdm_cp,[],N);%拉平序列


%% Pulse-Shpaing修改
% d_filt = fdesign.pulseshaping(round(round(fs/bw)*bw/bandwidth),'Raised Cosine','Ast,Beta',70,0.2);
% h_filt = design(d_filt);
h_filt = rcosdesign(0.2,70,round(round(fs/bw)*bw/bandwidth),"normal");%滤波器


ofdm_sig = [ofdm_sig(end-length(h_filt):end,:);ofdm_sig];
ofdm_sig_ps = filter(h_filt,1,ofdm_sig);

ofdm_sig_ps = ofdm_sig_ps(length(h_filt)+1:end,:);
%
% o = inf_ofdm_cp(:);
% [ofdm_sig_ps,o,~] = align_circular_signal(ofdm_sig_ps,o);
% a = reshape(ofdm_sig_ps(1:end),[],3);%塑形
% b = a(160+1:end,:);%裁剪
% c = fft(b,[],1);%逆变换
% d = c(:);%拉直
% scatterplot(d);%看图
%
%% Setting the frequency
ts = (0:size(ofdm_sig_ps,1)-1).'/(round(fs/bw)*bw);
t = (0:te*fs-1).'/fs;

ofdm_sig_ps = interp1(ts,ofdm_sig_ps,t,'spline');

%% Gating
gate = get_step(25e-6,length(ofdm_sig_ps)/fs,3/25e-6,fs)-get_step(length(ofdm_sig_ps)/fs-25e-6,length(ofdm_sig_ps)/fs,3/25e-6,fs);

%ofdm_sig_ps = ofdm_sig_ps.*repmat(gate,1,N);

%% Separate in frequency and add
% Generate time vector
t = (0:size(ofdm_sig_ps,1)-1)/fs;
carriers = repmat(exp(1i*2*pi*f_sep*t),N,1);
carriers = bsxfun(@power,carriers,(1:N)').';

x = sum(ofdm_sig_ps.*carriers,2).*exp(-1i*2*pi*(round(N/2)+1/2-rem(N,2)/2)*f_sep*t).';

if(size(x,1)>round(te*fs))
    x = x(1:round(te*fs));
elseif(size(x,1)<round(te*fs))
    x = [x;zeros(round(te*fs)-size(x,1),1)];
end
end

function [x,t] = get_step(t0,T,bw,fs)
    t = (0:1/fs:T-1/fs)';
    
    sigma_f = bw/3;
    sigma = 1/(2*pi*sigma_f);
    x = 1/2*(1+erf((t-t0)/(sqrt(2)*sigma)));
end

function [y]=dpd_FFTFilter(Fs,Fo,BWa,x,type)

    L=length(x);

    switch type
    case 'LP'

        BW=ceil(L*BWa/Fs);

        plantilla_fft=zeros(L,1);    
        plantilla_fft(1:1:BW,1)=ones(BW,1);
        plantilla_fft(L-BW+1:1:L,1)=ones(BW,1);

        tf=fft(x); 
        y=ifft(plantilla_fft.*tf);

        %disp(sprintf('input_filterLP= %f   output_filterLP= %f',max(abs(x)),max(abs(y)) ));

    case 'BP' 
        BW=ceil(L*BWa/Fs);
        Wo=round(L*Fo/Fs);
        W1=round(Wo-0.5*BW);

        plantilla_fft=zeros(L,1);    
        plantilla_fft(1:1:(W1+BW))=[zeros(W1,1);ones(BW,1)];
        plantilla_fft(end:-1:end+2-(W1+BW),1)=[zeros(W1-1,1);ones(BW,1)];

        tf=fft(x);
        y=ifft(plantilla_fft.*tf);


    case 'COMPLEX_BP' 
        BW=floor(L*BWa/Fs)+1;

            if (L*Fo/Fs>0) 
             Wo=ceil(L*Fo/Fs);
            else
             Wo=floor(L*Fo/Fs)+L; 
            end

        W1=round(Wo-0.5*BW); 

        plantilla_fft=zeros(L,1);
        plantilla_fft(W1+1:1:W1+BW,1)=ones(BW,1);

        tf=fft(x);
        y=ifft(plantilla_fft.*tf);
    otherwise
    end


    if (isreal(x)==1)
    y=real(y);    
    end

end

function [y]=dpd_normalization(x,type)

switch type

    case 'NONE' 
        y=x;
       
    case 'RMS'
        rms=sqrt( mean( abs(x).^2 ) );
        y=x/rms;
          
    case 'PEAK' 
        peak=max(abs(x));
        y=x/peak;
    otherwise
end
end

function [y_ofdm1,y_ofdm2]=demod_dual_BB(yBB,PARAM)
yBB1=dpd_FFTFilter(PARAM.CLK.CLKBB,PARAM.MOD.f_ofdm1,PARAM.MOD.bw_lte1*5,yBB,'COMPLEX_BP');
yBB2=dpd_FFTFilter(PARAM.CLK.CLKBB,PARAM.MOD.f_ofdm2,PARAM.MOD.bw_lte2*5,yBB,'COMPLEX_BP');
t=[0:1:PARAM.L.LBB-1]';
omega1=2*pi*PARAM.MOD.f_ofdm1/PARAM.CLK.CLKBB;
omega2=2*pi*PARAM.MOD.f_ofdm2/PARAM.CLK.CLKBB;
y_ofdm1 = yBB1.*exp(-1j*omega1*t);
y_ofdm2 = yBB2.*exp(-1j*omega2*t);
end

function [xsupply]=Supply_Envelope_Generation(gBB,PARAM)
[x_ofdm1,x_ofdm2]=demod_dual_BB(gBB,PARAM);
E=abs(x_ofdm1)+abs(x_ofdm2);
E=E/max(E);
Env_min=0.2; 
E_s=(Env_min^6+abs(E).^6).^(1/6);
xsupply=E_s/max(E_s)*1;
end


