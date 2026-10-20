function [x_a,y_a,del_t] = align_circular_signal(x,y)
% [x_a,y_a] = align_circular_signal(x,y)
%
% Aligns the signals x and y considering that they are circular (periodic).
% The procedure is iterative, default 10 iterations.


% Author: Filipe Barradas
% Date: November 2013

%% Coarse allignment
Hf = fft(x).*conj(fft(y));
H = ifft(Hf);
[~,ind_m] = max(abs(H));

x_a = circshift(x,-(ind_m-1));
h = (x_a'*y)./(x_a'*x_a);
h = exp(1j*angle(h));
x_a = h.*x_a;

%% Fine Alignment
N = length(x_a);
f = fftshift((-1/2+1/N:1/N:1/2).');

tol = 1e-5*norm(x)*norm(y);

for k1 = 1:10
    d = delay_estimator(x_a,y,tol);
    del = exp(-1j*2*pi*f*d);
    x_a = ifft(fft(x_a).*conj(del));
    del_t = ind_m + d;
end
h = (x_a'*y)./(x_a'*x_a);
h = exp(1j*angle(h));
x_a = h.*x_a;
y_a = y;
end