function [model] = model_gmp_domp_omp(ymed, xmed, model)

% Substracting DC
xmed = xmed(:)-mean(xmed);
ymed = ymed(:)-mean(ymed);

% Reading model configuration.
Ka = model.Ka; Kb = model.Kb; Kc = model.Kc; La = model.La; Lb = model.Lb; 
Lc = model.Lc; Mb = model.Mb; Mc = model.Mc;

% Calculation of maximum memory buffers (Qpmax) and (Qnmax)
Qa = max(La); Qb = max(Lb) + max(Mb); Qc = max(Lc) - min(Mc); 
Qpmax = max([Qa, Qb, Qc]); Qnmax = max(Mc);

N = length(ymed);

% Extraction of the signal to model
n = 1:N;
y = ymed(n(1+Qpmax:end-Qnmax));

% Building the measurement matrix
[X_mat, Rmat] = buildW(xmed, n, Qpmax, Qnmax, Ka, Kb, Kc, La, Lb, Lc, Mb, Mc);
if ~isempty(model.h)
    ymod = X_mat*model.h;
end

% If model.h is empty, perform the regression
if isempty(model.h)
    switch(model.calculo)
        case 'pinv'
            % --- Normalize columns ---
            % normX = vecnorm(X_mat);
            % Xn = X_mat * diag(1./normX);
            % 
            % h_n = (Xn' * Xn) \ (Xn' * y);
            % h_vec = diag(1./normX) * h_n;

            % Without normalization:
            h_vec = pinv(X_mat)*y;
            ymod = X_mat* h_vec;
        case 'OMP'
            [h_vec, s, nopt, h_full] = omp_domp(X_mat, y, Rmat, model.calculo);
            ymod = X_mat* h_vec;
        case 'DOMP'
            [h_vec, s, nopt, h_full] = omp_domp(X_mat, y, Rmat, model.calculo);
            ymod = X_mat* h_vec;      
    end
end

model.nmse = 20*log10(norm(ymod-y)/norm(y));

% Save results in the output structure
if exist('nopt','var')
    model.nopt = nopt;
end
if exist('s','var')
    model.s = s;
end
if exist('h_vec','var')
    model.h = h_vec;
    model.numcoef = length(h_vec);
end
if exist('h_full','var')
    model.h_full = h_full;
end

end

function [Xmat, Rmat] = buildW(x, n, Qpmax, Qnmax, Ka, Kb, Kc, La, Lb, Lc, Mb, Mc)
% This function builds the measurement matrix of the GMP model
Regr_a = []; Regr_b = []; Regr_c = [];
Xa = []; Xb = []; Xc = [];

indk = 0;
for k = 1:length(Ka)
    for l = 0:La(k)
        indk = indk + 1;
        
        Xa(:,indk) = x(n(1+Qpmax-l:end-Qnmax-l)).*(abs(x(n(1+Qpmax-l:end-Qnmax-l))).^(Ka(k)));
        if(Ka(k)~=0)
            Regr_a{indk} = sprintf('x(n-%d)*|x(n-%d)|^{%d}',l,l,Ka(k));
        else
            Regr_a{indk} = sprintf('x(n-%d)',l);
        end
    end
end


indk = 0;
for k = 1:length(Kb)
    for l = 0:Lb(k)
        for m = 1:Mb(k)
            indk = indk + 1;
            Xb(:,indk) = x(n(1+Qpmax-l:end-Qnmax-l)).*(abs(x(n(1+Qpmax-l-m:end-Qnmax-l-m))).^(Kb(k)));
            Regr_b{indk} = sprintf('x(n-%d)|x(n-%d)|^{%d}',l,l+m,Kb(k));
        end
    end
end

indk = 0;
for k = 1:length(Kc)
    for l = 0:Lc(k)
        for m = 1:Mc(k)
            indk = indk + 1;
            Xc(:,indk) = x(n(1+Qpmax-l:end-Qnmax-l)).*(abs(x(n(1+Qpmax-l+m:end-Qnmax-l+m))).^(Kc(k)));
            Regr_c{indk} = sprintf('x(n-%d)|x(n-%d)|^{%d}',l,l-m,Kc(k));
        end
    end
end

Xmat = [Xa,  Xb, Xc];
Rmat = [Regr_a,  Regr_b, Regr_c]';
end



