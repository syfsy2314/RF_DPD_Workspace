% Demo script for the Orthogonal Matching Pursuit (OMP) and the Doubly OMP
%   (DOMP) with the Generalized Memory polynomial model (GMP).
%
% Author: Juan A. Becerra.
% Date: Jan. 2019
% Version: 1.0
%
% References:
% [OMP] J. Reina-Tosina, M. Allegue-Martínez, C. Crespo-Cadenas, C. Yu, and
%   S. Cruces, ?Behavioral modeling and predistortion of power amplifiers 
%   under sparsity hypothesis,? IEEE Trans. Microw. Theory Techn., vol. 63,
%   no. 2, pp. 745?753, Feb 2015.
% [DOMP] J. A. Becerra, M. J. Madero-Ayora, J. Reina-Tosina, C. 
%   Crespo-Cadenas, J. García-Frías and G. Arce, "A Doubly Orthogonal 
%   Matching Pursuit Algorithm for Sparse Predistortion of Power 
%   Amplifiers," in IEEE Microw. Wireless Compon. Lett., vol. 28, no. 8, 
%   pp. 726-728, Aug. 2018.
% [GMP] D. Morgan, Z. Ma, J. Kim, M. Zierdt, and J. Pastalan, ?A 
%   generalized memory polynomial model for digital predistortion of RF 
%   power amplifiers,? IEEE Trans. Signal Process., vol. 54, no. 10, pp. 
%   3852? 3860, Oct 2006.

clear all; close all;

% Load a sample measurement
load demo_meas;

% Use 10% samples for identification
indices = sel_indices(x, y, 0.1);
x_id = x(indices);
y_id = y(indices);

% Model configuration (Generalized Memory Polynomial)
% Part A: 13th order (with even terms), memory fading from 15 to 1.
% Part B,C: 7th order (only odd terms). Memory in advance and delayed 1.
model.Ka = [0 1 2 3 4 5 6 7 8 9 10 11 12 ];
model.La = [15 10 5 5 5 5 5 1 1 1 1 1 1 1 ];
model.Kb = [2 4 6];
model.Lb = [1 1 1];
model.Mb = [1 1 1];
model.Kc = [2 4 6];
model.Lc = [1 1 1];
model.Mc = [1 1 1];

%% Least-squares 
% Identification
model.calculo = 'pinv'; model.h = [];
[model_LS] = model_gmp_domp_omp(y_id, x_id, model);
NMSE_LS_id = model_LS.nmse;

% Validation
[model_LS] = model_gmp_domp_omp(y, x, model_LS);
NMSE_LS_va = model_LS.nmse;


%% OMP
% Identification
model.calculo = 'OMP'; model.h = [];
[model_OMP] = model_gmp_domp_omp(y_id, x_id, model);
NMSE_OMP_id = model_OMP.nmse;

% Validation
[model_OMP] = model_gmp_domp_omp(y, x, model_OMP);
NMSE_OMP_va = model_OMP.nmse;


%% DOMP
% Identification
model.calculo = 'DOMP'; model.h = [];
[model_DOMP] = model_gmp_domp_omp(y_id, x_id, model);
NMSE_DOMP_id = model_DOMP.nmse;

% Validation
[model_DOMP] = model_gmp_domp_omp(y, x, model_DOMP);
NMSE_DOMP_va = model_DOMP.nmse;

fprintf('\n -- Results --- \n');
fprintf('LS identification (%d coeff). NMSE = %4.2f dB\n',length(model_LS.h),NMSE_LS_id);
fprintf('LS validation (%d coeff). NMSE = %4.2f  dB\n',length(model_LS.h),NMSE_LS_va);
fprintf('OMP identification (%d coeff). NMSE = %4.2f dB\n',model_OMP.nopt,NMSE_OMP_id);
fprintf('LS validation (%d coeff). NMSE = %4.2f  dB\n',model_OMP.nopt,NMSE_OMP_va);
fprintf('DOMP identification (%d coeff). NMSE = %4.2f dB\n',model_DOMP.nopt,NMSE_DOMP_id);
fprintf('DOMP validation (%d coeff). NMSE = %4.2f  dB\n',model_DOMP.nopt,NMSE_DOMP_va);