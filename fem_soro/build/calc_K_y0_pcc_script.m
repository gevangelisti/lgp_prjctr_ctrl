% CALC_K_Y0_PCC_SCRIPT   Generate MEX-function calc_K_y0_pcc_mex from
%  calc_K_y0_pcc.
% 
% Script generated from project 'calc_K_y0_pcc.prj' on 22-May-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'calc_K_y0_pcc'.
ARGS = cell(1,1);
ARGS{1} = cell(5,1);
ARGS{1}{1} = coder.typeof(0,[1000   8],[1 0]);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0,[4 4]);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg calc_K_y0_pcc -args ARGS{1}

