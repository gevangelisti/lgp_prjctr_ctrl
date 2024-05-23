% CALC_KYY_PCC_SCRIPT   Generate MEX-function calc_Kyy_pcc_mex from
%  calc_Kyy_pcc.
% 
% Script generated from project 'calc_Kyy_pcc.prj' on 22-May-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'calc_Kyy_pcc'.
ARGS = cell(1,1);
ARGS{1} = cell(11,1);
ARGS{1}{1} = coder.typeof(0,[1000   8],[1 0]);
ARGS{1}{2} = coder.typeof(0,[1000   4],[1 0]);
ARGS{1}{3} = coder.typeof(0,[4 4]);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0,[4 4]);
ARGS{1}{7} = coder.typeof(0,[4 4]);
ARGS{1}{8} = coder.typeof(0,[4 4]);
ARGS_1_9 = struct;
ARGS_1_9.mu = coder.typeof(0,[4 1]);
ARGS_1_9.L = coder.typeof(0,[4 1]);
ARGS_1_9.d = coder.typeof(0,[4 1]);
ARGS_1_9.k = coder.typeof(0,[4 1]);
ARGS_1_9.g = coder.typeof(0);
ARGS_1_9.Iz = coder.typeof(0,[4 1]);
ARGS{1}{9} = coder.typeof(ARGS_1_9);
ARGS{1}{10} = coder.typeof(0);
ARGS{1}{11} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg calc_Kyy_pcc -args ARGS{1}

