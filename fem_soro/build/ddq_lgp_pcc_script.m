% DDQ_LGP_PCC_SCRIPT   Generate MEX-function ddq_lgp_pcc_mex from ddq_lgp_pcc.
% 
% Script generated from project 'ddq_lgp_pcc.prj' on 22-May-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'ddq_lgp_pcc'.
ARGS = cell(1,1);
ARGS{1} = cell(14,1);
ARGS{1}{1} = coder.typeof(0,[4 1]);
ARGS{1}{2} = coder.typeof(0,[4 1]);
ARGS{1}{3} = coder.typeof(0,[4 1]);
ARGS{1}{4} = coder.typeof(0,[1000   8],[1 0]);
ARGS{1}{5} = coder.typeof(0,[1000   4],[1 0]);
ARGS{1}{6} = coder.typeof(0,[4 4]);
ARGS{1}{7} = coder.typeof(0);
ARGS{1}{8} = coder.typeof(0);
ARGS{1}{9} = coder.typeof(0,[4 4]);
ARGS{1}{10} = coder.typeof(0,[4000   1],[1 0]);
ARGS{1}{11} = coder.typeof(0);
ARGS_1_12 = struct;
ARGS_1_12.mu = coder.typeof(0,[4 1]);
ARGS_1_12.L = coder.typeof(0,[4 1]);
ARGS_1_12.d = coder.typeof(0,[4 1]);
ARGS_1_12.k = coder.typeof(0,[4 1]);
ARGS_1_12.g = coder.typeof(0);
ARGS_1_12.Iz = coder.typeof(0,[4 1]);
ARGS{1}{12} = coder.typeof(ARGS_1_12);
ARGS{1}{13} = coder.typeof(0);
ARGS{1}{14} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg ddq_lgp_pcc -args ARGS{1}
