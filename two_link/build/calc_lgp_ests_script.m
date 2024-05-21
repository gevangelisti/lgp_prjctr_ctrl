% CALC_LGP_ESTS_SCRIPT   Generate MEX-function calc_lgp_ests_mex from
%  calc_lgp_ests.
% 
% Script generated from project 'calc_lgp_ests.prj' on 21-May-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'calc_lgp_ests'.
ARGS = cell(1,1);
ARGS{1} = cell(16,1);
ARGS{1}{1} = coder.typeof(0,[1000   4],[1 0]);
ARGS{1}{2} = coder.typeof(0,[1000   2],[1 0]);
ARGS{1}{3} = coder.typeof(0,[2 2]);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0,[2 2]);
ARGS{1}{7} = coder.typeof(0,[2 2]);
ARGS{1}{8} = coder.typeof(0,[2 2]);
ARGS{1}{9} = coder.typeof(0,[1000   4],[1 0]);
ARGS{1}{10} = coder.typeof(0,[1000   2],[1 0]);
ARGS{1}{11} = coder.typeof(0,[1000   2],[1 0]);
ARGS_1_12 = struct;
ARGS_1_12.m1h = coder.typeof(0);
ARGS_1_12.m2h = coder.typeof(0);
ARGS_1_12.l1h = coder.typeof(0);
ARGS_1_12.l2h = coder.typeof(0);
ARGS_1_12.dh = coder.typeof(0,[1 2]);
ARGS_1_12.gh = coder.typeof(0);
ARGS_1_12.alphah = coder.typeof(0);
ARGS_1_12.betah = coder.typeof(0);
ARGS_1_12.deltah = coder.typeof(0);
ARGS{1}{12} = coder.typeof(ARGS_1_12);
ARGS{1}{13} = coder.typeof(0,[2000   1],[1 0]);
ARGS{1}{14} = coder.typeof(0,[3 1]);
ARGS{1}{15} = coder.typeof(0);
ARGS{1}{16} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg calc_lgp_ests -args ARGS{1}