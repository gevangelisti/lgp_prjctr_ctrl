% V_LGP_SCRIPT   Generate MEX-function V_lgp_mex from V_lgp.
% 
% Script generated from project 'V_lgp.prj' on 21-May-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'V_lgp'.
ARGS = cell(1,1);
ARGS{1} = cell(9,1);
ARGS{1}{1} = coder.typeof(0,[2 1]);
ARGS{1}{2} = coder.typeof(0,[1000   4],[1 0]);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[2 2]);
ARGS{1}{5} = coder.typeof(0,[2000   1],[1 0]);
ARGS{1}{6} = coder.typeof(0,[3 1]);
ARGS_1_7 = struct;
ARGS_1_7.m1h = coder.typeof(0);
ARGS_1_7.m2h = coder.typeof(0);
ARGS_1_7.l1h = coder.typeof(0);
ARGS_1_7.l2h = coder.typeof(0);
ARGS_1_7.dh = coder.typeof(0,[1 2]);
ARGS_1_7.gh = coder.typeof(0);
ARGS_1_7.alphah = coder.typeof(0);
ARGS_1_7.betah = coder.typeof(0);
ARGS_1_7.deltah = coder.typeof(0);
ARGS{1}{7} = coder.typeof(ARGS_1_7);
ARGS{1}{8} = coder.typeof(0);
ARGS{1}{9} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg V_lgp -args ARGS{1}