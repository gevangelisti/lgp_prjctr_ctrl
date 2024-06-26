% EVAL_PCC_LGP_VAR_NAT_PDP_SCRIPT   Generate MEX-function
%  eval_pcc_lgp_var_nat_pdp_mex from eval_pcc_lgp_var_nat_pdp.
% 
% Script generated from project 'eval_pcc_lgp_var_nat_pdp.prj' on 22-May-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'eval_pcc_lgp_var_nat_pdp'.
ARGS = cell(1,1);
ARGS{1} = cell(21,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0,[4 1]);
ARGS{1}{3} = coder.typeof(0,[4 1]);
ARGS{1}{4} = coder.typeof(0,[4 1]);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0,[4 4]);
ARGS{1}{7} = coder.typeof(0,[4 4]);
ARGS{1}{8} = coder.typeof(0,[4 4]);
ARGS{1}{9} = coder.typeof(0,[4 4]);
ARGS{1}{10} = coder.typeof(0,[4000 4000],[1 1]);
ARGS{1}{11} = coder.typeof(0,[1000   8],[1 0]);
ARGS{1}{12} = coder.typeof(0,[1000   4],[1 0]);
ARGS{1}{13} = coder.typeof(0,[4 4]);
ARGS{1}{14} = coder.typeof(0);
ARGS{1}{15} = coder.typeof(0);
ARGS{1}{16} = coder.typeof(0,[4 4]);
ARGS{1}{17} = coder.typeof(0,[4000   1],[1 0]);
ARGS{1}{18} = coder.typeof(0);
ARGS_1_19 = struct;
ARGS_1_19.mu = coder.typeof(0,[4 1]);
ARGS_1_19.L = coder.typeof(0,[4 1]);
ARGS_1_19.d = coder.typeof(0,[4 1]);
ARGS_1_19.k = coder.typeof(0,[4 1]);
ARGS_1_19.g = coder.typeof(0);
ARGS_1_19.Iz = coder.typeof(0,[4 1]);
ARGS{1}{19} = coder.typeof(ARGS_1_19);
ARGS{1}{20} = coder.typeof(0);
ARGS{1}{21} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg eval_pcc_lgp_var_nat_pdp -args ARGS{1}

