path_of_dir = fileparts(mfilename('fullpath'));
addpath(genpath(path_of_dir));

cd(fullfile(path_of_dir,'/two_link/functions'))

calc_K_yy_script;
calc_K_y0_script;
calc_lgp_ests_script;
tau_lgp_pdp_script;
tau_lgp_nat_pdp_script;
tau_lgp_var_nat_pdp_script;
V_lgp_script;
MD_lgp_script;
eval_lgp_w_sigma_script;

cd ..
cd ..