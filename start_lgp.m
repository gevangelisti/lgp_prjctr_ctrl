path_of_dir = fileparts(mfilename('fullpath'));
addpath(genpath(path_of_dir));

cd(fullfile(path_of_dir,'/two_link/functions'))

if exist('calc_K_yy_mex','file') ~= 3
    calc_K_yy_script;
end
if exist('calc_K_y0_mex','file') ~= 3
    calc_K_y0_script;
end
if exist('calc_lgp_ests_mex','file') ~= 3
    calc_lgp_ests_script;
end
if exist('tau_lgp_pdp_mex','file') ~= 3
    tau_lgp_pdp_script;
end
if exist('tau_lgp_nat_pdp_mex','file') ~= 3
    tau_lgp_nat_pdp_script;
end
if exist('tau_lgp_var_nat_pdp_mex','file') ~= 3
    tau_lgp_var_nat_pdp_script;
end
if exist('V_lgp_mex','file') ~= 3
    V_lgp_script;
end
if exist('MD_lgp_mex','file') ~= 3
    MD_lgp_script;
end
if exist('eval_lgp_w_sigma_mex','file') ~= 3
    eval_lgp_w_sigma_script;
end

cd(fullfile(path_of_dir,'/fem_soro/functions'))
if exist('eval_pcc_lgp_pdp_mex','file') ~= 3
    eval_pcc_lgp_pdp_script;
end
if exist('eval_pcc_lgp_nat_pdp_mex','file') ~= 3
    eval_pcc_lgp_nat_pdp_script;
end
if exist('eval_pcc_lgp_var_nat_pdp_mex','file') ~= 3
    eval_pcc_lgp_var_nat_pdp_script;
end
if exist('tau_pcc_lgp_pdp_mex','file') ~= 3
    tau_pcc_lgp_pdp_script;
end
if exist('tau_pcc_lgp_nat_pdp_mex','file') ~= 3
    tau_pcc_lgp_nat_pdp_script;
end
if exist('tau_pcc_lgp_var_nat_pdp_mex','file') ~= 3
    tau_pcc_lgp_var_nat_pdp_script;
end
if exist('calc_Kyy_pcc_mex','file') ~= 3
    calc_Kyy_pcc_script;
end
if exist('calc_K_y0_pcc_mex','file') ~= 3
    calc_K_y0_pcc_script;
end
if exist('ddq_lgp_pcc_mex','file') ~= 3
    ddq_lgp_pcc_script;
end

cd ..
cd ..