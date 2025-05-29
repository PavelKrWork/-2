function main_call_all_pim_functions_channel_yz

   warning off
   addpath('LTE_sig_gen_v3');
   addpath('extQuaziNewton');
   addpath('UNN');
   addpath('UNN\Current density');
   
   pkg load signal;

   file_name_prefix = 'data_files\Test_chest_NB9_1pim_fminunc_power_search_multi1_UNNLinearDipole_10\Step ';
   
   main_non_wizard_compensator_extQuaziNewton(file_name_prefix,7,'data_files\LTE_sb_2_9','data_files\LTE_sb_5_9');
   
   main_W_gen(file_name_prefix,7);

   chanl_est_first(file_name_prefix,8);

end