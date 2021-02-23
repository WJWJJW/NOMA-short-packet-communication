clc; clear variables; close all;
test_idx = 1;
path_str1 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\1221_UP_test\UP_1time_' num2str(test_idx) '_u_distribution'];
path_str2 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\1221_UP_test\UP_1time_' num2str(test_idx) '_check'];
path_str3 = ['E:\WeiJie\NOMA\Matlab\Thesis_data\1221_UP_test\UP_1time_' num2str(test_idx) '_u_combination'];

load (path_str1);
load (path_str2);
load (path_str3);

user_distance

Simulated_Anealing_Pairing
User_pre_grouping

% > thred = 1
SAP_thred1_check
SAP_thred2_check
UPG_thred1_check
UPG_thred2_check