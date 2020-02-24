%% analysing hit rate
load hit_rate_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
hit_rate2 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_4_safetylevel_2_samplesize_334.mat
hit_rate4 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_6_safetylevel_2_samplesize_334.mat
hit_rate6 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_8_safetylevel_2_samplesize_334.mat
hit_rate8 = hit_rate
load hit_rate_datatype_1_trans_2_transpar_10_safetylevel_2_samplesize_334.mat
hit_rate10 = hit_rate
plot(hit_rate2)
hold on
plot(hit_rate4)
plot(hit_rate6)
plot(hit_rate8)
plot(hit_rate10)
%% analysing mean value
load p_c_datatype_1_trans_2_transpar_2_safetylevel_2_samplesize_334.mat
pc2 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_4_safetylevel_2_samplesize_334.mat
pc4 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_6_safetylevel_2_samplesize_334.mat
pc6 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_8_safetylevel_2_samplesize_334.mat
pc8 = p_c_save_matrix;
load p_c_datatype_1_trans_2_transpar_10_safetylevel_2_samplesize_334.mat
pc10 = p_c_save_matrix;

p_true = 8.5333e-05;

clf
plot(mean(pc2))
hold on
plot(mean(pc4))
plot(mean(pc6))
plot(mean(pc8))
plot(mean(pc10))
line([0,10], [p_true,p_true],'LineStyle','--','color','r')
ylim([0 3e-3])

%%% analysing accuracy %%
pp = 0.8;

frac2 = pc2/p_true;
frac4 = pc4/p_true;
frac6 = pc6/p_true;
frac8 = pc8/p_true;
frac10 = pc10/p_true;

clf
plot(frac2(:,8))
hold on
%plot(frac4(:,1))
plot(frac6(:,8))
%plot(frac8(:,1))
%plot(frac10(:,1))
line([0,500], [1,1],'LineStyle','--','color','r')

accuracy2 = mean(abs(frac2 - 1)<pp);
accuracy4 = mean(abs(frac4 - 1)<pp);
accuracy6 = mean(abs(frac6 - 1)<pp)
accuracy8 = mean(abs(frac8 - 1)<pp)
accuracy10 = mean(abs(frac10 - 1)<pp)
clf
plot(accuracy2)
hold on
plot(accuracy4)
plot(accuracy6)
plot(accuracy8)
plot(accuracy10)



