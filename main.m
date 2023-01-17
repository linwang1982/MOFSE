
%% Initialization
% Clears the variables in the environment and add folders to the path
clc; clear all; close all;


%% ================== Section 1: Load Dataset and Partition  =======================

load('frequency_664.mat');
load('DrugSimMat1.mat');
load('similarity_se.mat');

y = R;
y_train = y;

beta1 = 2; % warm start
beta2 = 2; % warm start
gamma = 0;

kn = 20; % KNN
knew = 10;

kp = [200];% Dimension of embeddings
alphap = [0.05]; % Reciprocal of variance 
marginp = [1]; % Mean
comb{1} = [];
comb{2} = [];
comb{3} = [];
comb{4} = [];
comb{5} = [];

for i = 1:1
    for j = 1:1
        for h = 1:1
            k = kp(i);
            alpha = alphap(j);
            margin = marginp(h);
            [aupr_num, auc_num, rmse_num, spear_num,pcc_num] = warm_test(k, alpha, margin,beta1,beta2,gamma,y_train,DrugSimMat1, functionS, SHHPS, kn, knew); % warm start
            comb{1} = [comb{1},k];
            comb{2} = [comb{2},alpha];
            comb{3} = [comb{3},margin];
            comb{4} = [comb{4},auc_num];
            comb{5} = [comb{5},rmse_num];
        end
    end
end
