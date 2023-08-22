function [aupr_num, auc_num, rmse_num, spear_num,pcc_num] = warm_test(k, alpha, margin,beta1,beta2,gamma,y_train,DrugSimMat1, DrugfunctionS, SHHPS, kn,knew)

%% ================== Section 2: Cross-validation  ====================
% Our matrix decomposition model parameters
num_simulation = 1;
Kfold = 10;
rng(0);
tic;
auc_num = zeros(num_simulation,1);
aupr_num = zeros(num_simulation,1);
rmse_num = zeros(num_simulation,1);
pcc_num = zeros(num_simulation,1);
spear_num = zeros(num_simulation,1);
for r = 1 : num_simulation % Number of independent simulations 
    disp('===========================================================');
    disp(['============= r = ' num2str(r) ' ========== of ' num2str(num_simulation) ' =======']);
    disp('===========================================================');
    fprintf('\n');
    
    y_ori = y_train;
    y_oriPos = find(y_ori>0);
    [rowsPos,~] = find(y_ori>0);
    crossvalind('Kfold',y_oriPos,Kfold); % K-fold division
    auc_fold = zeros(Kfold,1);
    aupr_fold = zeros(Kfold,1);
    rmse_fold = zeros(Kfold,1);
    pcc_fold = zeros(Kfold,1);
    spearman_fold = zeros(Kfold,1);
    mae_fold = zeros(Kfold,1);
    
    RHAT = cell(1,10);
     for i = 1 : Kfold
        disp(['--- Run ' num2str(r) ' of ' num2str(num_simulation) ', kfold ' num2str(i) ' of ' num2str(Kfold) ' ---']);
        
        s = strcat('mask', num2str(i-1));
        mask = load('mask_mat_664.mat', s);
        mask = cell2mat(struct2cell(mask));
        a = mask(y_oriPos);
        
        test_idx = find(a==0);

        
        Test_pos_cv = y_oriPos(test_idx);
        rowsPos_cv = rowsPos(test_idx);

        Test_all_neg_cv = find(y_ori==0);
        [rowsNeg, ~] = find(y_ori==0);
        Test_res_cv = [Test_pos_cv;Test_all_neg_cv];
        rowsTest_res_cv = [rowsPos_cv;rowsNeg];

        y_train_cv = y_ori;
        y_train_cv(Test_pos_cv) = 0;
        
    
        fprintf('Training the model...\n\n');

        y_train_cv_temp = y_train_cv;
        news = find(sum(y_train_cv,1)==0);
        y_train_cv_temp(:,news) = eps;
        Ah = ones(size(y_train_cv,2)) - squareform(pdist(y_train_cv_temp','cosine'))-eye(size(y_train_cv,2));
        Ah(:,news) = 0;
        Ah(news,:) = 0;
        y_train_cv_temp = y_train_cv;

        newd = find(sum(y_train_cv,2)==0);
        y_train_cv_temp(newd,:) = eps;
        Aw = ones(size(y_train_cv,1)) - squareform(pdist(y_train_cv_temp,'cosine'))-eye(size(y_train_cv,1));
        Aw(newd,:) = 0;
        Aw(:,newd) = 0;
        
        Awt{1} = Aw;
        Awt{2} = DrugSimMat1-eye(size(y_train_cv,1));
        Awt{3} = DrugfunctionS-eye(size(y_train_cv,1));
        
        Aht{1} = Ah;
        Aht{2} = SHHPS-eye(size(y_train_cv,2));
        
       [ndrugs, nses] = size(y_train_cv);
        W0 = rand(ndrugs, k);
        W0 = W0/sqrt(trace(W0'*W0));
        H0 = rand(k, nses);
        H0 = H0/sqrt(trace(H0*H0'));

        [W, H, ~, wei, weih] = DecompositionAlgorithm_wei(y_train_cv, Awt, Aht, k, alpha, margin, beta1, beta2, gamma,kn,W0,H0);
        disp([wei,weih]);
        if knew>0
            newd = find(sum(y_train_cv,2)==0);% new drugs
            alrd = find(sum(y_train_cv,2)~=0);% known drugs
            news = find(sum(y_train_cv,1)==0);% new side effects
            alrs = find(sum(y_train_cv,1)~=0);% known side effects
            
            tempW = zeros(size(Awt{1}));
            pow = 2;
            for i1 = 1:length(Awt)
                tempW = tempW + (wei(i1)^pow)*Awt{i1};
            end
            if ~isempty(newd)
                for i1 = 1:length(newd)
                    [B,indexA] = sort(tempW(newd(i1),alrd),'descend');
                    if sum(B(1:knew))~=0
                        W(newd(i1),:) = B(1:knew)*W(alrd(indexA(1:knew)),:)/sum(B(1:knew));
                    end
                end
            end
            
            tempH = zeros(size(Aht{1}));
            pow = 2;
            for i1 = 1:length(Aht)
                tempH = tempH + (weih(i1)^pow)*Aht{i1};
            end
            if ~isempty(news)
                for i1 = 1:length(news)
                    [B,indexA] = sort(tempH(alrs,news(i1)),'descend');
                    if sum(B(1:knew))~=0
                        H(:,news(i1)) = H(:,alrs(indexA(1:knew)))*B(1:knew)/sum(B(1:knew));
                    end
                end
            end
        end
        % prediction model
        Rhat = W*H; % drug signatures x side effect signatures
        RHAT{i} = Rhat;
        
        yy=y_ori;  
        testPos = yy(Test_pos_cv);
        predPos = Rhat(Test_pos_cv);
        rmse1 = sqrt(mean((testPos-predPos).^2));

        pcc1 = corr(testPos, predPos,'type','Pearson');
        spear1 = corr(testPos, predPos,'type','Spearman');
        mae1 = mean(abs(testPos-predPos));
%         
        yy=y_ori;
        yy(yy>0) = 1;
        
        count = 0;
        auccount = 0;
        auprcount = 0;
        
        testScore = yy(Test_res_cv);
        pred = Rhat(Test_res_cv);
        for p = 1:length(unique(rowsTest_res_cv))
            testScore_cur = testScore(rowsTest_res_cv == p);
            pred_cur = pred(rowsTest_res_cv == p);
            if sum(testScore_cur) > 0
                count = count + 1;
                [auc1, aupr, ~, ~, ~, ~] = auc(testScore_cur, pred_cur);
                auccount = auccount + auc1;
                auprcount = auprcount + aupr;
            end             
        end
%         
        auc_fold(i) = auccount/count;
        aupr_fold(i) = auprcount/count;
        rmse_fold(i) = rmse1;
        
        pcc_fold(i) = pcc1;
        spearman_fold(i) = spear1;
        mae_fold(i) = mae1;
%         
        disp([rmse1,pcc1,auccount/count,auprcount/count]);
        
     end

    auc_num(r) = mean(auc_fold);
    aupr_num(r) = mean(aupr_fold);  
    rmse_num(r) = mean(rmse_fold);
    spear_num(r) = mean(spearman_fold);
    pcc_num(r) = mean(pcc_fold);
    %disp([rmse_num,pcc_num,auc_num,aupr_num]);
    disp(['RMSE=',num2str(rmse_num),'; ', 'PCC=',num2str(pcc_num),'; ', 'AUC=',num2str(auc_num),'; ', 'AUPR=',num2str(aupr_num)]);
end
end