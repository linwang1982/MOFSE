function [ W, H, objs, wei, weih] = DecompositionAlgorithm_wei( R, Aw, Ah, k, alpha, margin, beta1, beta2, gamma,kn,W0,H0)


    maxiter = 500; % default.


    % filter for clinical trials values
    CT = R > 0;
    % filter for unobserved associations
    UN = R == 0;
    

    Lw = zeros(size(Aw{1}));

    wei = ones(1,length(Aw))/length(Aw);
    
    Lh = zeros(size(Ah{1}));
    weih = ones(1,length(Ah))/length(Ah);
    
    pow = 2;

    for i1 = 1:length(Aw)
        Aw_knn{i1} = KNN_fun(Aw{i1},kn);
        Dw1{i1} = diag(sum(Aw_knn{i1},2));   
        Dw2{i1} = diag(sum(Aw_knn{i1},1));
        Lw1{i1} = (Dw1{i1}+Dw2{i1})-(Aw_knn{i1}+Aw_knn{i1}');
        Lw = Lw + (wei(i1)^pow)*Lw1{i1};
    end
    
    for i1 = 1:length(Ah)
        Ah_knn{i1} = KNN_fun(Ah{i1},kn);

        Dh1{i1} = diag(sum(Ah_knn{i1},2));   
        Dh2{i1} = diag(sum(Ah_knn{i1},1));

        Lh1{i1} = (Dh1{i1}+Dh2{i1})-(Ah_knn{i1}+Ah_knn{i1}');
        Lh = Lh + (weih(i1)^pow)*Lh1{i1};
    end


    Iter = 1;
    epsilon = 1e-4;
    relErr = 10e8;
    Jopt1 = 0.5*norm(CT.*(R - W0*H0), 'fro')^2 +...
            0.5*alpha*norm(UN.*(margin*ones(size(R)) - W0*H0), 'fro')^2 +...
            0.5*trace(W0'*(beta1*Lw+gamma)*W0)+...
            0.5*trace(H0*(beta2*Lh+gamma)*H0');
    objs = [Jopt1];
    while Iter <= maxiter && relErr > epsilon

        numer = (CT.*R+alpha*margin*UN)*H0';
        
        for i1 = 1:length(Aw)       
            numer = numer + beta1*(wei(i1)^pow)*(Aw_knn{i1}+Aw_knn{i1}')*W0;
        end   


        deno = (CT .* (W0*H0) + alpha*UN .* (W0*H0))*H0';
        
        for i1 = 1:length(Aw)
            deno = deno + (beta1*(wei(i1)^pow)*(Dw1{i1}+Dw2{i1})+gamma)*W0;
        end
        deno = deno + eps(numer);
        
        W = ...
            max(0,W0 .* (numer./deno));

        numer = W'*(CT.*R+alpha*margin*UN);
        for i1 = 1:length(Ah)       
            numer = numer + beta2*(weih(i1)^pow)*H0*(Ah_knn{i1}+Ah_knn{i1}');
        end
        
        deno = W'*(CT.*(W*H0) + alpha*UN.*(W*H0));
        
        for i1 = 1:length(Ah)
            deno = deno + H0*(beta2*(weih(i1)^pow)*(Dh1{i1}+Dh2{i1})+gamma);
        end
        deno = deno + eps(numer);
        
        H = ...
            max(0,H0 .* (numer./deno));
        

        deno1 = 0;
        for i1 = 1:length(Aw)
            numer1{i1} = (1/trace(W'*Lw1{i1}*W))^(1/(pow-1));
            deno1 = deno1 + numer1{i1};
        end

        Lw = zeros(size(Aw{1}));
        for i1 = 1:length(Aw)
            wei(i1) = numer1{i1}/deno1;
            Lw = Lw + (wei(i1)^pow)*Lw1{i1};
        end
        
        deno2 = 0;
        for i1 = 1:length(Ah)
            numer2{i1} = (1/trace(H*Lh1{i1}*H'))^(1/(pow-1));
            deno2 = deno2 + numer2{i1};
        end

        Lh = zeros(size(Ah{1}));
        for i1 = 1:length(Ah)
            weih(i1) = numer2{i1}/deno2;
            Lh = Lh + (weih(i1)^pow)*Lh1{i1};
        end
        
        % Compute cost function
        Jopt2 =  0.5*norm(CT.*(R - W*H), 'fro')^2 +...
            0.5*alpha*norm(UN.*(margin*ones(size(R)) - W*H), 'fro')^2 +...
            0.5*trace(W'*(beta1*Lw+gamma)*W)+...
            0.5*trace(H*(beta2*Lh+gamma)*H');

         
        relErr = Jopt1 - Jopt2;
        objs = [objs, Jopt2];
        Jopt1 = Jopt2;
        Iter = Iter + 1;

        % Remember previous iteration results
        W0 = W;
        H0 = H;   

    end

end
