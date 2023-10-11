function [x, k_hat] = preimageLE(psi, mapping)

    % From Alex Cloninger
    %  
    %Input
    % psi - new point
    % mapping - structure containing following embedding information
    %   vec - embedding coordinates Nxd
    %   val - embedding eigenvalues dx1
    %   X - original data Nxm
    %   sigma - bandwidth used in creating embedding
    %   noise - estimate of std of noise added to images (default to 0)
    %   k - number of nearest neighbors used in MDS step (default to 100)
    %   aff - for original affinity matrix K on training points, aff = sum(K,2)
    
    

    kNearest = mapping.k;
    if ~isfield(mapping,'noise')
        mapping.noise = 0;
    end
    if ~isfield(mapping,'k')
        mapping.k = 100;
    end
    if size(mapping.aff,1)>1 && size(mapping.aff,2)>1 
        mapping.aff = sum(mapping.aff,2);
    end
%Find K  
            
    A = mapping.vec*psi';
    [~,index] = sort(A,'descend');
    ind = index(1:mapping.k);
    weights = A(ind)/sum(A(ind));
    A = mapping.vec(ind,:)' * weights;
    tau = 1/norm(A) * 1.2; %slight leeway factor

    k_hat = spg_lasso(mapping.vec',psi'/norm(psi),tau);

    k_hat=k_hat*norm(psi);

    k_hat = k_hat*exp(size(mapping.X,2) * ...
        mapping.noise^2/(4*mapping.sigma^2));            
            
    %Simple K for test purposes
    % k_hat = mapping.vec * psi'/norm(psi);
    % k_hat(k_hat < .7*max(k_hat)) = 0;
    
    if sum(k_hat > (1e-3)*max(k_hat))<kNearest
        kNearest = sum(k_hat > (1e-3)*max(k_hat));
        if kNearest==0
            x=NaN(1,size(mapping.X,2));
            k_hat = NaN(size(mapping.X,1),1);
            return
        end
        display(kNearest);
    end


%Exp Dist
    Xorig = mapping.X;

    [k_val, index] = sort(k_hat,'descend');
    c_i = k_val(1:kNearest) .* ...
        sqrt(mapping.aff(index(1:kNearest)));
    e = c_i * sum(c_i);


%MDS
% Jerry here: adding a fix for the situation in which all neighbors are
% identical; in that case, just assign x to be that constant value of all
% the neighbors
d_nearest = -2*mapping.sigma^2 * log( e(1:kNearest) );

X = Xorig(index(1:kNearest),:)';
mean_x = mean(X, 2);
if isequal(X, diag(mean_x)*ones(size(X, 1), kNearest))
    x = mean_x';
else
    H = eye(kNearest) - 1/kNearest * ones(kNearest);
    [U,Sigma,V] = svd(X*H,'econ');
    totalE = sum(diag(Sigma));
    indexE = cumsum(diag(Sigma)/totalE) < .99;
    indexE(find(indexE==0,1,'first'))=1;

    U=U(:,indexE); Sigma=Sigma(indexE,indexE); V=V(:,indexE);

    Z = Sigma*V';
    z = -1/2 * diag(1./diag(Sigma)) * V' * (d_nearest - sum(Z.^2,1)');
    x_bar = 1/kNearest * sum(X,2);
    x = U*z + x_bar;
    x=x';
end

% %MDS
%     d_nearest = -2*mapping.sigma^2 * log( e(1:kNearest) );

%     X = Xorig(index(1:kNearest),:)';
%     H = eye(kNearest) - 1/kNearest * ones(kNearest);
%     [U,Sigma,V] = svd(X*H,'econ');
%     totalE = sum(diag(Sigma));
%     indexE = cumsum(diag(Sigma)/totalE) < .99;
%     indexE(find(indexE==0,1,'first'))=1;
    
%     U=U(:,indexE); Sigma=Sigma(indexE,indexE); V=V(:,indexE);

%     Z = Sigma*V';
%     z = -1/2 * diag(1./diag(Sigma)) * V' * (d_nearest - sum(Z.^2,1)');
%     x_bar = 1/kNearest * sum(X,2);
%     x = U*z + x_bar;
%     x=x';
