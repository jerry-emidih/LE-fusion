function [y, varargout] = LEPre_init(Psi, mapping)

    % make a function that does the preimage of an entire set of points, all
    % with the same embedding info 
    
    % inputs: same as preimageLE, except Psi is a matrix of points, add a flag
    % for using eigenvalues (default true)
    
    % step 1: Find all the weights in a loop, use those to figure out tau
    % step 2: Use cvx to solve for the kerner vectors,
    % hopefully quickly
    % step 3: Find all the points in a loop
    
    % step 2 stuff: it might(?) be an issue that the first point needs to work
    % well, but hopefully not
    
    %% Prelim (this stuff will change in the actual function!)
    
    % Bring in an embedding called mapping, and a matrix of points called Psi
    % (points are rows)
    
    % mapping = mapC; % test
    % Psi = mapping.vec; % just testing for now
    
    
%     kNearest = mapping.k;
    if ~isfield(mapping,'noise')
        mapping.noise = 0;
    end
    if ~isfield(mapping,'k')
        mapping.k = 100;
    end
    if size(mapping.aff,1)>1 && size(mapping.aff,2)>1 
        mapping.aff = sum(mapping.aff,2);
    end
    
    
    
    M = size(Psi, 1); % number of out-of-sample points
    
%     nout = max(nargout, 1) - 1;
    convCheck = nargout > 1;
    conv = false(M, 1); % checking the convergence of out-of-sample points
    multi = M > 1;
%     N = size(mapping.vec, 1); % number of in-sample points
    IL = diag(ones(size(mapping.val)) - mapping.val);
    invIL = diag(1./(ones(size(mapping.val)) - mapping.val));

    coeffs = invIL*mapping.vec'; % the 'A' in 'Ax - b'
    cfs = mapping.vec*IL; % for weights and for x0
    
    y = zeros(M, size(mapping.X, 2)); % preallocating preimages
    Xorig = mapping.X; % needed for part 3

    opts.verbosity = 1;
    opts.iterations = 30;

    
    %% Step 1:
    tau = zeros(1, M);
    thres = zeros(1, M);
    
    for i = 1:M
    A = Psi(i, :);
    A = A';
    A = cfs*A;
    [~,index] = sort(A,'descend');
    thres(i) = A(index(floor(1.5*mapping.k))); % record a threshold for making x0
    ind = index(1:mapping.k);
    weights = A(ind)/sum(A(ind));
    A = mapping.vec(ind,:)' * weights;
    % in paper, A here is given as alpha_hat
    
    tau(i) = 1/norm(A) * 1.2; %slight leeway factor
    end
    
    
    %% Step 2 (should be in a loop w/ step 3 to avoid having a dense kernel matrix all at once)
    

%     x0 = zeros(N, 1); % SPGL1 does a multiplication to set x0 as all zeros, so here it is instead
   
    for i = 1:1
    A = Psi(i, :);
    x0 = ((cfs*A' >= thres(i)).*(cfs*A'))/norm(A);

%     damp = exp(-(pdist2(A, mapping.vec).^2)/(2*mapping.sigma^2));
%     damp = damp/max(damp);
%     x0 = damp'.*x0;

    [k_hat,~, ~, info] = spgl1(coeffs, A'/norm(A),tau(i), [], x0, opts);
    if convCheck, conv(i) = info.stat < 5; end
    k_hat = (k_hat)*norm(A);
    
    k_hat = k_hat*exp(size(mapping.X,2) * ...
    mapping.noise^2/(4*mapping.sigma^2));    
    
    near = sum(k_hat > (1e-3)*max(k_hat));
    
    kNearest = min(near, mapping.k);
    % Solve for degrees:
    
    
    [k_val, index] = sort(k_hat,'descend');
    c_i = k_val(1:kNearest) .*sqrt(mapping.aff(index(1:kNearest)));
    e = c_i * sum(c_i);
    
    %MDS
    % Jerry here: adding a fix for the situation in which all neighbors are
    % identical; in that case, just assign x to be that constant value of all
    % the neighbors
    d_nearest = -2*mapping.sigma^2 * log( e(1:kNearest) );
    
    X = Xorig(index(1:kNearest),:)';
    mean_x = mean(X, 2);
    if isequal(X, diag(mean_x)*ones(size(X, 1), kNearest))
        y(i, :) = mean_x';
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
        y(i, :)=x';
    end
    if near == 0
        y(i, :) = NaN(1, size(mapping.X, 2));
    %         k_hat = NaN(size(mapping.X, 1), 1);
    end  
    end
    
    % in the case of more than one out-of-sample point:
    if multi
    for i = 2:M
        A = Psi(i, :);
        x0 = ((cfs*A' >= thres(i)).*(cfs*A'))/norm(A);
        k_hat = spgl1(coeffs, A'/norm(A),tau(i), [], x0, opts);
        k_hat = (k_hat)*norm(A);
        if convCheck, conv(i) = info.stat < 5; end
        
        k_hat = k_hat*exp(size(mapping.X,2) * ...
        mapping.noise^2/(4*mapping.sigma^2));    
        
        near = sum(k_hat > (1e-3)*max(k_hat));
        
        kNearest = min(near, mapping.k);
        % Solve for degrees:
            
        
        [k_val, index] = sort(k_hat,'descend');
        c_i = k_val(1:kNearest) .*sqrt(mapping.aff(index(1:kNearest)));
        e = c_i * sum(c_i);
        
        %MDS
        % Jerry here: adding a fix for the situation in which all neighbors are
        % identical; in that case, just assign x to be that constant value of all
        % the neighbors
        d_nearest = -2*mapping.sigma^2 * log( e(1:kNearest) );
        
        X = Xorig(index(1:kNearest),:)';
        mean_x = mean(X, 2);
        if isequal(X, diag(mean_x)*ones(size(X, 1), kNearest))
            y(i, :) = mean_x';
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
            y(i, :)=x';
        end
        
    end
    if convCheck
        varargout(1) = {conv};
    end
    

    end
%% 