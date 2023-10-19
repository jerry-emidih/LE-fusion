% this is to do Laplacian Eigenmaps, as described in the original paper by Belkin and Niyogi
% uses normLap function for a lot of the work

% Jerry Emidih
% should be compatible with older version of LE preimage

function [lapX, varargout] = lapEig(X, nodims, K, varargin)
% X is data in a matrix, so m points in n dimensions (e.g. m pixels in n spatially-registered modalities).

% K is the number of intensity nearest neighbors
% sigma (out) is a scale parameter in (0, Inf)
% sigma (in) is a scale parameter in [0, Inf) but is automated if 0


% lapX is an mt by nodims graph fourier transform (nodims-bandlimited) of the mt nontrivial points in X
% lambda is a vector of the nodims increasing frequencies
% deg is the kernel density estimate for each of the mt points
% aff is the affinity (weight) matrix
% for reference, the outputs: lapX, lambda, deg, kDensity, sig, aff
% ------- Load up data
    p = inputParser;
    addRequired(p, 'X', @isnumeric)
    % addRequired(p, 'sizevec', @isnumeric)
    addRequired(p, 'nodims', @isnumeric)
    % addRequired(p, 'Ks', @isnumeric)
    addRequired(p, 'K', @isnumeric)
    addOptional(p, 'sigma', 0, @isnumeric)
    % addParameter(p, 'numWorkers', 1, @isnumeric)
    addParameter(p, 'beQuiet', false, @islogical)
    parse(p, X, nodims, K, varargin{:})
    vars = p.Results;

    X = vars.X;
    % sizevec = vars.sizevec;
    % Ks = vars.Ks;
    K = vars.K;
    sigma = vars.sigma;
    % nW = vars.numWorkers;
    nodims = vars.nodims;
    bQ = vars.beQuiet;

    nout = max(nargout,1) - 1; %nout only counts the variables after lapX
    varargout = cell(1, nout); 
    outs = struct;
    saveKernelDensity = nout > 3;
    saveWeights = nout > 4;
    

  
    % if ~bQ
    %     fprintf('\nNo spatial information will be used. \n')
    % end
    [L, sig, deg, kDensity, aff] = normLap(X, K, sigma, saveKernelDensity, saveWeights, 'beQuiet', bQ);

    [lapX, lambda] = findLapEigs(L, nodims, 1, bQ); % do first pass at eigendecomp
    [lambda, il] = sort(diag(lambda), 'ascend');
    zeroM = find((lambda < 1e-4) & (~isnan(lambda)), 1, 'last'); % assume e-val of less than 1e-9 is just 0
    if nodims == size(L, 1)
        zeroM = 0; % if doing a full eigendecomp, we want the zero eigenvector included
    end

%         if (zeroM > 1) && (zeroM < (nodims - 1))
    if (zeroM > 1) && (zeroM < nodims )
        if ~bQ
            fprintf('Found %d connected components. \n', zeroM)
        end
        [lapX, lambda] = findLapEigs(L, nodims, zeroM, bQ);
        [lambda, il] = sort(diag(lambda), 'ascend');
%         elseif zeroM >= (nodims - 1)
    elseif zeroM >= nodims
        error('ConnComp:Disconnected', 'Too many connected components (%d) in graph for requested dimension (%d).', zeroM, nodims)
    end

    lapX = lapX(:, il((zeroM + 1):end)); % cutoff first zeroM vecs and vals
    lambda = lambda((zeroM + 1):end);


    % output parsing
    outs(1).A = lambda; % second output
    outs(2).A = deg; % third, etc.
    outs(3).A = sig;
    % outs(4).A = trunc;
    if saveKernelDensity
        outs(4).A = kDensity;
        if saveWeights
            outs(5).A = aff;
        end
    end

    for iVar = 1:nout
        varargout{iVar} = outs(iVar).A;
    end
% ---

end



function [eVecs, eVals] = findLapEigs(lapMat, numOfEigs, final0, isQuiet)
    if numOfEigs == size(lapMat, 1)
        [eVecs, eVals] = eig(full(lapMat));
    else
        [eVecs, eVals] = eigs(lapMat, numOfEigs + final0, 'smallestreal', 'Display', ~isQuiet, 'SubspaceDimension',...
            min(max(2*(numOfEigs + final0), 40), 2 + max(size(lapMat, [1 2]))));
    end
    if any(isnan(eVals), 'all')
        [eVecs, eVals] = eigs(lapMat, numOfEigs + final0, 'smallestabs', 'Display', ~isQuiet);
    end
end


% old code graveyard:

%         trunc = logical(sum(X, 2));
%         truncInd = find(trunc);
% %         mt = numel(truncInd);
%         % keep track of which voxels are just 0
%         % truncInd is a vector of length mt <= m
% %         Xt = X(trunc, :);
%         Xt = X;
%         mt = size(X, 1);
%         fprintf('Finding %d nearest neighbors... \n', K)
%         tic;
%         [Ak, distA] = knnsearch(Xt, Xt, 'K', K + 1, 'NSMethod', 'exhaustive');
%         if isequal(Ak(:, 1), (1:mt)') % Check that each point is its own NN
%             Ak = Ak(:, 2:end);
%         else % Swap around to make sure each point is its own NN
%             % but before swapping, double the neighborhood size to increase
%             % the likelyhood of the point being in its own neighborhood
%             [Ak, distA] = knnsearch(Xt, Xt, 'K', min(ceil((3*mt)/2), 10*(K + 1)), 'NSMethod', 'exhaustive');
%             diffA = find(Ak(:, 1) ~= (1:mt)');
%             for i = 1:size(diffA)
%                 if ~ismember(diffA(i), Ak(diffA(i), :)) % if each point is not even in its own NNeigborhood then idk
%                     error('Distances were miscalculated.')
%                 else
%                     j = find(diffA(i) == Ak(diffA(i), :));
%                     Ak(diffA(i), 1:j) = flip(Ak(diffA(i), 1:j)); % all the distances from 1 to j should be 0 here so flipping only affects indices of points which are identical to the point we care about
%                 end
%             end
%             % return things to N + 1 size neighborhood
%             Ak = Ak(:, 1:K+1);
%             distA = distA(:, 1:K+1);
%             % then delete each point from its neighborhood
%             Ak = Ak(:, 2:end);
%         end
%         distA = distA(:, 2:end);
%         tt(2) = toc;
%         fprintf('Neighbors found in %.3f sec. \n', tt(2))
%         % sort is order big O n^2*log(n), so...
% 
% 
%     % ------- Compute graph Laplacian
%         sig = sigma;
%         if sigma == 0
%             sigma = sqrt(sum(sum(distA)) / nnz(distA));
%             sig = sigma;
%         end
%         countA = nnz(distA);
%         distA = exp(-(distA.^2)./(2*(sigma.^2))); 
%         if nnz(distA) < countA
%             warning('Bad sigma choice reduced unsymmetric edges from %d to %d.', countA, nnz(distA))
%         end
%         ix = repmat((1:mt)', 1, K);
%         W = sparse(ix(:), Ak(:), distA(:), mt, mt); % put dist into a sparse matrix        
%         
% %         W = spfun(@(x) exp(-x.^2 / 2*(sigma)^2), W);
%         kDensity = sum(W, 2); % record and normalize the kernel density estimate
%         kDensity = kDensity/(sum(kDensity));
%         W = max(W, W'); % now W is an adjacency matrix of an undirected graph
%         aff = W;
%         D = diag(sum(W, 2)); % degree matrix
%         deg = diag(D);
% %         deg = deg/(sum(deg));
%         L = diag(diag(D).^(-1/2))*(D - W)*diag(diag(D).^(-1/2)); % normalized graph Laplacian
%         L = max(L, L'); % make sure no symmetry is lost in normalization
        % We use JDQR to do the eigendecomposition -- not anymore
%         [lapX, lambda] = eigs(L, nodims + 1, 'smallestreal', 'Display', 1, 'SubspaceDimension', max(2*(nodims + 1), 40));
%         if any(isnan(lambda)) % if there are any issues with the eigendecomp, use a different one
%             [lapX, lambda] = eigs(L, nodims + 1, 'smallestabs', 'Display', 1);
%         end

% if ~isempty(sizevec)
%     m = prod(sizevec); % m is number of voxels
%     trunc = logical(sum(X, 2));
%     truncInd = find(trunc);
%     mt = numel(truncInd);
%     % keep track of which voxels are just 0
%     % truncInd is a vector of length mt <= m
%     [~, truncInv] = ismember(1:m, truncInd);
%     % truncInv maps from m-space to mt-space, truncInd returns
%     % P is the voxel lattice
%     P = zeros(m, length(sizevec));
%     [P(:,1), P(:,2), P(:,3)] = ind2sub(sizevec, 1:m);
%     if sizevec(3) == 1
%         P = P(:, 1:2); % calculations on 2d grid are more accurate sometimes, I think
%     end
% % ------- Create graphs
%     Xt = X(trunc, :);
%     Ak = zeros(mt, K); % K-nn of points determined by A
%     distA = Ak; % distances to K-nn
%     fprintf('\n Finding spatial nearest neighbors... \n')
%     tic; labs = knnsearch(P, P(trunc, :), 'K', Ks + 1, 'NSMethod', 'kdtree');
%     labs = truncInv(labs); tt(1) = toc; % get labels in truncated set
%     fprintf('Spatial neighbors found in %.3f sec. \n', tt(1))
%     if isequal(labs(:, 1), (1:mt)')
%         labs = labs(:, 2:Ks + 1);
%         fprintf('Finding nearest intensity neighbors... \n')
%     else
%         error('Spatial distances were miscalculated.')
%     end
%     tic;
%     if nW == 1 % no parallelism
%         for i = 1:mt
%             nlabs = nonzeros(labs(i, :))';
%             [tmpL, tmpD] = knnsearch(Xt(nlabs,:), Xt(i, :), 'K', K, 'NSMethod', 'exhaustive');
%             tmpL = nlabs(tmpL);
%             % first, tmpL is an index of nlabs but now it's a vector of values of nlabs which are labels of Xt
%             tmpL = padarray(tmpL, [0, max(K - numel(tmpL), 0)], 0, 'post');
%             tmpD(tmpD == 0) = 1e-12;
%             distA(i, :) = padarray(tmpD, [0, max(K - numel(tmpD), 0)], 0, 'post');
%             tmpL(tmpL == 0) = i;
%             Ak(i, :) = tmpL;
%             if i == ceil(mt/2)
%                 fprintf('Halfway done finding neighbors... \n')
%             end
%         end
%     else
%         if isempty(gcp('nocreate'))
%             parpool('local', nW);
%         else
%             fprintf('Parallel Pool is already running... \n')
%         end
%         nlabs = cell(mt, 1);
%         parfor i = 1:mt
%             nlabs{i} = nonzeros(labs(i, :))';
%         end
%         parfor i = 1:mt
%             [tmpL, tmpD] = knnsearch(Xt(nlabs{i}, :), Xt(i, :), 'K', K, 'NSMethod', 'exhaustive');
%             tmpL = nlabs{i}(tmpL);
%             % first, tmpL is an index of nlabs{i} but now its values of nlabs which are labels of Xt
%             tmpL = padarray(tmpL, [0, max(K - numel(tmpL), 0)], 0, 'post');
%             tmpD(tmpD == 0) = 1e-12;
%             distA(i, :) = padarray(tmpD, [0, max(K - numel(tmpD), 0)], 0, 'post');
%             tmpL(tmpL == 0) = i;
%             Ak(i, :) = tmpL;
%         end
%     end
%     tt(2) = toc;
%     fprintf('Neighbors found in %.3f sec. \n', tt(2))
%     if m ~= mt
%         fprintf('Some points have fewer spatial neighbors than desired. \nFinding eigendecomposition... \n')
%     else
%         fprintf('All points have the full number of nearest neighbors. \nFinding eigendecomposition... \n')
%     end
%     % this works because any 0s left will be gone after sparsifying
%     ix = repmat((1:mt)', 1, K);
%     W = sparse(ix(:), Ak(:), distA(:), mt, mt); % put dist into a sparse matrix
%     % W = max(W, W'); % symmetrize and possibly add neighbors
%     % not yet...
%     % right now only the rows of W give accurate info
% % ------- Compute graph Laplacian
%     sig = sigma;
%     if sigma == 0
%         sigma = sqrt(sum(sum(distA)) / nnz(distA));
%         sig = sigma;
%     end
%     W = spfun(@(x) exp(-x.^2 / 2*(sigma)^2), W);
%     kDensity = sum(W, 2); % record and normalize the kernel density estimate
%     kDensity = kDensity/(sum(kDensity));
%     W = max(W, W'); % now W is an adjacency matrix of an undirected graph
%     aff = W;
%     D = diag(sum(W, 2)); % degree matrix
%     deg = diag(D);
%     deg = deg/(sum(deg));
%     L = diag(diag(D).^(-1/2))*(D - W)*diag(diag(D).^(-1/2)); % normalized graph Laplacian
%     L = max(L, L'); % make sure no symmetry is lost in normalization
%     % We use JDQR to do the eigendecomposition  -- not anymore
%     [lapX, lambda] = eigs(L, nodims + 1, 'smallestreal', 'Display', 1);
%     if any(isnan(lambda)) % if there are any issues with the eigendecomp, use a different one
%         [lapX, lambda] = eigs(L, nodims + 1, 'smallestabs', 'Display', 1);
%     end
%     [lambda, il] = sort(diag(lambda));
%     lapX = lapX(:, il(2:end));
%     lambda = lambda(2:end);