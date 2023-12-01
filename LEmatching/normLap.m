% This takes in a data matrix and makes a normalized Laplacian matrix 


% the inputs are:
% Xt (data matrix, points as rows), Kt (number of neighbors), sg (sigma), getAff, getKDen

% outputs are Lap (Laplacian matrix), dg (degree), wt (weights, or affinity, or adjacency matrix),
% kDen, sig

function varargout = normLap(Xt, Kt, varargin)
    p = inputParser;
    addRequired(p, 'Xt', @isnumeric)
    addRequired(p, 'Kt', @isnumeric)
    addOptional(p, 'sg', 0, @isnumeric)
    addOptional(p, 'getKDen', false, @islogical)
    addOptional(p, 'getAff', false, @islogical)
    addParameter(p, 'beQuiet', false, @islogical)
    parse(p, Xt, Kt, varargin{:})
    vars = p.Results;

    Xt = vars.Xt;
    Kt = vars.Kt;
    sg = vars.sg;
    getKDen = vars.getKDen;
    getAff = vars.getAff;
    bQ = vars.beQuiet;

    % output order: [Lap, sig, dg, kDen, wt]
    nout = max(nargout, 1);
    varargout = cell(1, nout);
    outs = struct;


    % tnc = logical(sum(Xt, 2));
    mt = size(Xt, 1);
    if ~bQ
        fprintf('Finding %d nearest neighbors... \n', Kt)
    end
    tic;
    if mt < 200000
        try
            [Ak, distA] = knnsearch(Xt, Xt, K=Kt + 1, NSMethod='exhaustive', Distance='fasteuclidean');
        catch
            [Ak, distA] = knnsearch(Xt, Xt, 'K', Kt + 1, 'NSMethod', 'exhaustive');
        end
    else
        [Ak, distA] = knnsearch(Xt, Xt, 'K', Kt + 1, 'NSMethod', 'kdtree');
    end
    if isequal(Ak(:, 1), (1:mt)') % Check that each point is its own NN
        Ak = Ak(:, 2:end);
    else % Swap around to make sure each point is its own NN
        % but before swapping, double the neighborhood size to increase
        % the likelyhood of the point being in its own neighborhood
        if mt < 200000
            try
                [Ak, distA] = knnsearch(Xt, Xt, K=min(ceil((3*mt)/2), 10*(Kt + 1)), NSMethod='exhaustive', Distance='fasteuclidean');
            catch
                [Ak, distA] = knnsearch(Xt, Xt, 'K', min(ceil((3*mt)/2), 10*(Kt + 1)), 'NSMethod', 'exhaustive');
            end
        else
            [Ak, distA] = knnsearch(Xt, Xt, 'K', min(ceil((3*mt)/2), 10*(Kt + 1)), 'NSMethod', 'kdtree');
        end
        diffA = find(Ak(:, 1) ~= (1:mt)');
        for i = 1:size(diffA)
            if ~ismember(diffA(i), Ak(diffA(i), :)) % if each point is not even in its own NNeigborhood then idk
                error('Neighbors:NotReflex', 'Distances were miscalculated.')
            else
                j = find(diffA(i) == Ak(diffA(i), :));
                Ak(diffA(i), 1:j) = flip(Ak(diffA(i), 1:j)); % all the distances from 1 to j should be 0 here so flipping only affects indices of points which are identical to the point we care about
            end
        end
        % return things to N + 1 size neighborhood
        Ak = Ak(:, 1:Kt+1);
        distA = distA(:, 1:Kt+1);
        % then delete each point from its neighborhood
        Ak = Ak(:, 2:end);
    end
    distA = distA(:, 2:end);
    tt(2) = toc;
    if ~bQ
        fprintf('Neighbors found in %.3f sec. \n', tt(2))
    end
    % sort is order big O n^2*log(n), so...
    
    
    % ------- Compute graph Laplacian
    sig = sg;
    if sg == 0
%         sg = sum(sum(distA)) / nnz(distA);
%         [st, me] = std(distA, 1, 'all');
%         sg = st + me;
        sg = prctile(distA,80, 'all', Method="approximate");
        sig = sg;
    end
    countA = nnz(distA);
    distA = exp(-(distA.^2)./(2*(sg.^2))); 
    if nnz(distA) < countA
        warning('Bad sigma choice reduced unsymmetric edges from %d to %d.', countA, nnz(distA))
    end
    ix = repmat((1:mt)', 1, Kt);
    W = sparse(ix(:), Ak(:), distA(:), mt, mt); % put dist into a sparse matrix        
    
    % W = spfun(@(x) exp(-x.^2 / 2*(sigma)^2), W);
    if getKDen
        kDen = sum(W, 2); % record and normalize the kernel density estimate
        kDen = kDen/(sum(kDen));
    else
        kDen = [];
    end

    W(W < eps) = 0; % dropping sufficiently weak connections
    W = max(W, W'); % now W is an adjacency matrix of an undirected graph
    % W = max(W, W');
    if getAff
        wt = W;
    else
        wt = [];
    end
    
    D = diag(sum(W, 2)); % degree matrix
    dg = diag(D);
    % deg = deg/(sum(deg));
    Lap = diag(diag(D).^(-1/2))*(D - W)*diag(diag(D).^(-1/2)); % normalized graph Laplacian
    Lap = max(Lap, Lap'); % make sure no symmetry is lost in normalization

    % output parsing --
    % output order: [Lap, sig, dg, kDen, wt]
    outs(1).A = Lap;
    outs(2).A = sig;
    outs(3).A = dg;
    % outs(4).A = tnc;
    outs(4).A = kDen;
    outs(5).A = wt;

    for iVar = 1:nout
        varargout{iVar} = outs(iVar).A;
    end
end