


function [mapSource, lapCustomS, idxS, idxT, cSInd, cTInd, proTG] = clustAndMatch(mapSource, mapTarget, Y, clustNum, clustType, rotateMethod)

    
    % hopefully it isn't too bad?
    % inputs are mapSource, mapTarget, Y, clustNum, clustType

    % mapSource and mapTarget should be structs like the output of lapEig
    % Y is the matrix (points as rows) of out-of-sample points to be embedded in the source LE embedding
    % clustNum is number of clusters
    % clustType is one of these options: ["k-means", "dbscan", "linkage", "curve"]
    % "curve" refers to an experimental curvature-based clustering method
    % rotateMethod is "procrustes" or "scalerot", but scalerot doesn't work

    
    % how does matching work in the case of points which aren't sampled?
    % answer: the same way. we do the matching of the rep points with the
    % target set, then do OOS of the source points based on the new rep
    % locations
    

    % outputs are lapCustomS, idxS, idxT, cSInd, cTInd, mapSource.vecX, proTG

    % lapCustomS is the struct containing the embedded + matched OOS source points, now in the target LE space
    % idxS is cluster value for the source points, idxT is cluster value for the target points
    % cSInd and cTInd are binary versions of those (e.g. Row17Col4 being 1 means that the 17th point is in cluster 4)
    % mapSource.vecX has the coordinates of the embedded OOS points in the source space, before getting matched/rotated
    % proTG (can't remember why this is the name tbh, G for grayscale?) is a struct with the rotation info for each cluster (shift, scale, rotate matrix)

    %  -- basically whatever is needed for plotting and/or preimaging
    
    
    % mapSource = [];
    % mapTarget = [];
    % clustNum = 1;
    % clustType = 'k-means';
%     Y = [];
    % X = mapTarget.X;
    % Z = mapSource.X;
    
    m = size(mapTarget.vec, 1);
    n = size(Y, 1);
    nodims = size(mapTarget.vec, 2);
    
    % first, do OOS embedding into source space:
    mapSource.vecX = zeros(n, nodims);
    mapSource.vecX = out_of_sample(Y, mapSource); % everything else here
    
    % part 1: matching

    % repLocations = zeros(m, nodims);
    % for nn = 1:m
    %     src = mapSource.vec(nn, :);
    %     tar = mapTarget.vec(nn, :);
    % %     tar = mapC.vec(sol.P(nn), :);
    %     sc = norm(tar)/norm(src);
    %     [u, ~, v] = svd(tar'*src);
    %     r = v*u';
    %     repLocations(nn, :) = sc*src*r;
    % end
    repLocations = mapTarget.vec;
    
    mapOos = mapSource;
    mapOos.X(~mapSource.SampleSet, :) = Inf(nnz(~mapSource.SampleSet), size(mapSource.X, 2));
    % the above allows for control over which points form the basis of the
    % OOS extension
    mapOos.vec = repLocations;
    lapCustomS = out_of_sample(Y, mapOos); % everything else here
    
    % part 2: clustering 
    cN = clustNum;
    txtClust = clustType;
    idxS = zeros(n, 1);
    idxT = zeros(m, 1);
    rng(11+15+2022)
    
    
    % switch for clustering the rep points
    switch txtClust
        case "k-means"
            idxT = kmeans(mapSource.vec, cN, 'Replicates', 15); % clustering the rep set
            [idxS, ~] = clustVote(mapSource.vecX, mapSource.vec, idxT, 5);
        case "dbscan"
            idxT = dbscan(mapSource.vecX, mapSource.sigma/3, cN);
            oldCN = cN;
            cN = max(idxT); % the actual number of realized clusters
            
            % importantly, this also will refuse to cluster noisy points
            unclusteredCt = sum(idxG < 0);
            
            % as a start, create a cluster containing all the noisy points
            if unclusteredCt > 0
                idxT(idxT < 0) = (cN + 1)*ones(unclusteredCt, 1);
                cN = max(idxT);
            end
            [idxS, ~] = clustVote(mapSource.vecX, mapSource.vec, idxT, 5);
        case "linkage"
            txtLinkClust = "complete";
            lnk = linkage(pdist(mapSource.vec), txtLinkClust);
            idxT = cluster(lnk, 'maxclust', cN);
            oldCN = cN;
            cN = max(idxT); % the actual number of realized clusters
            [idxS, ~] = clustVote(mapSource.vecX, mapSource.vec, idxT, 5);
        case "curve"
            [idxT, mapSource.curvature] = createCurveClusters(mapSource.vec, cN);
            [idxS, ~] = clustVote(mapSource.vecX, mapSource.vec, idxT, 5);
    end
    % now assign clusters to the other points, generate binary indices
    
    cSInd = false(numel(idxS), cN);
    for nn = 1:cN
        cSInd(:, nn) = (idxS == nn);
    end
    
    cTInd = false(numel(idxT), cN);
    for nn = 1:cN
        cTInd(:, nn) = (idxT == nn);
    end
    
    % finally, do rotations
    [mapSource, proTG] = doRotation(mapSource, mapTarget, cN, cSInd, cTInd);
    
    
    %%
    function [mapS, proTG] = doRotation(mappingS, mappingT, cN, cSInd, cTInd)
        mappingS.rotX = zeros(size(mappingS.vecX)); 
        cNames = cell(1, cN);
        for nn = 1:cN
            cNames{nn} = char(sprintf('c%d', nn));
        end
        [nk, nodimsk] = size(mappingS.vecX, [1 2]); 
        
        proTG = struct();
        sh = zeros(nk, nodimsk); % shift matrix, probably inefficient
        if rotateMethod == "scalerot"
        % for nn = 1:cN
        %     na = cNames{nn};
            
        %     proTG.(na) = scaleRot(mappingT.vec(cTInd(:, nn), :), mappingS.vec(cTInd(:, nn), :));
        %     %     % try without the medoids
        %     muS = mean(mappingS.vec(cTInd(:, nn), :), 1);
        %     muOOS = mean(mappingS.vecX(cSInd(:, nn), :), 1);
        %     sh(cSInd(:, nn), :) = repmat(muS - muOOS, [nnz(cSInd(:, nn)), 1]); % shifting the current clusters' oos points only
        %     mappingS.rotX(cSInd(:, nn), :) = ( mappingS.vecX(cSInd(:, nn), :) + sh(cSInd(:, nn), :) )*proTG.(na).T +...
        %     proTG.(na).c(1, :); % then apply to each cluster in full
        % end
            
            return % scaled rotation doesn't really work currently
        elseif rotateMethod == "procrustes"
            for nn = 1:cN
                na = cNames{nn};
                [~, ~, proTG.(na)] = procrustes(mappingT.vec(cTInd(:, nn), :), mappingS.vec(cTInd(:, nn), :), 'Reflection', 'best', 'Scaling', false);
                muS = mean(mappingS.vec(cTInd(:, nn), :), 1);
                muOOS = mean(mappingS.vecX(cSInd(:, nn), :), 1);
                sh(cSInd(:, nn), :) = repmat(muS - muOOS, [nnz(cSInd(:, nn)), 1]); % shifting the current clusters' oos points only
                mappingS.rotX(cSInd(:, nn), :) = ( mappingS.vecX(cSInd(:, nn), :) + sh(cSInd(:, nn), :) )*proTG.(na).T +...
                    proTG.(na).c(1, :); % then apply to each cluster in full
            end
        end
        
        mappingS.txtRot = txtClust + " clustered ";
        mappingS.cNames = cNames;
        mapS = mappingS;
        % this is also using the new scaleRot function
    end
    
    %%
    function [idxS, crv] = createCurveClusters(meshPts, cN)
        % need to enforce uniqueness
        [~, uniPts, ~] = unique(meshPts, 'rows', 'stable');
    %     meshps = bigIndex(meshPts);
    %     uniPts = meshps(uniPts);
        nk = size(meshPts, 1);
        meshIdx = ismember((1:nk)', uniPts);
        X = meshPts(meshIdx, :);
        vD = [1 2 3];
        bound = boundary(X(:, vD(1:2)), .5);
        bound = bound(1:end -1); % drop the loop point
        if size(bound, 2) == 1
            bound = [bound, circshift(bound, -1)];
        end
        crv = struct('L', 0, 'R', 0, 'K', -Inf);
        path = X([bound(:, 1); bound(1, 1)], vD(1:2));
        [crv.L, crv.R, crv.K] = curvature(path);
    
        [~, maxJump] = max(diff(crv.L));
        bound = circshift(bound, -maxJump);
        
        % now repeat the above
        crv = struct('L', 0, 'R', 0, 'K', -Inf);
        path = X([bound(:, 1); bound(1, 1)], vD(1:2));
        [crv.L, crv.R, crv.K] = curvature(path);
    
        bdCt = numel(unique(bound)) + 1; % add one to close the loop
        [~, radIdx] = sort(crv.R, 'ascend');
        halfSmall = crv.R < crv.R(radIdx(floor(bdCt/2)));
        posSecond = [false; false; diff(crv.R, 2) > 0];
        gradCurve = [0; diff(crv.R, 1)];
        zeroCross = gradCurve.*circshift(gradCurve, 1) < 0;
        clustSep = intersect(find(halfSmall), find(zeroCross & posSecond) - 1);
    %     cN = numel(clustSep);
        
        % for manual cluster size:
    %     cN = 7;
        [~, subSep] = sort(crv.R(clustSep(2:end)), 'ascend');
        clustSep = sort(clustSep(subSep(1:cN)));
        clustSep = [0; clustSep]; % padding a zero in front for indexing ease
        pathInd = bound(:, 1);
        idxS = zeros(nk, 1);
        % first, cluster the boundary
        % cluster n ends at clustSep(n)
        for nn = 1:cN - 1
            segInd = pathInd(clustSep(nn) + 1 : clustSep(nn+1));
            idxS(uniPts(segInd)) = nn;
        end
        idxS(uniPts(pathInd(clustSep(cN) + 1 : end))) = cN; % assign everything after the last separation to the last cluster
        % maybe this is bad, idk idc
        bPathInd = ismember((1:nk)', uniPts(pathInd));
        [idxS(~bPathInd), ~] = clustVote(meshPts(~bPathInd, vD(1:2)), meshPts(bPathInd, vD(1:2)), idxS(bPathInd), 5); % assign clusters to the other points
        % then, every other point is clustered based on nearest neighbors
    end

end