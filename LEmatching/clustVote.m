function [winner, nbhd] = clustVote(newPoints, training, clusters, voteSize)
%% Assigns cluster labels to points in X (newpoints) by using the points in Y (training)

    if size(training, 1) ~= size(clusters, 1)
        error('Not all of the training points are labeled, or there is a size mismatch. \n')
    end
    Ak = knnsearch(training, newPoints, 'K', voteSize, 'NSMethod', 'exhaustive'); % nearest neighbor indices
    C = reshape(clusters(Ak(:)), size(Ak)); % cluster labels of those neighbors, in order of closeness
    [M, F] = mode(Ak, 2); % get most frequent label for each point
    ties = F <= voteSize/2; % is the most frequent label a majority?
    winner = zeros(size(newPoints, 1), 1);
    winner(~ties) = M(~ties); % if a point has a majority label, then assign it
    winner(ties) = C(ties, 1); % if there's not a majority, then assign the label of the nearest neighbor
    nbhd = Ak;