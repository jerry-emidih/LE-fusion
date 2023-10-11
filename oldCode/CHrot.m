% Based on "Diffusion Maps for Changing Data" by Coifman and Hirn
% Code written by Jerry Emidih from a conversation with Alex Cloninger


function [rotA, R] = CHrot(A, B)
    % this rotates A into B using the knowledge that A and B have
    % orthonormal columns.
    R = A'*B;
    [u, ~, v] = svd(R);
    R = u*v';
    rotA = A*R;
end