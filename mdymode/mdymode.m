function [ eigenvalues, energy, modes ] = mdymode( S )
%PRINT2H5 Summary of this function goes here
%   Detailed explanation goes here

format longg 

%% SVD
[U, Sig, V] = svd(S(:,1:end-1), 0);
Sig = diag(Sig);

disp('The first 5 singular values are:')
disp(Sig');

%% Invert Sig
Sigp = zeros(min(size(Sig, 1), size(Sig, 2)), 1);
threshold = eps('double') * size(Sigp, 1) * Sig(1, 1);
Sigp(Sig > threshold) = 1 ./ Sig(Sig > threshold);


%% Create B
B = U' * S(:, 2:end) * V * diag(Sigp);


%% Eigen problem
[X, eigenvalues] =  eig(B);


%% Weights
w = X\(U'*S(:,1));


%% Create the modes
modes = U * X;
modes = modes * diag(w);


%% 
energy = sqrt(sum(modes' * modes, 1));


end