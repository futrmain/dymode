function [ eigenvalues, modes, energy, Sig ] = dmd_core( S )
%DMD_CORE Computes the DMD of snapshot matrix S
%   Detailed explanation goes here

%% SVD
[U, Sig, V] = svd(S(:,1:end-1), 0);
Sig = diag(Sig);


%% Invert Sig
Sigp = zeros(min(size(Sig, 1), size(Sig, 2)), 1);
threshold = eps('double') * size(Sigp, 1) * Sig(1, 1);
Sigp(Sig > threshold) = 1 ./ Sig(Sig > threshold);


%% Create B
B = U' * S(:, 2:end) * V * diag(Sigp);


%% Eigen problem
[X, eigenvalues] =  eig(B);
eigenvalues = diag(eigenvalues);

%% Weights
w = X\(U'*S(:,1));


%% Create the modes
modes = U * X;
modes = modes * diag(w);

% R = modes * fliplr(vander(eigenvalues)) - S(:,1:end-1);
% disp(max(max(abs(R))))


%% Compute energy
energy = sum(abs(modes).^2, 1);


% %% Get rid of half of the conjugate pairs
% modes = modes(:, imag(eigenvalues) >= 0);
% e_reduced = energy(imag(eigenvalues) >= 0);

%% Sort the modes by energy
[energy, i] = sort(energy, 'descend');
modes = modes(:, i);
eigenvalues = eigenvalues(i);


end

