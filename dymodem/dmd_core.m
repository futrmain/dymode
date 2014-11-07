function [ eigenvalues, modes, energy, Sig ] = dmd_core( S, varargin )
%DMD_CORE Computes the DMD of snapshot matrix S
%
%   dmd_core(S, 'singulars', 5, 'residual', 'true');
%
%   INPUT
%   S               Snapshot matrix. Must have at least 2 columns.
%   
%   OPTION      DEFAULT     DESCRIPTION
%   singulars   0           Number of singular values to display.
%   residual    'false'     Setting 'true' displays the residuals.
%
%   OUTPUT
%   eigenvalues     Eigenvalues associated with the modes.
%   modes           Matrix of modes, columns are sorted by mode energy.
%   energy          Vector containing the energy of each mode.
%   Sig             Singular values.


%% Parse input
p = inputParser;

addRequired(p,'snaps');
addParameter(p,'singulars',0,@isnumeric);
addParameter(p,'residuals','false',@ischar);

parse(p,S, varargin{:});


%% SVD
[U, Sig, V] = svd(S(:,1:end-1), 0);

if strcmp(p.Results.residuals, 'true') 
    r = max(max(U*Sig*V' - S(:,1:end-1)));
    disp (['Residual from SVD: ' num2str(r)]);
end

Sig = diag(Sig);

NSING = min(length(Sig), p.Results.singulars);
if NSING > 0
    disp (['First' num2str(NSING) ' singular values: ']);
    disp (Sig(1:NSING)');
end


%% Invert Sig
Sigp = zeros(min(size(Sig, 1), size(Sig, 2)), 1);
threshold = eps('double') * size(Sigp, 1) * Sig(1, 1);
Sigp(Sig > threshold) = 1 ./ Sig(Sig > threshold);


%% Create B
B = U' * S(:, 2:end) * V * diag(Sigp);

clear V Sigp 


%% Eigen problem
[X, eigenvalues] =  eig(B);
eigenvalues = diag(eigenvalues);

if strcmp(p.Results.residuals, 'true') 
    r = max(max(abs(B*X - X*diag(eigenvalues))));
    disp (['Residual from eigen problem: ' num2str(r)]);
end

clear B

%% Weights
w = X\(U'*S(:, 1));

clear S


%% Create the modes
modes = U * X;
clear U X

modes = modes * diag(w);
clear w

if strcmp(p.Results.residuals, 'true') 
    r = max(max(abs(modes*fliplr(vander(eigenvalues)) - S(:, 1:end-1))));
    disp (['Residual from modes: ' num2str(r)]);
end
% R = modes * fliplr(vander(eigenvalues)) - S(:,1:end-1);
% disp(max(max(abs(R))))




%% Compute energy
energy = sum(modes .* conj(modes), 1);

% %% Get rid of half of the conjugate pairs
% modes = modes(:, imag(eigenvalues) >= 0);
% e_reduced = energy(imag(eigenvalues) >= 0);


end

