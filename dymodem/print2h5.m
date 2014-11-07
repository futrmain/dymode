function [ last_index ] = print2h5( S, filename, dataset, snaps_per_file, start_index )
%PRINT2H5 Summary of this function goes here
%   Detailed explanation goes here

% Default values
if ~exist('dataset','var') || isempty(dataset)
  dataset = 'snapshots_T';
end

if ~exist('snaps_per_file','var') || isempty(snaps_per_file)
  snaps_per_file = size(S,2);
end

if ~exist('start_index','var') || isempty(start_index)
  start_index = 0;
end


nfiles = ceil(size(S,2) / snaps_per_file);
for fnum = 1:nfiles
    fname = sprintf('%s%04i.h5', filename, fnum + start_index);
    
    s = S(:, (f-1) * snaps_per_file + 1 : min(f * snaps_per_file, end));
    hdf5write(fname, dataset, s);
end

last_index = nfiles + start_index;

end


