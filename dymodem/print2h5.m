function [ last_index ] = print2h5( S, filename, dataset, snaps_per_file, start_index )
%PRINT2H5 Saves a Matlab Matrix into a set of HDF5 files.
%   This function saves a given Matrix into a set of hdf5 files compatible 
%   with dymode and dymodem. 
%
%   INPUT
%   S                   The matrix to be saved on disk
%   filename            A string containing the path and name of the files
%                       to be saved; e.g. '~/data/mydata'
%   dataset             The name of the dataset in which to save the data;
%                       e.g. '/snapshots_T'
%   snaps_per_file      The number of columns of the matrix (snapshots) to
%                       save in one file before creating a new file
%   start_index         The number of the first file that will be created.
%                       The default is 1. 
%
%   OUTPUT
%   last_index          The number of the last file that was saved on disk.
%%

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


