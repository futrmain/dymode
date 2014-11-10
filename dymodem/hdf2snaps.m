function [ snaps ] = hdf2snaps( path, fileroot, n, s, dataset, var )
%hdf2snaps Create a snapshot matrix from an hdf5 file sequence
%   hdf2snaps reads the data contained inside a series of hdf5 files and
%   returns it as a matrix output. This function reads the same type of
%   hdf5 files as those used by dymode.
%
%   INPUT
%   path            The path to the directory containing the hdf5 files
%   fileroot        The name of the hdf5, without the file number
%   n               The number of files to be read; i.e. the files
%                   <fileroot>0001.h5 to <fileroot><n>.h5 will be read
%   s               The stride between the snapshots to read
%   dataset         The name of the dataset containing the snapshot data,
%                   e.g. '/snapshots_T'
%   var             The list of variables making up the rows of the matrix.
%                   Use "null" to skip reading blocks of rows; e.g. use 
%                   'a,null' to read only the first half of the rows
%
%   OUTPUT
%   snaps           A MAtlab matrix containing the snapshos read from disk.
%%

nvar = 0;
var = strsplit(var, ',');
for v = 1:length(var)
    if ~strncmp(var{v}, 'null', 4)
        nvar = nvar +1;
    end
end

snaps = [];
coffset = 1;
for f = 0+1:n+0
    filename = sprintf('%s\\%s%04i.h5', path, fileroot, f);
    disp(filename)
    
    info = h5info(filename, dataset);
    size = info.Dataspace.Size;
    
    varsize = size(1) / length(var);
    filedata = zeros(nvar * varsize, ceil((size(2) - coffset+1) / s));
    i = 1;
    for v = 1:length(var)
        if ~strncmp(var{v}, 'null', 4)
            start = [(v - 1) * varsize + 1, coffset];
            count = [varsize, Inf];
            stride = [1, s];
            filedata((i - 1) * varsize + 1:i * varsize, :) = ...
                h5read(filename, dataset, start, count, stride);
            i = i + 1;
        end
    end
    snaps = [snaps, filedata];
    
    off = mod(size(2) - (coffset - 1), s);
    if off == 0
        coffset = 1;
    else
        coffset = 1 + s - off;
    end
end


end

