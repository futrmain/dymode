function [ snaps ] = hdf2snaps( path, fileroot, n, s, dataset, var )
%READH5 Summary of this function goes here
%   Detailed explanation goes here

% clear all
% path = 'D:\DMD\DMD\x64\NNDEB';
% fileroot = 'Re350_oscillating';
% n = 2;
% s = 7;
% dataset = '/snapshots_T';
% var = {'u', 'null', 'null', 'p'};

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
    filedata = h5read(filename, dataset, [1 1], [327680 inf]);
    snaps = [snaps, filedata];
    
%     info = h5info(filename, dataset);
%     size = info.Dataspace.Size;
%     
%     varsize = size(1) / length(var);
%     filedata = zeros(nvar * varsize, ceil((size(2) - coffset+1) / s));
%     i = 1;
%     for v = 1:length(var)
%         if ~strncmp(var{v}, 'null', 4)
%             start = [(v - 1) * varsize + 1, coffset];
%             count = [varsize, Inf];
%             stride = [1, s];
%             filedata((i - 1) * varsize + 1:i * varsize, :) = ...
%                 h5read(filename, dataset, start, count, stride);
%             i = i + 1;
%         end
%     end
%     snaps = [snaps, filedata];
%     
%     off = mod(size(2) - (coffset - 1), s);
%     if off == 0
%         coffset = 1;
%     else
%         coffset = 1 + s - off;
%     end
end


end

