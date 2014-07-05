function [ Modes ] = modesFromEnsight( path, var )
%MODESFROMENSIGHT Summary of this function goes here
%   Detailed explanation goes here

% clear all
% path = 'D:\DMD\DMD\x64\NNDEB\out';
% var = {'u', 'p'};

nvars = 0;
for v = 1:length(var)
    if ~strncmp(var{v}, 'null', 4)
        nvars = nvars + 1;
    end
end
filelist = dir( path );

%Parse filenames in directory
n = 0;
for f = 1:size(filelist,1)
    [i, ~, e] = sscanf(filelist(f).name, ['mode%06i.' var{1} '.%*s']);
    if isempty(e)
        n = n + 1;
        num(n) = i; 
%         disp(num); disp(var); 
    end   
end

% Get list of variables and mode numbers
num = unique(num);

A = zeros(327680*nvars, length(num));
P = zeros(327680*nvars, length(num));
for m = 0:length(num)-1
    i = 1;
    for v = 1:length(var)
        if ~strncmp('null', var{v}, 4)
            filename = sprintf('%s\\mode%06i.%s.abs', path, m, var{v});
            A((i-1)*327680+1:i*327680,m+1) = readEnsight_variable(filename, {327680});
            filename = sprintf('%s\\mode%06i.%s.ang', path, m, var{v});
            P((i-1)*327680+1:i*327680,m+1) = readEnsight_variable(filename, {327680});
            
            i = i + 1;
        end
    end
end

Modes = A .* exp(1i .* P);
end