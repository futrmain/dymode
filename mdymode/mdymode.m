function [ eigenvalues, energy, modes, Sig ] = mdymode( filename, outdir, varargin )
%PRINT2H5 Summary of this function goes here
%   Detailed explanation goes here

%% Parse input
p = inputParser;

addRequired(p,'filename',@ischar);
addRequired(p,'outdir',@ischar);

addParameter(p,'geo',@ischar);
addParameter(p,'dataset','/snapshots_T',@ischar);
addParameter(p,'nfiles',1,@isnumeric);
addParameter(p,'stride',1,@isnumeric);
addParameter(p,'variables','u',@ischar);
addParameter(p,'modes',1,@isnumeric);
addParameter(p,'singulars',0,@isnumeric);
addParameter(p,'residuals','false',@ischar);

parse(p,filename,outdir, varargin{:});


%% Load data
[pathstr,name,~] = fileparts(p.Results.filename);
snaps = hdf2snaps( pathstr, name, p.Results.nfiles, p.Results.stride, p.Results.dataset, p.Results.variables );


%% Compute DMD
[ eigenvalues, modes, energy, Sig ] = dmd_core(snaps, ...
                                        'residuals', p.Results.residuals, ...
                                        'singulars', p.Results.singulars);

%% Save data
variables = strsplit(p.Results.variables, ',');
[success, mess, messid] = mkdir(p.Results.outdir);
if success == 1
    % Save light data
    % Not implemented yet
    
    % Save modes
    m = 1;
    msaved = 0;
    variable_file_list = '';
    
    while msaved < p.Results.modes && m <= size(modes, 2)
        while eigenvalues(m) < 0 && m < length(eigenvalues)
            m = m + 1;
        end
        if eigenvalues(m) >= 0
            for v = 1:length(variables)
                if strncmp(variables{v}, 'null', 4) == false
                    disp(['Writing mode' num2str(msaved)]);
                    
                    varfile = [p.Results.outdir '\mode' num2str(msaved, '%06i') '.' variables{v} '.abs']; 
                    fid = fopen(varfile, 'w+');
                    
                    writeEnsightHeader(fid, 'Module', msaved, variables{v});
                    % FIXME This should only print the part of modes
                    % corresponding to the right variable
                    fwrite(fid, abs(modes(:, m)), 'single');
                    
                    fclose(fid);
                    variable_file_list = [variable_file_list ...
                        'scalar per element: ' ...
                        variables{v} num2str(msaved) 'abs ' ...
                        'mode' num2str(msaved, '%06i') '.' variables{v} '.abs' sprintf('\n')];
                    
                    varfile = [p.Results.outdir '/mode' num2str(msaved, '%06i') '.' variables{v} '.ang']; 
                    fid = fopen(varfile, 'w+');
                    
                    writeEnsightHeader(fid, 'Angle', msaved, variables{v});
                    fwrite(fid, angle(modes(:, m)), 'single');
                    
                    fclose(fid);
                    variable_file_list = [variable_file_list ...
                        'scalar per element: ' ...
                        variables{v} num2str(msaved) 'ang ' ...
                        'mode' num2str(msaved, '%06i') '.' variables{v} '.ang' sprintf('\n')];
                    
                    msaved = msaved + 1;
                end
            end
        else
            m = m + 1;
        end
        
    end
    
    fid = fopen( [p.Results.outdir '/dmd.case'], 'w+' ); 
    
    fprintf(fid, 'FORMAT\n');
    fprintf(fid, 'type: ensight gold\n');
    fprintf(fid, 'GEOMETRY\n');
    fprintf(fid, 'model: dmd.geo\n');
    fprintf(fid, 'VARIABLE\n');
    fprintf(fid, '%s', variable_file_list);
    fprintf(fid, 'TIME\n');
    fprintf(fid, 'time set: 1 \n');
    fprintf(fid, 'number of steps: 1 \n');
    fprintf(fid, 'filename start number: 0 \n');
    fprintf(fid, 'filename increment: 1 \n');
    fprintf(fid, 'time values: \n');
    fprintf(fid, '0\n');
        
    fclose(fid);
else
    error(['Error, could not create ' p.Results.outdir]); 
end


end


function string = pad80(string)
if length(string) <= 79
    string = [string repmat(' ',1,80 - length(string)) ];
elseif length(string) > 79
    string = [string(1:79)  char(10)];
end
end

function writeEnsightHeader(fid, part, mode, var)

line = sprintf('%s of Mode %06d for %s', part, mode, var);
line = pad80(line);
% fprintf(fid, '%s\n', line);
fwrite(fid, line, 'char*1');

line = 'part';
line = pad80(line);
fwrite(fid, line, 'char*1');

fwrite(fid, 1, 'int16');
fwrite(fid, 0, 'int8');
fwrite(fid, 0, 'int8');

line = 'hexa8';
line = pad80(line);
fwrite(fid, line, 'char*1');

end