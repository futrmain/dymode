function [ eigenvalues, energy, modes, Sig ] = dymodem( filename, outdir, varargin )
%MDYMODE Matlab equivalent to dymode
%
%   Example:
%   mdymode('/scratch/my_file', '/scratch/results/', 'variables', 'u,null,null,p');
%
%   INPUT
%   filename        Path and rootname of the HDF5 files containing the snapshot matrix.
%   outdir          Path to the directory where to save results.
%
%   OPTION      DEFAULT         DESCRIPTION
%   geo                         Ensight geometry file.
%   datasert    'snapshots_T'   Dataset name inside the HDF5 files.
%   nfiles      1               Number of HDF5 files to read from.
%   stride      1               Stride between the snapshots to read from
%                               on disk.
%   variables   'u'             Name of the variables contained inside the
%                               data on disk, seprated by commas. Use 'null'
%                               to not read a variable from disk.
%   modes       1               Number of modes to save in Ensight format.
%   singulars   0               Number of singular values to display.
%   residuals   'false'         Set to 'true' to display residuals.
%   conj        'false'         Set to true to keep modes in the bottom half plane.
%   order       'descend'       Sort the modes in 'ascend'-ing or
%                               'descend'-ing order.
%   sorting     -11             Way to sort the modes:
%                               >= 0: Use the modes' energies at specified
%                                 timestep;
%                               -10:  Use the modes' mean energy over the sampled
%                                 period;
%                               -11:  Use the modes' median energy over the
%                                 sampled period.
%
%   OUTPUT
%   eigenvalues     Eigenvalues associated with the modes.
%   modes           Matrix of modes.
%   energy          Vector containing the energy of each mode.
%   Sig             Singular values.

%% Parse input
p = inputParser;

addRequired(p,'filename',@ischar);

addParameter(p,'outdir','',@ischar);

addParameter(p,'geo',@ischar);
addParameter(p,'dataset','/snapshots_T',@ischar);
addParameter(p,'nfiles',1,@isnumeric);
addParameter(p,'stride',1,@isnumeric);
addParameter(p,'variables','u',@ischar);
addParameter(p,'modes',1,@isnumeric);
addParameter(p,'singulars',0,@isnumeric);
addParameter(p,'residuals','false',@ischar);
addParameter(p,'sorting',-11,@isnumeric);
addParameter(p,'conj','false',@ischar);
addParameter(p,'order','descend',@ischar);

parse(p,filename,outdir, varargin{:});


%% Load data
[pathstr,name,~] = fileparts(p.Results.filename);
snaps = hdf2snaps( pathstr, name, p.Results.nfiles, p.Results.stride, p.Results.dataset, p.Results.variables );


%% Compute DMD
[ eigenvalues, modes, energy, Sig ] = dmd_core(snaps, ...
    'residuals', p.Results.residuals, ...
    'singulars', p.Results.singulars ...
    );


%% Sort the modes by energy
if p.Results.sorting >= 0
    % Use mode energy at specified time step
    energ_compensate = eigenvalues .^ p.Results.sorting;
    
elseif p.Results.sorting == -10
    % Use mode energy averaged over sampling period
    energ_compensate = abs(eigenvalues);
    energ_compensate(abs(eigenvalues) < 1 | abs(eigenvalues) > 1) ...
        =  (1 - abs(eigenvalues(abs(eigenvalues) < 1 | abs(eigenvalues) > 1)) .^ (2 * length(eigenvalues))) ...
        ./ ( (1 - abs(eigenvalues(abs(eigenvalues) < 1 | abs(eigenvalues) > 1)) .^ 2) * length(eigenvalues) );
    
elseif p.Results.sorting == -11
    % Use median mode energy
    energ_compensate = abs(eigenvalues) .^ (0.5 * (length(eigenvalues) - 1));
    
end

[~, i] = sort(energy .* energ_compensate', p.Results.order);
modes = modes(:, i);

if strncmp(p.Results.conj, 'false', 5)
    ev_sorted = eigenvalues(i);
    modes = modes(:, imag(ev_sorted) >= 0);
    clear ev_sorted
end


%% Save data
if ~isempty(p.Results.outdir)
    variables = strsplit(p.Results.variables, ',');
    [success, ~, ~] = mkdir(p.Results.outdir);
    if success == 1
        % Save light data
        % Not implemented yet
        
        % Save modes
        variable_file_list = '';
        
       NMODES = min(size(modes, 2), p.Results.modes);
        for k = 1:NMODES
            disp(['Writing mode' num2str(k-1)]);
            for v = 1:length(variables)
                if strncmp(variables{v}, 'null', 4) == false
                    varfile = [p.Results.outdir '\mode' num2str(k-1, '%06i') '.' variables{v} '.abs'];
                    fid = fopen(varfile, 'w+');
                    
                    writeEnsightHeader(fid, 'Module', k-1, variables{v});
                    % FIXME This should only print the part of modes
                    % corresponding to the right variable
                    fwrite(fid, abs(modes(:, k)), 'single');
                    
                    fclose(fid);
                    variable_file_list = [variable_file_list ...
                        'scalar per element: ' ...
                        variables{v} num2str(k-1) 'abs ' ...
                        'mode' num2str(k-1, '%06i') '.' variables{v} '.abs' sprintf('\n')];
                    
                    varfile = [p.Results.outdir '/mode' num2str(k-1, '%06i') '.' variables{v} '.ang'];
                    fid = fopen(varfile, 'w+');
                    
                    writeEnsightHeader(fid, 'Angle', k-1, variables{v});
                    fwrite(fid, angle(modes(:, k)), 'single');
                    
                    fclose(fid);
                    variable_file_list = [variable_file_list ...
                        'scalar per element: ' ...
                        variables{v} num2str(k-1) 'ang ' ...
                        'mode' num2str(k-1, '%06i') '.' variables{v} '.ang' sprintf('\n')];
                end
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