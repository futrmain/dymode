clear all

check_only = 0;

dymodeexec = 'D:\Dymode_intel\x64\YDEBUG\Dymode_intel.exe';
inputpath = 'D:\DMD\channeldata';
fileroot = 'Re350_oscillating';
dataset = '/snapshots_T';
outputdir = 'D:\DMD\DMD\matlabbench';
geofile = 'D:\DMD\dmd.geo';

m = 5;
vars = {'u', 'null', 'null', 'null'};

stride = fliplr(floor(linspace(1, 20, 20)))
nfiles = floor(linspace(1, 10, 10))
block = floor(linspace(6, 256, 15));
block = [block 8 16 32 64 128 256];
block = unique(block)

% stride = (floor(linspace(3,3, 1)))
% nfiles = floor(linspace(3,3,1))
% block = floor(linspace(6, 60, 4));
% block = [block 8 16 32 64];
% block = unique(block)

for s = stride
    for n = nfiles
        s_mat = hdf2snaps(inputpath, fileroot, n, s, dataset, vars);
        [eiv_mat, energy_mat, modes_mat, sing_mat] = mdymode(s_mat);
        mmin = min(m, size(modes_mat,2));
        modes_mat = modes_mat(:, 1:mmin);
        
        for b = block
            % Full output path
            fulloutpath = sprintf('%s\\%s\\%i\\%i\\%i', ...
                outputdir, [sprintf('%s,', vars{1:end-1}) vars{end}], s, n, b);
            
            if ~check_only
                % Create the folder if it doesn't exist already.
                if ~exist(fulloutpath, 'dir')
                    mkdir(fulloutpath);
                    
                    copyfile(geofile, sprintf('%s\\%s\\%i\\%i\\%i\\dmd.geo', ...
                        outputdir, [sprintf('%s,', vars{1:end-1}) vars{end}], s, n, b));
                    
                    % create system command
                    command = sprintf('%s -n %i -s %i --block %i -e EigHess -m %i -f %s\\%s -o %s -i %s ', ...
                        ...np, ...
                        dymodeexec, ... dymode.exe
                        n, s, b, mmin, ... n s b m
                        inputpath, fileroot, ...
                        fulloutpath, ...
                        [sprintf('%s,', vars{1:end-1}) vars{end}] ...
                        );
                    disp(command);
                    [status, result] = system(command);
                    resultfile = sprintf('%s\\%s\\%i\\%i\\%i\\output.txt', ...
                        outputdir, [sprintf('%s,', vars{1:end-1}) vars{end}], s, n, b);
                    fid = fopen(resultfile, 'w+');
                    fprintf(fid, '%s', result);
                    fclose(fid);
                end
            end
            
            modes_c = modesFromEnsight( fulloutpath, vars );
            modes_c = modes_c(:, 1:mmin);
            energy_c = load([fulloutpath '\\spectrum.txt']);
            singulars_c = load([fulloutpath '\\singulars.txt']);
            eigenvalues_c = readEigencomplex([fulloutpath '\\eigenvalues.txt']);
            [eigenvalues_c, ie] = sort(eigenvalues_c);
            energy_c = energy_c(ie, 2);
            
            res_modes = max(abs(modes_c - modes_mat));
            res_energy = max(abs(energy_c' - sqrt(energy_mat)));
            res_eig = max(abs(eigenvalues_c - eiv_mat));
            res_sings = max(singulars_c - sing_mat);
            
            disp(['stride:               ' num2str(s) ', nfiles: ' num2str(n) ', block: ' num2str(b)]);
            disp(['modes max diff:       ' num2str(res_modes)]);
            disp(['energy max diff:      ' num2str(res_energy)]);
            disp(['eigenvalues max diff: ' num2str(res_eig)]);
            disp(['singulars max diff:   ' num2str(res_sings)]);
            
            if max(res_modes) > 1e-6 ...
                    || res_energy > 1e-5 ...
                    || res_eig > 1e-10 ...
                    || res_sings > 1e-10
            error('Criteria not matched')
            end
            
        end
    end
end

