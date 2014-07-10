timings = [];

benchdir = 'D:\DMD\DMD\matlabbench\u,null,null,null';

m = 5;
vars = {'u', 'null', 'null', 'null'};

dymodeexec = 'D:\Dymode_intel\x64\Release\Dymode_intel.exe';
inputpath = 'D:\DMD\channeldata';
fileroot = 'Re350_oscillating';
dataset = '/snapshots_T';
outputdir = 'D:\DMD\DMD\matlabbench';
geofile = 'D:\DMD\dmd.geo';



filelist = dir( benchdir );

nchecked = 0;

% Parse stride level
stride = [];
for i = 1:size(filelist,1)
    if isdir([benchdir '/' filelist(i).name])
        num = sscanf(filelist(i).name, '%i');
        if ~isempty(num)
            stride = [stride num]; %#ok<*AGROW>
        end
    end
end

for is = 1:length(stride)
    currdir = [benchdir '/' num2str(stride(is))];
    
    % Parse nfiles level
    filelist = dir( currdir );
    
    nfiles = [];
    for i = 1:size(filelist,1)
        if isdir([currdir '/' filelist(i).name])
            num = sscanf(filelist(i).name, '%i');
            if ~isempty(num)
                nfiles = [nfiles num];
            end
        end
    end
    
    for in = 1:length(nfiles)
        currdir = [benchdir '/' num2str(stride(is)) '/' num2str(nfiles(in))];
        
        
        tic;
        s_mat = hdf2snaps(inputpath, fileroot, nfiles(in), stride(is), dataset, vars);
        [eiv_mat, energy_mat, modes_mat, sing_mat] = mdymode(s_mat);
        mmin = min(m, size(modes_mat,2));
        modes_mat = modes_mat(:, 1:mmin);
        time_mat = toc;
        
        
        % Parse np level
        filelist = dir( currdir );
        
        np = [];
        for i = 1:size(filelist,1)
            if isdir([currdir '/' filelist(i).name])
                num = sscanf(filelist(i).name, '%i');
                if ~isempty(num)
                    np = [np num];
                end
            end
        end
        
        for p = 1:length(np)
            currdir = [benchdir '/' num2str(stride(is)) '/' num2str(nfiles(in)) '/' num2str(np(p))];
            % Parse block level
            filelist = dir( currdir );
            
            blocks = [];
            for i = 1:size(filelist,1)
                currdir = [benchdir '/' num2str(stride(is)) '/' num2str(nfiles(in)) '/' num2str(np(p))];
                if isdir([currdir '/' filelist(i).name])
                    num = sscanf(filelist(i).name, '%i');
                    if ~isempty(num)
%                         blocks = [blocks num];
                        currdir = [benchdir '/' num2str(stride(is)) '/' num2str(nfiles(in)) '/' num2str(np(p)) '/' num2str(num)];
                        
                        
                        try
                            modes_c = modesFromEnsight( currdir, vars );
                            modes_c = modes_c(:, 1:mmin);
                            energy_c = load([currdir '/spectrum.txt']);
                            singulars_c = load([currdir '/singulars.txt']);
                            eigenvalues_c = readEigencomplex([currdir '/eigenvalues.txt']);
                            [eigenvalues_c, ie] = sort(eigenvalues_c);
                            energy_c = energy_c(ie, 2);
                            
                            res_modes = max(abs(modes_c - modes_mat));
                            res_energy = max(abs(energy_c' - sqrt(energy_mat)));
                            res_eig = max(abs(eigenvalues_c - eiv_mat));
                            res_sings = max(singulars_c - sing_mat);
                            
                            disp(['NP: ' num2str(np(p)) ', stride: ' num2str(stride(is)) ', nfiles: ' num2str(nfiles(in)) ', block: ' num2str(num)]);
                            disp(['modes max diff:       ' num2str(res_modes)]);
                            disp(['energy max diff:      ' num2str(res_energy)]);
                            disp(['eigenvalues max diff: ' num2str(res_eig)]);
                            disp(['singulars max diff:   ' num2str(res_sings)]);
                            
                            nchecked = nchecked + 1;
                            if max(res_modes) > 1e-6 ...
                                    || res_energy > 1e-5 ...
                                    || res_eig > 1e-10 ...
                                    || res_sings > 1e-10
                                error('Criteria not matched')
                            end
                            
                            yamlfile = [currdir '/profiler-0.yml'];
                            if exist(yamlfile, 'file')
                                yamldata = ReadYaml(yamlfile);
                                time_c = mean(yamldata.Dymode.samples{:});
                                timings = [timings newdata];
                                disp(['Time Matlab: ' num2str(time_mat) ', time C++: ' num2str(time_c)]);
                            end
                        catch
                        end
                    end
                end
            end
        end
        
    end
end

nchecked

