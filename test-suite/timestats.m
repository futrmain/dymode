timings = [];

benchdir = 'D:\DMD\DMD\matlabbench\u,null,null,null';

filelist = dir( benchdir );

% Parse stride level 
stride = [];
for i = 1:size(filelist,1)
    if isdir([benchdir '\' filelist(i).name])
        num = sscanf(filelist(i).name, '%i');
        if ~isempty(num)
            stride = [stride num]; %#ok<*AGROW>
        end
    end   
end

for is = 1:length(stride)
    currdir = [benchdir '\' num2str(stride(is))];
    
    % Parse nfiles level
    filelist = dir( currdir );
    
    nfiles = [];
    for i = 1:size(filelist,1)
        if isdir([currdir '\' filelist(i).name])
            num = sscanf(filelist(i).name, '%i');
            if ~isempty(num)
                nfiles = [nfiles num];
            end
        end
    end
    
    for in = 1:length(nfiles)
        currdir = [benchdir '\' num2str(stride(is)) '\' num2str(nfiles(in))];
       
        % Parse block level
        filelist = dir( currdir );
        
        blocks = [];
        for i = 1:size(filelist,1)
            if isdir([currdir '\' filelist(i).name])
                num = sscanf(filelist(i).name, '%i');
                if ~isempty(num)
                    blocks = [blocks num];
                    yamlfile = [currdir '\' num2str(num) '\profiler-0.yml'];
                    if exist(yamlfile, 'file')
                        yamldata = ReadYaml(yamlfile);
                        newdata = [stride(is); nfiles(in); num; ...
                            mean(yamldata.Dymode.subFunctions.Computations.samples{:})];
                        timings = [timings newdata];
                    end
                end
            end
        end
        
        
    end
end

%%
close all
nsize = timings(2, :) ./ timings(1, :);

scatter(timings(3, :), timings(4, :)./nsize.^1.5 , 2+timings(2, :).^2, timings(1, :))
a = gca;
set(a, 'YScale', 'log')

%%
figure 
t = timings(4, :);
scatter(nsize, timings(4, :) , 2+timings(2, :).^2, timings(1, :))