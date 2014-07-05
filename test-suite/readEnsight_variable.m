%
% readEnsight_var.m
% Read ensight variable file
% Ensight variable files can be either per_node or per_element
% Here assume the file to be read is per_element
% For now, only deals with single part
%
% usage: sca_elem = readEnsight_var(varfile,ne,partNumToRead);
% varfile: string containing the path to the variable file
% ne:      number of nodes/elements per element type (array of size ntypes)
%

function sca_elem = readEnsight_variable(varfile,ne)

fid = fopen(varfile);

s = strtrim(char(fread(fid, 80, 'char*1')')); %disp(s);                          % Description line 1

% maxElemTypes = 1;                                                               % max. number of element types to deal with
% numParts     = size(ne,2);
% sca_elem     = cell([maxElemTypes numParts]);                                   % define connectivity vector array
% sca_elem = single([]);
             
    s = strtrim(char(fread(fid, 80, 'char*1')')); %disp(s);                      % part
    part_nb = fread(fid, 1, 'int'); %disp(['Part number: ' num2str(part_nb)]);   % part number
    
%     for k=1:length(ne{countParts})
        s = strtrim(char(fread(fid, 80, 'char*1')')); %disp(s);                  % Element type
        sca_elem = fread(fid, ne{1}, 'float');
%     end
    

%sca_elem = sca_elem(~cellfun('isempty',sca_elem)); % Remove empty elements of the connectivity cell array

fclose(fid);

end
