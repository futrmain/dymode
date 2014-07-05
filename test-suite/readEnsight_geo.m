%
% readEnsight_geoFile.m
% Original code by R. Furtrzynski
% Purpose: Read Ensight data from STARCCM+
%
% Usage: [Nnodes,coord,conn,np_nodesPerElement] = readEnsight_geoFile(geofile,partNumToRead);
%
% INPUTS
% geofile: string containing path to Ensight geometry file
% OUTPUTS (in order)
% Nnodes : total number of nodes in geometry part
% coord  : coordinates of nodes
% conn   : connectivity matrix
% np_nodesPerElement : array containing number of nodes for each element (only relevant for n-sided)
%
% See p.9 of Ensight 7 User Manual for format of Ensight Gold Geometry file (.geo) 
%
% Last modif.: DESC.                        		|DATE      |	WRITTEN BY
%	       Extended comments			14/11/2013     J. DAHAN		
%	       Deal with several element types		15/11/2013     J. DAHAN	
%              Works for n-sided elements    		09/12/2013     J. DAHAN
%
%

function [Nnodes,coord,conn,np_nodesPerElement] = readEnsight_geo(geofile)

element_types = {'point', 'g_point', ...
    'bar2', 'g_bar2', ...
    'bar3', 'g_bar3', ...
    'tria3', 'g_tria3', ...
    'tria6', 'g_tria6', ...
    'quad4', 'g_quad4', ...
    'quad8', 'g_quad8', ...
    'tetra4', 'g_tetra4', ...
    'tetra10', 'g_tetra10', ...
    'pyramid5', 'g_pyramid5', ...
    'pyramid13', 'g_pyramid13', ...
    'penta6', 'g_penta6', ...
    'penta15', 'g_penta15', ...
    'hexa8', 'g_hexa8', ...
    'hexa20', 'g_hexa20', ...
    'nsided', 'g_nsided', ...
    'nfaced', 'g_nfaced'};
element_points = [1 1 2 2 3 3 3 3 ...
    6 6 4 4 8 8 4 4 10 10 5 5 13 13 ...
    6 6 15 15 8 8 20 20 inf inf inf inf];

% Open .geo file
fid = fopen(geofile);

fprintf(1,'******* HEADER ******* \n');

% Read header
% strtrim removes insignificant whitespace
s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s); % C Binary
s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s); % description line 1
s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s); % description line 2
s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s); % node id <off/given/assign/ignore>
s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s); % element id <off/given/assign/ignore> 

% Read [extents] if they are present
s = char(fread(fid, 80, 'char*1')');    % [extents                                         
if strcmp(s(1:7), 'extents')            %  xmin xmax ymin ymax zmin zmax] (optional)
    extents = fread(fid, 6, 'float32');	% Contains xmin xmax ymin ymax zmin zmax for the part    
else
    fseek(fid, -80, 'cof');
end


%
% Need to include here a mechanism for multiple parts !!! Not yet
% implemented
%
maxElemTypes = 3; % max. number of element types to deal with
maxParts     = 2; % max. number of parts to deal with
conn               = cell([maxElemTypes maxParts]); % define connectivity vector array
np_nodesPerElement = cell([maxElemTypes maxParts]); % for n-sided data
conn_nFaced        = []; % define connectivity vector array for n-faced elements
np_nodesPerFace    = []; % for n-faced data

moreParts  = 1;
countParts = 0;
while moreParts

    countParts=countParts+1;
    
    s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s);                    % part
    if strcmp(s(1:4), 'part') ~= true
       disp('Error. It was expected to find ''part'' here.');
       fseek(fid, -80, 'cof');   % Go back the 80 char that were just read
       break;
    end

    part_nb = fread(fid, 1, 'int'); disp(['Part number: ' num2str(part_nb)]); % part number
    s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s);                    % description line
    s = strtrim(char(fread(fid, 80, 'char*1')')); disp(s);                    % coordinates

    Nnodes(countParts) = fread(fid, 1, 'int'); disp(['Number of nodes: ' num2str(Nnodes(countParts))]); % total number of nodes in the part

    % Now follows [id*Nnodes] x_coord*Nnodes y_coord*Nnodes z_coord*Nnodes
    % We read the first 3*Nnodes values. For now, assume there are no node IDs.
    coord{1,countParts} = fread(fid, 3*Nnodes(countParts), 'float'); % read X,Y,Z coords
    coord{1,countParts} = reshape(coord{1,countParts}, [Nnodes(countParts) 3]);    % Arrange X,Y,Z coordinates, by column

    % LOOP OVER ELEMENT TYPES FOR GIVEN PART
    moreElementTypes   = 1;  % boolean, true if file contains more element types that haven't been read yet
    countElTypes       = 0;  % counts number of different element types
   
    
    while moreElementTypes

        countElTypes = countElTypes + 1;

        elem_type = strtrim(char(fread(fid, 80, 'char*1')')); disp(['Type of elements: ' elem_type]);               % name of element type
        np = element_points(strcmp(elem_type, element_types)); disp(['Number of nodes per element: ' num2str(np)]); % num. of nodes per element (look-up)
        ne = fread(fid, 1, 'int'); disp(['Number of elements: ' num2str(ne)]);                                      % num. of elements of type elem_type
              
        % Follows [element_IDs*ne] connectivity_matrix*(ne x np)
        % (i.e. Number of elements x number of nodes per element)
        % Let's assume there are no IDs in the file
        if (np ~= inf)
            conn{countElTypes,countParts} = fread(fid, ne*np, 'int=>int');   % Read int, and don't default-convert to Matlab's double
            conn{countElTypes,countParts} = reshape(conn{countElTypes,countParts}, [np ne])';
        elseif strcmp(elem_type,'nsided');
            fprintf(1,'N-sided elements... \n');
            np_nodesPerElement{countElTypes,countParts} = fread(fid,ne,'int');
            for k=1:ne
               conn{countElTypes,countParts} = [conn{countElTypes,countParts}; fread(fid, np_nodesPerElement{countElTypes,countParts}(k), 'int=>int')];  
            end

            np_nPE_max = max(np_nodesPerElement{countElTypes,countParts});

            cursor=0; % Pad and reshape connectivity array for n-sided case
            for k=1:ne
                cursor=cursor+np_nodesPerElement{countElTypes,countParts}(k);
                if np_nodesPerElement{countElTypes,countParts}(k) < np_nPE_max
                    np_diff = np_nPE_max - np_nodesPerElement{countElTypes,countParts}(k); 
                    for q=1:np_diff
                        conn{countElTypes,countParts} = [conn{countElTypes,countParts}(1:cursor); conn{countElTypes,countParts}(cursor); conn{countElTypes,countParts}(cursor+1:end)];
                        cursor = cursor + 1;
                    end
                end
            end
            conn{countElTypes,countParts} = reshape(conn{countElTypes,countParts}, [np_nPE_max ne])';           
        elseif strcmp(elem_type, 'nfaced') % NOT FULLY FUNCTIONAL FOR NOW
            fprintf('N-faced elements... \n');
            nf_facesPerElement = fread(fid,ne,'int');
            NfacesTot          = sum(nf_facesPerElement);
            for k=1:ne
                np_nodesPerFace = [np_nodesPerFace; fread(fid, nf_facesPerElement(k), 'int=>int')];
                %conn = [conn; fread(fid, np_nsided_list(k), 'int=>int')];
            end
            for k=1:NfacesTot
                conn_nFaced = [conn_nFaced; fread(fid, np_nodesPerFace(k), 'int=>int')];
            end
     
        else
            fprintf(1,'Error, case not dealt with...\n');
            break;
        end

        % Check if there is a further element type in the geo file
        s = strtrim(char(fread(fid, 80, 'char*1')')); % read 80 characters
        if ~max(strcmp(s, element_types));
           disp('No more element types to read.');
            moreElementTypes = 0;
        end
        fseek(fid, -80, 'cof');   % Go back the 80 char that were just read
        
        % Avoid infinite loops
        if countElTypes > maxElemTypes
            disp('Error, too many element types, exiting while loop...');
            break;
        end

    end % end of looping over element types

    %conn = conn(~cellfun('isempty',conn)); % Remove empty elements of the connectivity cell array for vector cell
    %conn = mycell(~all(cellfun(@numel,conn)==0,2),:); % Remove empty elements of the connectivity cell array for matrix cell

    % Populate output argument cell array
    %indexOut = 1+(countParts-1)*4;
    %varargout{indexOut}   = Nnodes;
    %varargout{indexOut+1} = coord;
    %varargout{indexOut+2} = conn;
    %varargout{indexOut+3} = np_nodesPerElement;
    
    s = strtrim(char(fread(fid, 80, 'char*1')')); % part
    if strcmp(s(1:4), 'part') ~= true
      moreParts = 0; 
    end
    fseek(fid, -80, 'cof');   % Go back the 80 char that were just read    
    
end % end of looping over parts

% Let's assume we have everything. Count how many bytes are left in the
% file, just in case
[s, cnt] = fread(fid, 'char*1');
disp(['Warning: ' num2str(cnt) ' bytes (equivalent to ' num2str(cnt/4) ' floats) were ignored.']);

fclose(fid);
