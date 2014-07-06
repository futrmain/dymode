function [ encas_file_struct ] = loadEncasFile( varargin )
%LOADENCASFILE Loads case/ encas file and puts all variables found in the
% encas file in a single tree-like structure.
%
% USAGE: (Loading files)
% [ encas_file_struct ] = loadEncasFile( filelocation )
% [ encas_file_struct ] = loadEncasFile; % asks for encas/ case file.
%
%
% EXAMPLES: (To acces the data)
%            geometry = encas_file_struct.geometry
%     var_wss_on_wall = encas_file_struct.wall_shear.surface3.timestep1
%  velocity_timesteps = encas_file_struct.x_velocity.timesteps
% timing_of_timesteps = encas_file_struct.encas_file.timings
%
%
% Please note that these function will probably not work on ALL case/encas
% files.  (Only tested with selected encas files exported by Ansys Fluent: 
%           transient simulation (>1 TIMESTEPS) & (single) static geometry)
%
% Intended as a guideline for loading encas files in matlab; feel free to
% edit/ improve / add features/ ask questions.
% 
%
% Author: Wouter Potters, Radiology dpt, Academic Medical Center, Amsterdam
% Email: w.v.potters@amc.nl
% Date January 29, 2013
%
% Many thanks to Paul Groot (AMC, Amsterdam) for his tips
% regarding the binary file formats.
%
%
%
% #####
% VISUALISATION EXAMPLE:
% #####
% show triangles from surface in red
% f1 = figure(1); set(f1,'color','w'); cla; % make new figure window; empty axes
% patch('vertices',encas_contents.geometry.wall.vertices,... % select wall vertices
%       'faces',encas_contents.geometry.wall.faces_tria3,... % select wall triangle faces
%       'facecolor','r',... % red color of triangle
%       'edgecolor','k');  % with black edge color
% view(3); axis equal vis3d off; rotate3d on % show it in a 3d fashion.
%  
% %% show triangles from surface with WSS colors
% f2 = figure(2); set(f2,'color','w'); cla; % make new figure window; empty axes
% patch('vertices',encas_contents.geometry.wall.vertices,... % select wall vertices
%       'faces',encas_contents.geometry.wall.faces_tria3,... % select wall triangle faces
%       'facecolor','interp',... % interp face color from vertices
%       'cdata',(encas_contents.wall_shear.wall.timestep1),... %select ColorData as double datatype from wal_shear timestep1
%       'edgecolor','k');  % with black edge color
% colorbar;
% view(3); axis equal vis3d off; rotate3d on % show it in a 3d fashion.
%  
%  
% %% show QUAD ELEMENT from INT_BLOOD
% f3 = figure(3); set(f3,'color','w'); cla; % make new figure window; empty axes
% patch('vertices',encas_contents.geometry.int_blood.vertices,... % select wall vertices
%       'faces',encas_contents.geometry.int_blood.faces_quad4,... % select wall QUAD4 faces
%       'facecolor','r',... % red color
%       'facealpha',0.5,... % with transparency 50%
%       'edgecolor','k');  % with black edge color
% view(3); axis equal vis3d off; rotate3d on % show it in a 3d fashion.
%  
% hold on; % add another patch on top of it
% patch('vertices',encas_contents.geometry.int_blood.vertices,... % select wall vertices
%       'faces',encas_contents.geometry.int_blood.faces_tria3,... % select wall TRIANGLE faces
%       'facecolor','b',... % blue color
%       'facealpha',0.3,... % with transparency 30%
%       'edgecolor','k');  % with black edge color



version = ver('matlab');
if str2double(version.Version) > 8
    narginchk(0,1); nargoutchk(1,1);
end

    switch nargin
        case 0, validfile = 0;
        case 1, filelocation = varargin{1}; if exist(filelocation,'file') == 2, validfile = 1; else validfile = 0; end
    end
    if ~validfile % if no (valid) filepath is provided; let the user select one
        [filename,pathtofile] = uigetfile({'*.encas;*.case','encas/case files'},'Select case/encas file...',filelocation);
        filelocation = fullfile(pathtofile,filename);
    end
    [pathtofile] = fileparts(filelocation);

    % initialize struct and add encas file
    encas_file_struct                     = struct();
    encas_file_struct.encas_file.path     = filelocation;
    istimeresolved_simulation = checkTimeResolvedSimulation(encas_file_struct.encas_file.path);
    [encas_file_struct.encas_file.filelist, encas_file_struct.encas_file.timings] = loadCaseFile(encas_file_struct.encas_file.path,istimeresolved_simulation);


    % Load all variable-files found in the encas/case file
    fnames = (fieldnames(encas_file_struct.encas_file.filelist));
    for fname_i = 1:length(fnames) % load all datafiles

        % Set path and type for current file in encas
        encas_file_struct.(fnames{fname_i}).path = fullfile( pathtofile, encas_file_struct.encas_file.filelist.(fnames{fname_i}).file );

        % Load contents of current file in encas
        switch encas_file_struct.encas_file.filelist.(fnames{fname_i}).type

            case 'model'  %geometry
                encas_file_struct.(fnames{fname_i})      = loadGeoFile( encas_file_struct.(fnames{fname_i}).path, istimeresolved_simulation );
                encas_file_struct.(fnames{fname_i}).type = encas_file_struct.encas_file.filelist.(fnames{fname_i}).type;

            case 'scalar' %pressure, x_wall_shear, x_velocity, wall_shear, velocity_magnitude
                encas_file_struct.(fnames{fname_i})      = loadScalarFile( encas_file_struct.(fnames{fname_i}).path, encas_file_struct.geometry, istimeresolved_simulation );
                encas_file_struct.(fnames{fname_i}).type = encas_file_struct.encas_file.filelist.(fnames{fname_i}).type;

            case 'vector' %velocity
                encas_file_struct.(fnames{fname_i})      = loadVectorFile( encas_file_struct.(fnames{fname_i}).path, encas_file_struct.geometry, istimeresolved_simulation );
                encas_file_struct.(fnames{fname_i}).type = encas_file_struct.encas_file.filelist.(fnames{fname_i}).type;

            otherwise % TENSOR, TENSOR9, COMPLEX SCALAR and COMPLEX VECTOR - UNSUPPORTED // feel free to add them yourself :)
                encas_file_struct.(fnames{fname_i}) = load_OTHER_file(pathtofile);
                disp([fnames{fname_i} ' not loaded'])

        end
    end
    
    checkTimings(encas_file_struct);


    function [casefilecontents, timings] = loadCaseFile(pathtofile,istimeresolved_simulation)
        % LOADCASEFILE This subfunction loads the contents of the casefile
        % USAGE:  [casefilecontents_struct] = loadCaseFile(pathtofile)
        casefilecontents = struct; variable = ''; %#ok<NASGU>
        fid = fopen(pathtofile);
        while ~feof(fid)
            txtline = fgetl(fid);
            switch txtline
                case 'GEOMETRY'
                    geometry = fgetl(fid);
                    istimeresolved_simulation = 0;
                    switch istimeresolved_simulation
                        case 1
                            res = textscan(geometry,'%s %f  %f %s');
                            casefilecontents.geometry.file = res{4}{1};
                            casefilecontents.geometry.type = res{1}{1}(1:end-1);
                        case 0
                            res = textscan(geometry,'%s %s');
                            casefilecontents.geometry.file = res{2}{1};
                            casefilecontents.geometry.type = res{1}{1}(1:end-1);
                    end
                    
                case 'VARIABLE'
                    while ~feof(fid)
                        curpos = ftell(fid); % save position
                        variable = fgetl(fid);
                        if strcmp(variable,'TIME'), fseek(fid,curpos,-1); break; end % goto position if TIME is found
                        switch istimeresolved_simulation
                            case 1
                                res = textscan(variable,'%s per node: %f  %f %s %s');
                                casefilecontents.(res{4}{1}).file = res{5}{1};
                                casefilecontents.(res{4}{1}).type = res{1}{1};
                            case 0
                                res = textscan(variable,'%s per node: %s %s');
                                casefilecontents.(res{2}{1}).file = res{3}{1};
                                casefilecontents.(res{2}{1}).type = res{1}{1};
                        end
                    end
                    
                case 'TIME'
                    nr_of_models     = sscanf(fgetl(fid),'time set: %f Model');
                    if nr_of_models > 1
                        error('Only single timestep (static) geometries are supported by this function.');
                    end
                    nr_of_timesteps  = sscanf(fgetl(fid),'number of steps: %f');
                    timestep_values = [];
                    fgets(fid,13); %ignore 13 characters: 'time values: '
                    nrit = 1;
                    while length(timestep_values) < nr_of_timesteps && ~feof(fid) && nrit <= nr_of_timesteps  % find all timesteps
                        add_these_timesteps = sscanf(fgetl(fid),'%f');
                        timestep_values = [timestep_values; add_these_timesteps]; %#ok<AGROW>
                        nrit = nrit +1; % safety check; avoid infinite looping for finding timesteps...
                    end
                    timings = timestep_values;
                    
                otherwise
                    % ignore line
                    % disp(['otherwise: ' txtline])
                    
            end
        end
        if istimeresolved_simulation == 0
            timings = [];
        end
        fclose(fid);
    end                    % loads .case/.encas file

    function geofilecontents = loadGeoFile(pathtofile,istimeresolved_simulation)
        % http://vis.lbl.gov/NERSC/Software/ensight/docs8/UserManual.pdf
        % load geometry file; see below for details on the file format
        fid = fopen(pathtofile,'r');
        geofilecontents.type   = fread(fid,[1 80],'*char'); % C Binary                                            80 chars
        if istimeresolved_simulation == 1
            [~]                         = fread(fid,[1 80],'*char') ; % BEGIN TIME STEP                               80 chars % IGNORE
        end
        geofilecontents.description = fread(fid,[80 2],'*char')'; % description lines                           2*80 chars
        
        [~]                         = fread(fid,[1 80],'*char'); % node id <off/given/assign/ignore>                   80 chars % IGNORE
        [~]                         = fread(fid,[1 80],'*char'); % element id <off/given/assign/ignore>                80 chars % IGNORE
        
        if ~strcmp(geofilecontents.type(1:8),'C Binary'), error('This code only works for C binary data'); end
        
        index = 1;
        while true && ~feof(fid)
            part   = fread(fid,[1 80],'*char');              % [part                                               80 chars
            if ~isempty(part) && strcmp(part(1:4),'part')  % loop over all parts of the geometry; surface mesh and elements
                
                % STEP 1: LOAD ALL DATA FOR THE CURRENT ELEMENT
                part_number    = fread(fid,                 1, '*uint32','l'); % #                                                    1 int
                part_descr     = fread(fid,            [1 80], '*char'      ); % description line                                    80 chars
                [~]            = fread(fid,            [1 80], '*char'      ); % coordinates                                         80 chars % IGNORE
                nr_of_vertices = fread(fid,                 1, '*uint32','l'); % nn                                                   1 int
                vertices_xyz   = fread(fid,[nr_of_vertices 3], '*float' ,'l'); % x y z coordinates // index is skipped
                faces_type     = strcat(fread(fid,            [1 80], '*char'      )); % get the type of element
                nr_of_faces    = fread(fid,                 1, '*uint32','l'); % ne                                                   1 int
                switch faces_type(1:5)  % compare first 5 characters in faces_type
                    case {'tetra'}, cornerpts = 4;    % tetraeder; 4 cornerpoints
                    case 'tria3', cornerpts = 3;    % triangle;  3 cornerpoints
                    case 'penta', cornerpts = 6;
                    case 'quad4', cornerpts = 4;
                    otherwise, error('Unknown elementtype, only tetra4 and tria3 are implemented.');
                end
                element_faces = fread(fid,[cornerpts nr_of_faces],'*uint32','l')';
                
                % STEP 2: SAVE DATA FOR CURRENT ELEMENT IN STRUCT
                geofilecontents.(part_descr).vertices     =  vertices_xyz;
                geofilecontents.(part_descr).(['faces_' faces_type]) = element_faces;
                geofilecontents.(part_descr).part_nr      =   part_number;
                
                index = index + 1; % next element
            elseif ~isempty(part) && strcmp(part(1:8),'END TIME') % disp('END TIME STEP reached; check if this the end of the file?')
                curpos  = ftell(fid); %save current position
                fread(fid,1); %try reading 1 more line
                if feof(fid) % check if eof is reached
                    geofilecontents.nr_parts = index; % save nr of parts
                    break; %eof reached; break while loop;
                else % No; return to previous curpos
                    if (fseek(fid,curpos,'bof') ~= 0)
                        error('could not go to defined position...')
                    end
                end
            elseif length(part) >= 6 && strcmp(part(1:6),'penta6') % sometimes multiple element face_types are used; load PENTA6 elements too!
                faces_type = strcat(part);
                cornerpts = 6;
                nr_of_pentas    = fread(fid,                 1, '*uint32','l');
                geofilecontents.(part_descr).(['faces_' faces_type]) = fread(fid,[cornerpts nr_of_pentas],'*uint32','l')';
                
            elseif length(part) >= 5 && strcmp(part(1:5),'quad4') % sometimes multiple element face_types are used; load QUAD4 elements too!
                faces_type = strcat(part);
                cornerpts = 4;
                nr_of_quads    = fread(fid,                 1, '*uint32','l');
                geofilecontents.(part_descr).(['faces_' faces_type]) = fread(fid,[cornerpts nr_of_quads],'*uint32','l')';

            elseif isempty(part)
                [~] = fread(fid,1,'*char');
                if feof(fid)
                    break;
                else
                    error('Unexpected error; did not reach end of file after empty part');
                end

            else
                warning('This statement should not be reached. Most likely the GEO file was saved as multiple timesteps.. ignoring next timesteps.');
                break; % break out of while loop....
            end
        end
        fclose(fid); % close geo file
        
        %%% SOURCE OF GEO FILE FORMAT SPECIFIER BELOW
        %%% source: https://visualization.hpc.mil/wiki/EnSight_Gold_Geometry_Casefile_-_C_Binary_Example
        % All Data is plainly assumed to be exported as C binary Little Endian!
        %
        % C Binary                                            80 chars
        %+% ADDITIONAL_PARAMETER BEGIN TIME STEP                80 chars     % apparently added in fluent12 exported encas/geo files
        % description line 1                                  80 chars
        % description line 2                                  80 chars
        % node id <off/given/assign/ignore>                   80 chars
        % element id <off/given/assign/ignore>                80 chars
        %#% [extents                                            80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% xmin xmax ymin ymax zmin zmax]                       6 floats    % apparently skipped in fluent12 exported encas/ geo files
        
        % part                                                80 chars     % Multiple parts may exist within 1 file; loop over these until end of file
        % #                                                    1 int       % Index of file
        % description line                                    80 chars     % Name of current part
        % coordinates                                         80 chars
        % nn                                                   1 int       % Count of xyz coordinates
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        % x_n1 x_n2 ... x_nn                                  nn floats
        % y_n1 y_n2 ... y_nn                                  nn floats    % Loaded at once into vertices
        % z_n1 z_n2 ... z_nn                                  nn floats
        % element type                                        80 chars     % Determines np (nr of cornerpoints
        % ne                                                   1 int       % Count of faces
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        % n1_e1 n2_e1 ...                                     np_e1
        % n1_e2 n2_e2 ...                                     np_e2
        % .
        % .
        % n1_ne n2_ne ... np_ne                            ne*np ints
        % element type                                        80 chars     % apparently skipped in fluent12 exported encas/ geo files
        % .
        % .
        % part                                                80 chars     % Go BACK to PART above and repeat all elements
        
        %#% .                                                                % apparently skipped in fluent12 exported encas/ geo files
        %#% .                                                                % apparently skipped in fluent12 exported encas/ geo files
        %#% part                                                80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% #                                                    1 int       % apparently skipped in fluent12 exported encas/ geo files
        %#% description line                                    80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% block [iblanked] [with_ghost] [range]               80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% i j k # nn = i*j*k, ne = (i-1)*(j-1)*(k-1)           3 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [imin imax jmin jmax kmin kmax] # if range used:     6 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% nn = (imax-imin+1)* (jmax-jmin+1)* (kmax-kmin+1)                 % apparently skipped in fluent12 exported encas/ geo files
        %#% ne = (imax-imin)*(jmax-jmin)*(kmax-kmin)                         % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% x_n1 x_n2 ... x_nn                                  nn floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% y_n1 y_n2 ... y_nn                                  nn floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% z_n1 z_n2 ... z_nn                                  nn floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% [ib_n1 ib_n2 ... ib_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [ghost_flags]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [gf_e1 gf_e2 ... gf_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [node_ids]                                          80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [element_ids]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% part                                                80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% #                                                    1 int       % apparently skipped in fluent12 exported encas/ geo files
        %#% description line                                    80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% block rectilinear [iblanked] [with_ghost] [range]   80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% i j k # nn = i*j*k, ne = (i-1)*(j-1)*(k-1)           3 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [imin imax jmin jmax kmin kmax] # if range used:     6 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% nn = (imax-imin+1)* (jmax-jmin+1)* (kmax-kmin+1)                 % apparently skipped in fluent12 exported encas/ geo files
        %#% ne = (imax-imin)*(jmax-jmin)*(kmax-kmin)                         % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% x_1 x_2 ... x_i                                      i floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% y_1 y_2 ... y_j                                      j floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% z_1 z_2 ... z_k                                      k floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% [ib_n1 ib_n2 ... ib_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [ghost_flags]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [gf_e1 gf_e2 ... gf_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [node_ids]                                          80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [element_ids]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% part                                                80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% #                                                    1 int       % apparently skipped in fluent12 exported encas/ geo files
        %#% description line                                    80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% block uniform [iblanked] [with_ghost] [range]       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% i j k # nn = i*j*k, ne = (i-1)*(j-1)*(k-1)           3 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [imin imax jmin jmax kmin kmax] # if range used:     6 ints      % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% nn = (imax-imin+1)* (jmax-jmin+1)* (kmax-kmin+1)                 % apparently skipped in fluent12 exported encas/ geo files
        %#% ne = (imax-imin)*(jmax-jmin)*(kmax-kmin)                         % apparently skipped in fluent12 exported encas/ geo files
        %#%
        %#% x_origin y_origin z_origin                           3 floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% x_delta y_delta z_delta                              3 floats    % apparently skipped in fluent12 exported encas/ geo files
        %#% [ib_n1 ib_n2 ... ib_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [ghost_flags]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [gf_e1 gf_e2 ... gf_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [node_ids]                                          80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_n1 id_n2 ... id_nn]                             nn ints      % apparently skipped in fluent12 exported encas/ geo files
        %#% [element_ids]                                       80 chars     % apparently skipped in fluent12 exported encas/ geo files
        %#% [id_e1 id_e2 ... id_ne]                             ne ints      % apparently skipped in fluent12 exported encas/ geo files
    end                            % load .geo file

    function   scalarfilecontents = loadScalarFile(pathtofile,geometry_struct,istimeresolved_simulation)
        % load scalar file; see below for details on the file format
        version = ver('matlab'); 
        if str2double(version.Version) > 8
            narginchk(3,3); nargoutchk(1,1);
        end
        fid = fopen(pathtofile,'r');
        timepoint = 1;
        
        if istimeresolved_simulation == 1
            [first_line]                   = fread(fid,[1 80],'*char' ); % BEGIN TIME STEP % IGNORE
        end

        scalarfilecontents.description = fread(fid,[1 80],'*char' ); % description line 1          80 chars
        
        while ~feof(fid)
            part                       = fread(fid,[1 80],'*char');  % part                        80 chars
            
            if ~isempty(part) && strcmp(part(1:4),'part')
                
                %%% STEP 1: LOAD ALL DATA FOR THE CURRENT ELEMENT
                part_nr                = fread(fid,     1,'*uint32','l' );    % #                            1 int
                [~]                    = fread(fid,[1 80],'*char');
                [nrofelements_toread,elementdescr]    = getCountAndNameFromGeo( geometry_struct, part_nr );
                read_data              = fread(fid,nrofelements_toread,'*float','l' );
                
                %%% STEP 2: SAVE ALL DATA FOR THE CURRENT ELEMENT
                scalarfilecontents.(elementdescr).(['timestep' num2str(timepoint)]) = read_data;
                scalarfilecontents.(elementdescr).part_nr                           =   part_nr;
                
            elseif ~isempty(part) && strcmp(part(1:13),'END TIME STEP') % Check if NEXT TIME STEP IS STARTING OR IF FILE ENDS
                BEGIN_TIME_STEP = fread(fid,[1 80],'*char');      % Try reading a 80char line
                if feof(fid) % END OF FILE REACHED; STOP READING
                    break;
                elseif strcmp(BEGIN_TIME_STEP,first_line)   % NEW TIME STEP FOUND, restart while loop after reading the duplicate description
                    [~]       = fread(fid,[1 80],'*char');  % description line 1          80 chars % IGNORE description after saving 1st time
                    timepoint = timepoint + 1;              % Increase Time Step with 1
                end
                
            elseif isempty(part)
                [~] = fread(fid,1,'*char');
                if feof(fid)
                    break;
                else
                    error('Unexpected error; did not reach end of file after empty part')
                end
                
            else % code should not arrive here... throw error!
                disp(['lastfile: ' pathtofile ])
                error(['Unexpected error; Last element read: ' num2str(part_nr) ' (timepoint: ' num2str(timepoint) ')'])
                
            end % if
        end % while
        
        scalarfilecontents.timesteps = timepoint; % set number of loaded timepoints
        fclose(fid); % close file
        
        
        %%% SCALAR - from: http://www-vis.lbl.gov/NERSC/Software/ensight/docs/OnlineHelp/UM-C11.pdf
        %+% ADDED: BEGIN TIME STEP    80 chars %%% ADDED
        % description line 1          80 chars
        % part                        80 chars
        % #                            1 int
        % coordinates                 80 chars
        % s_n1 s_n2 ... s_nn          nn floats
        % part                        80 chars
        % .
        % .
        % part                        80 chars
        % #                            1 int
        % block # nn = i*j*k          80 chars
        % s_n1 s_n2 ... s_nn          nn floats
        
    end         % load .scl# file
   
    function   vectorfilecontents = loadVectorFile( pathtofile, geometry_struct, istimeresolved_simulation )
        % load scalar file; see below for details on the file format
        version = ver('matlab');
        if str2double(version.Version) > 8
            narginchk(3,3); nargoutchk(1,1);
        end
        
        fid = fopen(pathtofile,'r');
        timepoint = 1;
        
        if istimeresolved_simulation == 1
            [first_line]                   = fread(fid,[1 80],'*char' ); % BEGIN TIME STEP % IGNORE
        end
        vectorfilecontents.description = fread(fid,[1 80],'*char' ); % description line 1          80 chars
        
        while ~feof(fid)
            part                       = fread(fid,[1 80],'*char');  % part                        80 chars
            
            if ~isempty(part) && strcmp(part(1:4),'part')
                
                %%% STEP 1: LOAD ALL DATA FOR THE CURRENT ELEMENT
                part_nr                = fread(fid,     1,'*uint32','l' );    % #                            1 int
                [~]                    = fread(fid,[1 80],'*char');
                [nrofelements_toread,elementdescr]    = getCountAndNameFromGeo( geometry_struct, part_nr );
                read_data              = fread(fid,[nrofelements_toread 3],'*float','l' );
                
                %%% STEP 2: SAVE ALL DATA FOR THE CURRENT ELEMENT
                vectorfilecontents.(elementdescr).(['timestep' num2str(timepoint)]) = read_data;
                vectorfilecontents.(elementdescr).part_nr                           =   part_nr;
                
            elseif ~isempty(part) && strcmp(part(1:13),'END TIME STEP') % Check if NEXT TIME STEP IS STARTING OR IF FILE ENDS
                BEGIN_TIME_STEP = fread(fid,[1 80],'*char');      % Try reading a 80char line
                if feof(fid) % END OF FILE REACHED; STOP READING
                    break;
                elseif strcmp(BEGIN_TIME_STEP,first_line)   % NEW TIME STEP FOUND, restart while loop after reading the duplicate description
                    [~]       = fread(fid,[1 80],'*char');  % description line 1          80 chars % IGNORE description after saving 1st time
                    timepoint = timepoint + 1;              % Increase Time Step with 1
                end
            elseif isempty(part)
                [~] = fread(fid,1,'*char');
                if feof(fid)
                    break;
                else
                    error('Unexpected error; did not reach end of file after empty part');
                end

            else % code should not arrive here... throw error!
                error(['Unexpected error; Last element read: ' num2str(part_nr) ' (timepoint: ' num2str(timepoint) ')'])
                
            end % if
        end % while
        
        vectorfilecontents.timesteps = timepoint; % set number of loaded timepoints
        fclose(fid); % close file
        
        
        %%% VECTOR - from: http://www-vis.lbl.gov/NERSC/Software/ensight/docs/OnlineHelp/UM-C11.pdf
        % description line 1            80 chars
        % part                          80 chars
        % #                              1 int
        % coordinates                   80 chars
        % vx_n1 vx_n2 ... vx_nn         nn floats
        % vy_n1 vy_n2 ... vy_nn         nn floats
        % vz_n1 vz_n2 ... vz_nn         nn floats
        % part                          80 chars
        % .
        % .
        % part                          80 chars
        % #                              1 int
        % block # nn = i*j*k            80 chars
        % vx_n1 vx_n2 ... vx_nn         nn floats
        % vy_n1 vy_n2 ... vy_nn         nn floats
        % vz_n1 vz_n2 ... vz_nn         nn floats
        
    end      % load .vel file

    function    otherfilecontents = load_OTHER_file(pathtofile) %#ok<INUSD>
        otherfilecontents = 'UNSUPPORTED; please read, copy and edit loadVectorFile or loadScalarFile function in loadEncasFile.m to support this filetype.';
        warning('TENSOR, TENSOR9, COMPLEX SCALAR and COMPLEX VECTOR files are currently unsupported.')
        %%% TENSOR
        % description line 1            80 chars
        % part                          80 chars
        % #                              1 int
        % coordinates 80 chars
        % v11_n1 v11_n2 ... v11_nn      nn floats
        % v22_n1 v22_n2 ... v22_nn      nn floats
        % v33_n1 v33_n2 ... v33_nn      nn floats
        % v12_n1 v12_n2 ... v12_nn      nn floats
        % v13_n1 v13_n2 ... v13_nn      nn floats
        % v23_n1 v23_n2 ... v23_nn      nn floats
        % part                          80 chars
        % .
        % .
        % part                          80 chars
        % # 1 int
        % block     # nn = i*j*k        80 chars
        % v11_n1 v11_n2 ... v11_nn      nn floats
        % v22_n1 v22_n2 ... v22_nn      nn floats
        % v33_n1 v33_n2 ... v33_nn      nn floats
        % v12_n1 v12_n2 ... v12_nn      nn floats
        % v13_n1 v13_n2 ... v13_nn      nn floats
        % v23_n1 v23_n2 ... v23_nn      nn floats
        
        %%% TENSOR9
        % description line 1                    80 chars
        % part                                  80 chars
        % #                                      1 int
        % coordinates                           80 chars
        % v11_n1 v11_n2 ... v11_nn              nn floats
        % v12_n1 v12_n2 ... v12_nn              nn floats
        % v13_n1 v13_n2 ... v13_nn              nn floats
        % v21_n1 v21_n2 ... v21_nn              nn floats
        % v22_n1 v22_n2 ... v22_nn              nn floats
        % v23_n1 v23_n2 ... v23_nn              nn floats
        % v31_n1 v31_n2 ... v31_nn              nn floats
        % v32_n1 v32_n2 ... v32_nn              nn floats
        % v33_n1 v33_n2 ... v33_nn              nn floats
        % part                                  80 chars
        % .
        % .
        % part                                  80 chars
        % #                                      1 int
        % block # nn = i*j*k                    80 chars
        % v11_n1 v11_n2 ... v11_nn              nn floats
        % v12_n1 v12_n2 ... v12_nn              nn floats
        % v13_n1 v13_n2 ... v13_nn              nn floats
        % v21_n1 v21_n2 ... v21_nn              nn floats
        % v22_n1 v22_n2 ... v22_nn              nn floats
        % v23_n1 v23_n2 ... v23_nn              nn floats
        % v21_n1 v21_n2 ... v21_nn              nn floats
        % v22_n1 v22_n2 ... v22_nn              nn floats
        % v23_n1 v23_n2 ... v23_nn              nn floats
        
        %%% COMPLEX SCALAR
        % description line 1              80 chars
        % part                            80 chars
        % #                                1 int
        % coordinates                     80 chars
        % s_n1 s_n2 ... s_nn              nn floats
        % part                            80 chars
        % .
        % .
        % part                            80 chars
        % #                                1 int
        % block # nn = i*j*k              80 chars
        % s_n1 s_n2 ... s_nn              nn floats
        
        %%% COMPLEX VECTOR
        % description line 1                       80 chars
        % part                                     80 chars
        % #                                         1 int
        % coordinates 80 chars
        % vx_n1 vx_n2 ... vx_nn                    nn floats
        % vy_n1 vy_n2 ... vy_nn                    nn floats
        % vz_n1 vz_n2 ... vz_nn                    nn floats
        % part                                     80 chars
        % .
        % .
        % part                                     80 chars
        % #                                         1 int
        % block # nn = i*j*k                       80 chars
        % vx_n1 vx_n2 ... vx_nn                    nn floats
        % vy_n1 vy_n2 ... vy_nn                    nn floats
        % vz_n1 vz_n2 ... vz_nn                    nn floats
    end            % not yet implemented

    function [count,elementdescr] = getCountAndNameFromGeo( geometry_struct, part_nr )
        % find how much elements the vertices of element #part_nr have in the geometry
        for fname = fieldnames(geometry_struct)'
            if ~any(strcmpi({'type','description','nr_parts'},fname))
                if (geometry_struct.(fname{1}).part_nr == part_nr)
                    count = size(geometry_struct.(fname{1}).vertices,1);
                    elementdescr = fname{1};
                    return;
                end
            end
        end
    end % getCountAndNameFromGeo

    function checkTimings(encas_file_struct)
        try
            for fnm = fieldnames(encas_file_struct)'
                if any(strcmp(fnm,{'velocity','total_pressure','velocity_magnitude','wall_shear'}))
                    fieldname_thingy = fnm{1};
                    break; % exit for loop
                end
            end
            
            for fnm = fieldnames(encas_file_struct.(fieldname_thingy))'
                if ~any(strcmp(fnm,{'description','timesteps','type'}))
                    fieldname_thingy2 = fnm{1};
                end
            end
            
            if isempty(encas_file_struct.encas_file.timings) && length(fieldnames(encas_file_struct.(fieldname_thingy).(fieldname_thingy2)))-1 == 1
                % correct; no timesteps, thus 1 timestep!
                
            elseif length(encas_file_struct.encas_file.timings)  == length(fieldnames(encas_file_struct.(fieldname_thingy).(fieldname_thingy2)))-1 %-1 to get rid of fieldname part_nr
                % everything looks ok; continue
                
            else
                % not 100% ok; throw non-fatal warning to warn user
                warning(['The number of timesteps described in the encas file (' num2str(length(encas_file_struct.encas_file.timings)) ') does not match the number of timesteps in the datafiles (' num2str(length(fieldnames(encas_file_struct.x_velocity.surface3))-1) ')']);
            end
        catch err
            % THIS IS WEIRD THROW ERROR IN Warning message
             warning(['checktimings failed; but is not essential. (error: ' err.message ')']);
        end
    end %check timings

    function istimeresolved_simulation = checkTimeResolvedSimulation(encas_inputfile)
        fid = fopen(encas_inputfile);
        istimeresolved_simulation = 0;
        while ~feof(fid) && (istimeresolved_simulation == 0)
            istimeresolved_simulation = strcmp(fgetl(fid),'TIME');
        end
    end
end % main function