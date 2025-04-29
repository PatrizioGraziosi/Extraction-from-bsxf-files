% Written and developed by                                                %
% Patrizio Graziosi, patrizio.graziosi@gmail.com, during the              %  
% Marie Curie - Individual Fellowships  GENESIS - project ID 788465       %
% Generic transport simulator for new generation thermoelectric materials %
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License. See the file `LICENSE' in  the root directory   %
% of the present distribution.                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source when publishing results obtained            %
% using the present code                                                  %
%                                                                         %
% ----------------------------------------------------------------------- %

function bxsf_to_ELECTRA_2D(fileName,material_name,alat,reduce_flag,interpolation_factor) %codegen
if strcmp(reduce_flag,'y')
    reduce_SOC = 'SOC_not_magn';
else
    reduce_SOC = 'no';
end
if interpolation_factor > 1
    bands_interpolation =   'yes'; % if you wants numerical bands interpolation, it uses griddata for sparse points
else
    bands_interpolation =   'no'; % if you wants numerical bands interpolation, it uses griddata for sparse points
end
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


% extraction codes run

% this sub-function takes the data from the bxsf file and compose a 4D matrix
[points_in_axis_kx, points_in_axis_ky , num_of_bands,...
BandMatrix, Fermi, a, b, extra_point_flag, points_in_axis_kz ]...
= Taking_data_from_bxsf(fileName);

% now we compose the matrixes of the coordinates in k-space

Ek = zeros(points_in_axis_kx,points_in_axis_ky,num_of_bands);

if strcmp(extra_point_flag,'no')

for id_band = 1:num_of_bands %#ok<*FXUP>
    id_k=1;
    for id_x=1:points_in_axis_kx
        for id_y=1:points_in_axis_ky
                Ek(id_x,id_y,id_band) = BandMatrix(id_k,id_band)-Fermi;                
                id_k=id_k+1;
        end
    end
end

elseif strcmp(extra_point_flag,'yes')

Ek_temp = zeros(points_in_axis_kx,points_in_axis_ky,points_in_axis_kz);
for id_band = 1:num_of_bands %#ok<*FXUP>
    id_k=1;
    for id_x=1:points_in_axis_kx
        for id_y=1:points_in_axis_ky
            for id_z=1:points_in_axis_kz
                Ek_temp(id_x,id_y,id_z) = BandMatrix(id_k,id_band)-Fermi;                
                id_k=id_k+1;
            end
        end
    end
    Ek(:,:,id_band) = Ek_temp(:,:,1);
end

end

% alat shall be inputted      
    blat = alat;

B_matrix = [ a*2*pi/(alat*1e-9) ; b*2*pi/(blat*1e-9) ] ;

kx_matrix = zeros(points_in_axis_kx, points_in_axis_ky);
ky_matrix = kx_matrix; 
    for id_x = (points_in_axis_kx - 1) : -1 : 0
        for id_y = (points_in_axis_ky - 1) : -1 : 0 % 0 : points_in_axis_ky - 1
                
                k_vector_not_norm = [id_x id_y ]*B_matrix; 
                
                kx_matrix(id_x+1,id_y+1) = 1/points_in_axis_kx * k_vector_not_norm(1);
                ky_matrix(id_x+1,id_y+1) = 1/points_in_axis_ky * k_vector_not_norm(2);
                
        end
    end
    

% additional features
if strcmp(reduce_SOC,'SOC_not_magn')
    Ek_full = Ek;
    Ek=zeros(points_in_axis_kx,points_in_axis_ky,num_of_bands/2);
    for i = 1:num_of_bands/2
        Ek(:,:,i) = Ek_full(:,:,2*i-1);
    end    
end
if strcmp(bands_interpolation,'yes')
    [ Ek, kx_matrix, ky_matrix ] = Interpolation_2D(Ek,kx_matrix,ky_matrix,a,b,alat,interpolation_factor) ;
end

save_filename=['Ek_',material_name,'.mat'];
save(save_filename, 'Ek', 'kx_matrix', 'ky_matrix', 'a', 'b', 'alat', 'save_filename')  


% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

function [points_in_axis_kx, points_in_axis_ky , num_of_bands,...
BandMatrix, Fermi, a, b, extra_point_flag, points_in_axis_kz] ...
= Taking_data_from_bxsf(fileName)

inputFile = strcat(fileName);
fid = fopen(inputFile);
% searching for the point where the bands start
for i = 1:60
    temp = fgetl(fid);
    k = strfind(temp,'Fermi Energy');
    if k >= 1
        numstr = regexp(temp,'(-)?\d+(\.\d+)?(e(-|+)\d+)?','match') ;
        Fermi = str2double(numstr) ; 
    end
    k = strfind(temp,'BEGIN_BLOCK_BANDGRID');
    if k >= 1
        num_lines_to_skip = i+9;
        break % i is kept
    end
end
% skip two lines
for i=1:2
    fgetl(fid);
end

% number of bands in the file
temp = fgetl(fid);
num_of_bands = str2num(temp) ; 

% number of points in the k space, that corresponds to the number of values in each band 
temp = fgetl(fid);
points_array = str2num(temp) ; 
points_in_axis_kx = points_array(1);
points_in_axis_ky = points_array(2);
if max(size(points_array)) == 3
    points_in_axis_kz = points_array(3);
    num_of_points = points_in_axis_kx * points_in_axis_ky * points_in_axis_kz;
    extra_point_flag = 'yes';
else
    num_of_points = points_in_axis_kx * points_in_axis_ky ;
    extra_point_flag = 'no';
    points_in_axis_kz = 1;
end


% the a, b, c vectors that ate in from 3 to 1 lines above the number of
% line to skip, which indiactes the origin of the grid, usually Gamma
fgetl(fid);
% vectors
temp = fgetl(fid); 
a = str2num(temp); %#ok<*NASGU>
temp = fgetl(fid); 
b = str2num(temp); 
temp = fgetl(fid); 
c = str2num(temp); 


fgetl(fid); % thus we arrive to the num_line_to_skip, extraction starts
temp = fgetl(fid);
first_eigenvalues = str2num(temp);
num_elements_each_row = size(first_eigenvalues,2) ;

if rem(num_of_points,num_elements_each_row) ~= 0
    num_lines_to_read = floor(num_of_points/num_elements_each_row)+1;
elseif rem(num_of_points,num_elements_each_row) == 0
    num_lines_to_read = floor(num_of_points/num_elements_each_row);
end

fid = fopen(inputFile);
for i = 1:num_lines_to_skip
    fgetl(fid);
end
for id_band = 1:num_of_bands
    if id_band == 1
        starting_line = num_lines_to_skip+1;
    else
        starting_line = num_lines_to_skip+(num_lines_to_read+1)*(id_band-1)+1;
    end
    end_line = starting_line+num_lines_to_read-1;

    
    Band_temp=[];
    for id_line = starting_line:end_line
        
        temp = fgetl(fid);
        
        Band_temp = [Band_temp str2num(temp)]; %#ok<*AGROW>

    end
    BandMatrix(:,id_band) = (Band_temp(:)); %#ok<*ST2NM>
    num_lines_inbetween = 1; % number of line in between two bands, 1 is the typical values in bsxf files
    % skip the lines
    for i=1:num_lines_inbetween
        temp=fgetl(fid);
    end
    
end

end

    function [ Ek_i, kx_matrix_i, ky_matrix_i ] = ...
            Interpolation_2D(Ek,kx_matrix,ky_matrix,a,b,alat,interpolation_factor)
        nkx = size(Ek,1);
        nky = size(Ek,2);
        n_bands = size(Ek,3);

        if isempty(gcp('nocreate')) == 0
            delete(gcp('nocreate'))
        end

        nlc = feature('NumCores') ;
        if nlc > n_bands
            nlc = n_bands;
        end
        if nlc > 3
            nlc = nlc-2; % to keep the computer alive
        elseif nlc > 1
            nl = nlc-1;
        end
        ppp = parpool('local',nlc);
        
        if exist('blat','var') == 0 
            blat = alat;
        end
        if exist('nk_new_x','var') == 0 
            nk_new_x = floor( nkx*interpolation_factor ) + 1;
        end
        if exist('nk_new_y','var') == 0 
            nk_new_y = floor( nky*interpolation_factor ) + 1;
        end
        B_matrix = [ a*2*pi/(alat*1e-9) ; b*2*pi/(blat*1e-9) ] ;
        for id_y=0:nk_new_y-1
            for id_x=0:nk_new_x-1

                k_vector_not_norm = [id_x id_y ]*B_matrix; 

                Kx_interp(id_x+1,id_y+1) = 1/nk_new_x * k_vector_not_norm(1);
                Ky_interp(id_x+1,id_y+1) = 1/nk_new_y * k_vector_not_norm(2);

            end
        end
        
        Ek_m=ones(nk_new_x,nk_new_y,n_bands);
        
        parfor i = 1:n_bands

            Ek_temp = Ek(:,:,i);
            Ek_m_temp = griddata(kx_matrix,ky_matrix,Ek_temp,Kx_interp,Ky_interp,'nearest');

            n = isnan(Ek_m_temp);
            n_pos = find(n);

            kx_n = Kx_interp(n_pos)';
            ky_n = Ky_interp(n_pos)';

            % ptCloud = pointCloud( [Kx_interp(:),Ky_interp(:), 0*Ky_interp(:)] ) ;

            for i_p = 1:size(kx_n,2)

                dist_array = sqrt( (Kx_interp-kx_n(i_p)).^2 + (Ky_interp-ky_n(i_p)).^2 + (Kz_interp-kz_n(i_p)).^2 );

                [~,indexes] = mink(dist_array,20) ;

                % [indexes,~] = findNearestNeighbors(ptCloud,[kx_n(i_p),ky_n(i_p),0],6) ;

                nn = abs( isnan(Ek_m_temp(indexes)) -1) ;
                nn_pos = nn;
                indexes = double(indexes);

                E_nn = Ek_m_temp(nonzeros(indexes.*nn_pos))';        
                Ek_m_temp(n_pos(i_p)) = mean(E_nn); %#ok<*UDIM>

            end

            Ek_m(:,:,i) = Ek_m_temp;
        end
        
        delete(gcp('nocreate'))

        Ek_i = Ek_m;
        kx_matrix_i=Kx_interp; ky_matrix_i=Ky_interp;

    end



end