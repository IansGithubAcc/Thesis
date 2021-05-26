clear; close all; clc;

addpath(genpath('../utilities/'))
addpath(genpath('../datasets/'))
addpath(genpath('../utilities/m_map/'))
addpath(genpath('../m_map/'))

%% STEP 1: set mesh extents and set parameters for mesh.
bbox = [0.06 0.42;		% lon_min lon_max
        50.69490 50.92]; 		% lat_min lat_max 50.56
min_el    = 20;  		% minimum resolution in meters.
max_el    = 100; 		% maximum resolution in meters. 
max_el_ns = 100;        % maximum resolution nearshore in meters.
grade     = 0.35; 		% mesh grade in decimal percent.
R         = 100;    		% number of elements to resolve feature width.
%% STEP 2: specify geographical datasets and process the geographical data 
%% to be used later with other OceanMesh classes...
dem       = 'DEFRA_DEM_2016_1_Second_South.nc';
coastline = 'GSHHS_f_L1';
gdat = geodata('shp',coastline,'dem',dem,'bbox',bbox,'h0',min_el);
plot(gdat,'dem')
% data = geodata('DEM','DEFRA_DEM_2016_1_Second_South.nc','h0',100);
%% STEP 3: create an edge function class
fh = edgefx('geodata',gdat,...
            'fs',R,'max_el_ns',max_el_ns,...
            'max_el',max_el,'g',grade);
%% STEP 4: Pass your edgefx class object along with some meshing options and
% build the mesh...
mshopts = meshgen('ef',fh,'bou',gdat,'plot_on',1,'nscreen',5,'proj','trans');
mshopts = mshopts.build; 

%% STEP 5: Plot it and write a triangulation fort.14 compliant file to disk.
% Get out the msh class and put on nodestrings
m = mshopts.grd;
m = make_bc(m,'auto',gdat,'distance'); % make the boundary conditions
m = interp(m,dem,'type','depth');
m.b = m.b +3.09115637461;
m.b(isnan(m.b))=0;
plot(m,'bd',0);

% if you want to write into fort.14...
% write(m,'../../../Schenarios/Scenario_1/SWAN/Input/fort','14');