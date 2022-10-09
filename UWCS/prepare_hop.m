% Prepare the run of hop with from top to bottom: 
% 
% - source
% - surface 
% - sound speed profile 
% - bottom 
% - array 

%==================================================================
%  
%  Source properties 
%  
%==================================================================
zbox = 1.5*max( zbty(:) );
rbox = 1.5*source_array_range/1000;
source_data.box      = [source_ray_step zbox rbox];
source_data.f        = fc;
source_data.thetas   = [source_nrays -source_aperture source_aperture];
source_data.p        = [];
source_data.comp     = '';
source_data.ds       = source_array_range/1500; % Adapt ray step according to source/array distance

%==================================================================
%  
%  Surface properties 
%  
%==================================================================
surface_data.itype = '''L''';      % Sea surface interpolation type
surface_data.x     = [];           % Surface coordinates 
surface_data.p     = [];           % Surface properties 
%==================================================================
%  
%  Sound speed profile properties
%  
%==================================================================
ssp_data.r   = [];
ssp_data.z   = z;
ssp_data.c   = c;
%==================================================================
%  
%  Bottom properties 
%  
%==================================================================
nbty = 101; % number of spatial samples along one cut %
bottom_data.itype = '''L'''; % Bottom interpolation type
bottom_data.x     = [];         % Bottom coordinates 
bottom_data.p     = [bottom_properties(1) 0.0 bottom_properties(2:3)];     % Bottom properties 
bottom_data.aunits = '''W''' ; % Bottom absorption units: (dB/m)*kHZ 
%==================================================================
%  
%  Array properties 
%  
%==================================================================
 narr0 = 50; 
 if flag_flat == 1 
 options.options1    = '''SVWT'''; 
 else
 options.options1    = '''SVWT*'''; 
 end 
 options.options2    = '''A*'''; 
 options.options3    = '''A''';
 options.options4    = []; 
