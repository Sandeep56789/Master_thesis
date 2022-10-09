%======================================================================
% source
%======================================================================
source_x = [1900 900 5]; % source initial position (x,y,z) (m)
source_v = [0 0 0];   % source velocity (vx,vy,vz) (m/s)
source_nrays =  2001;   % number of propagation rays considered
source_aperture = 75;   % maximum launching angle (degrees)
source_ray_step = 10;   % ray calculation step (meters), modified later 

%======================================================================
% transmitted signal
%======================================================================
fc   = 25600;  % input signal carrier           frequency (Hz)
fs_x = 10000;  % input baseband signal sampling frequency (Hz)

%======================================================================
% Bottom
%=====================================================================
bottom_properties  = [1465 1.5 0.06];
%bottom_properties(1) = compressional speed (m/s)
%bottom_properties(2) = bottom density      (g/cm3)
%bottom_properties(3) = bottom attenuation  (dB/wavelength)

%======================================================================
% Array
%======================================================================
array_x = [2800 400 0]; % receivers array initial position (x,y,z) (m)
array_v = [0 0 0];      % receiver array velocity (vx,vy,vz) (m/s)
first_hyd = 30;         % depth of the first hydrophone  (m)
 last_hyd = 60;         % depth of the last hydrophone   (m)
delta_hyd =  2;         % separation between hydrophones (m)

%======================================================================
% Wind induced sea surface wave 
%======================================================================
U         =       0; % wind speed 
theta     =       0; % direction of propagation in degrees
spreading =  'none'; % options: 'none', 'mits', 'hass'
