% get the current time index %
indicetime = round(t(tt)*fs_x+1);

% get source position %
source_x = source_position(1,indicetime);
source_y = source_position(2,indicetime);
source_z = source_position(3,indicetime);
source_data.position = source_z;
source_data.c = interp1(ssp_data.z,ssp_data.c,source_z);
source_data.source_position=source_position;
% get receiver position %
receiver_x = receiver_position(1,indicetime);
receiver_y = receiver_position(2,indicetime);
receiver_z = hydrophones + receiver_position(3,indicetime);
source_data.receiver_position=receiver_position;
source_data.source_nrays=source_nrays;
% define source-receiver range %
dx = receiver_x - source_x;
dy = receiver_y - source_y;
rangemax = sqrt( dx*dx + dy*dy );

% define transect coordinates %
if flag_flat == 0
   if flag_mosaic == 1 
   np = fix(rangemax/Rs);
   else
   np = nati;
   end
else 
   np = 2; 
end 
xtransect = linspace(source_x,receiver_x,np);
ytransect = linspace(source_y,receiver_y,np);
rtransect = sqrt( ( xtransect - source_x ).^2 + ( ytransect - source_y ).^2 );
% figure(5)
% mesh(xu,yu,hu), hold on
% plot3(xtransect,ytransect,zeros(size(xtransect)))
% plot3(xtransect,ytransect,zeros(size(xtransect)),'o')
% hold off 
% get surface elevation along the transect %
if flag_flat == 0
  if flag_mosaic == 0 % Temporal shift of sea surface:
     hshifted  = shiftsurf(xu,yu,hu,wind.thetaU,cg,t(tt),tmax,1,1,1,1);
     htransect = interp2(xu,yu,hshifted,xtransect,ytransect);
  else                % Temporal and space shifts of sea surface:
    [htransect,rtransect] = mosaici(xu,yu,hu,xtransect,ytransect,rangemax,rangemax,wind.theta,wind.v,t(tt),tmax);
  end 
else 
  htransect = [0 0];
end 

% Non-flat surface: be sure altimetry starts from min to max at range 0: 

if flag_flat == 0 
altimetry(1,:) = rtransect;
altimetry(2,:) = htransect;

altimetry(1,:) = altimetry(1,:) - min(altimetry(1,:));
[altimetry(1,:),indexes] = sort(altimetry(1,:)); 
 altimetry(2,:) = altimetry(2,indexes);
 altimetry(1,:) = altimetry(1,:)/1000; % Bellhop expects this in km...
 surface_data.x = altimetry;
end 

% define bottom bathymetry %
xbtyi = linspace(source_x,receiver_x,nbty);
ybtyi = linspace(source_y,receiver_y,nbty);
bathymetry(1,:) = linspace(0,1.1*rangemax,nbty);
bathymetry(2,:) = interp2(xbty,ybty,zbty,xbtyi,ybtyi);
%{
figure(12),plot(bathymetry(1,:),bathymetry(2,:),'linewidth',2), view(0,-90);grid on
ylabel('Travelling Distance of nodes');
xlabel('Number of nodes');
title('Tracking of motion from Tx to Rx in UWOC');
%}
bathymetry(1,:) = bathymetry(1,:)/1000;
bottom_data.x   = bathymetry;

options.r = rangemax/1000; options.z = receiver_z;


amp = (receiver_z(1:nhyd)') * exp( 1i *(linspace(1,fc,5) )* pi/180);
amplitude=abs(real(amp));
dist=((source_position(1,1:5)-receiver_position(1,1:5)).^2-(source_position(2,1:5)-receiver_position(2,1:5)).^2);          
dg=300e3;
travel_time=hydrophones(1:nhyd)'*dist/dg;
% format the output data %
phi         = angle(squeeze(amp));

