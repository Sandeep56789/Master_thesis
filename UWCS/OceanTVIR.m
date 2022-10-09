function OceanTVIR(Scenario_info)

%  
% OceanTVIR - Generate Ocean Time-Varying Impulse Responses
%
% the results are saved in 'Results.mat'
%

%====================================================================== 
% load the waveguide, source and array properties
%======================================================================

xati = Scenario_info.xati;
yati = Scenario_info.yati;
zati = Scenario_info.zati;
ru   = Scenario_info.ru; nati = length( ru );
hu   = Scenario_info.hu;
xu   = Scenario_info.xu;
yu   = Scenario_info.yu;
Rs    = Scenario_info.Rs;
c    = Scenario_info.c;
z    = Scenario_info.z;
xbty = Scenario_info.xbty;
ybty = Scenario_info.ybty;
zbty = Scenario_info.zbty;
wind = Scenario_info.wind;
tmax = Scenario_info.tmax;
cg   = Scenario_info.cg;

flags= Scenario_info.flags; flag_flat = flags(1); flag_mosaic = flags(2); 
bottom_properties = Scenario_info.bottomp ;
rangemax          = Scenario_info.rangemax;

source_nrays      = Scenario_info.source_nrays;
source_aperture   = Scenario_info.source_aperture;
source_ray_step   = Scenario_info.source_ray_step;
source_position   = Scenario_info.source_position;
source_array_range= Scenario_info.source_array_range;
receiver_position = Scenario_info.receiver_position;
hydrophones       = Scenario_info.hydrophones;
nhyd =3% length(hydrophones);

tmax = Scenario_info.tmax;
fs_x = Scenario_info.fs_x ;
fc   = Scenario_info.fc;

clear Scenario_info

%======================================================================
% complete the waveguide definition and prepare hop estimation
%======================================================================
prepare_hop;

%======================================================================
%                                              
% run Bellhop and evaluate the impulse response
%
%======================================================================
% external parameters (to change if needed)
%======================================================================
criterium_threshold = 0.1;      % threshold applied to the scattering test %
%criterium_threshold = 0.5;     % threshold applied to the scattering test %
%fs_time = 5;                   % initial sampling frequency along time for the channel
fs_time = 10;                   % initial sampling frequency along time for the channel
securitas = 0.01;               % interval of time kept before the first and after the last arrival (in seconds) %
plotSC = 1;                     % eventual plot of the doppler content of the IRs for each iteration %

%======================================================================
% loop parameters initialization
%======================================================================
t = 0:1/fs_time:tmax;   % the first realization time axis %
ntimes = length(t);
ntimes=2;
h = []; save('Results','h','-APPEND'); clear h;
vt_h = [];
scatterok = 0;
iteration = 0;

%======================================================================
% deal with figures
%======================================================================
%while((~scatterok)&&(iteration<10))
% this will be up-to and including iteration 4
while((~scatterok)&&(iteration<4))
    iteration = iteration+1;
    disp(['iteration ',num2str(iteration),', fs_time = ',num2str(fs_time),', ntimes = ',num2str(ntimes)]);
    fileadvance = fopen('computation_advance.txt','a+');
    logrecord = strcat('iteration ',num2str(iteration),', fs_time = ',num2str(fs_time),', ntimes = ',num2str(ntimes));
    fprintf(fileadvance,logrecord);
    fclose(fileadvance);
    tic;
    % IR calculation %
    %================%
    if(iteration==1) % evaluate first the arrival times for optimization of the size of the impulse response %
        narr = source_nrays;
        min_travel_time = zeros(1,ntimes);
        max_travel_time = zeros(1,ntimes);
        store_travel_time = zeros(nhyd,5,ntimes);
        store_amplitude = zeros(nhyd,5,ntimes);
        store_phi = zeros(nhyd,5,ntimes);
        for tt = 1:ntimes
            run_hop;
            store_travel_time(:,:,tt) = travel_time;
            store_amplitude(:,:,tt) = amplitude;
            store_phi(:,:,tt) = phi;
            indexesvalid = (travel_time > 0);
            min_travel_time(tt) = min(travel_time(indexesvalid));
            max_travel_time(tt) = max(travel_time(indexesvalid));
        end
        tstart = min(min_travel_time) - securitas;
        tend = max(max_travel_time) + securitas;
        vtau = tstart:1/fs_x:tend;
        ntau = length(vtau);
        if(rem(ntau,2)) % force ntau to be even %
            vtau = vtau(1:end-1);
            ntau = ntau-1;
        end
        save('Results','vtau','-APPEND'); clear vtau;
        hit = zeros(ntimes,ntau,nhyd);
        for tt = 1:ntimes
            travel_time = squeeze(store_travel_time(:,:,tt));
            amplitude = squeeze(store_amplitude(:,:,tt));
            phi = squeeze(store_phi(:,:,tt));
            SerialiseIR;
        end
%         narr = max(narr0,max(narr));
        clear min_travel_time max_travel_time store_travel_time store_amplitude store_phi indexesvalid;
    else
        for tt = 1:ntimes
            run_hop;
            SerialiseIR;
        end
    end
    clear ranges_array depths_array travel_time amplitude phi
    
    % save the results %
    %==================%
    vt_h = [vt_h,t]; clear t;
    [vt_h,indices] = sort(vt_h,'ascend');
    load('Results','h');
    h = [h;hit]; clear hit;
%     h = h(indices,:,:);
    save('Results','h','-APPEND');

    hsize = prod(size(h));

    % check the scattering function %
    %===============================%
    SCsum = sum(sum(abs(fftshift(fft(h,size(h,1),1),1)).^2,3),2); clear h;
    SCsum = SCsum/mean(SCsum);
    criterium = mean(SCsum([1:((ntimes+rem(ntimes,2))/2),end-(ntimes+rem(ntimes,2))/2+1:end]));
    scatterok = (criterium<criterium_threshold);
     colormappe = colormap;
    if(plotSC)
        vd = [-(length(vt_h)-1)/2:(length(vt_h)-1)/2]/length(vt_h)*fs_time;
        %figure(13), hold on; plot(SCsum,'Color',colormappe(5*iteration,:),'linewidth',3);
       
%xlabel('doppler frequency (Hz)');
%ylabel('amplitude');
%title('UWOC scattering amplitude')
%grid on;
    end
    disp(['criterium = ',num2str(criterium)]);
    fileadvance = fopen('computation_advance.txt','a+');
    logrecord = strcat(', criterium = ',num2str(criterium), ', size of H = ', num2str(hsize),' \n ');
    fprintf(fileadvance,logrecord);
    
    % initialise the realization time for the next iteration %
    %========================================================%
    if(~scatterok)
        fs_time = 2*fs_time;
        tnew = 0:1/fs_time:tmax;
        t = tnew(2:2:end);
        ntimes = length(t);
        hit = zeros(ntimes,ntau,nhyd);
        clear tnew
    end
    toc;
end
save('Results','vt_h','-APPEND');

fileadvance = fopen('computation_advance.txt','a+');
logrecord = strcat('evaluation of impulse responses completed \n ');
fprintf(fileadvance,logrecord);
fclose(fileadvance);
