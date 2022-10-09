clc
clear;
close all
warning off all
NrRepetitions = 1000;%24;                     % Number of Monte Carlo repetitions
M_SNR_dB_OFDM = -5:5:30;                % Signal-to-noise ratio

PerfectChannelKnowledge = false;          % Perfect channel information or estimation

ModulationOrder = 4;                 % Modulation order, 4,16,64,...

L = 12;                                   % Number of subcarriers
K = 2^5;                                  % Number of symbols per Block = Spreading Length
Kall = K*3+3;                             % Total Number of symbols (3 guard symbols!)

SubcarrierSpacing = 15e3;                 % Subcarrier spacing (15kHz, same as LTE)
SamplingRate =  SubcarrierSpacing*L*K/4;  % Sampling rate (must be larger than the subcarrier spacing times L)

%% OFDM Object
OFDM = OFDM(...
    L,...                                 % Number of subcarriers
    K/2*3,...                             % Number of OFDM Symbols
    SubcarrierSpacing,...                 % Subcarrier spacing (Hz)
    SamplingRate,...                      % Sampling rate (Samples/s)
    0,...                                 % Intermediate frequency first subcarrier (Hz)
    false,...                             % Transmit real valued signal
    1/(SubcarrierSpacing*K), ...          % Cyclic prefix length (s) 1/SubcarrierSpacing/(K/2-1)
    (8-1/2)*1/SubcarrierSpacing*1/2 ...   % Zero guard length (s)
    );


%% Alamouti Object
Alamouti = SpaceCoding(...
    'Alamouti2x1',...                     % Space time coding method
    1 ...                                 % Frequency spreading = 0; time spreading = 1
    );

%% Modulation Object
ModObjecct = SignalConstellation(ModulationOrder,'CPM');
% For ML Detection
ML_MapIndex1 = reshape(repmat((1:ModObjecct.ModulationOrder),ModObjecct.ModulationOrder,1),1,ModObjecct.ModulationOrder^2);
ML_MapIndex2 = reshape(repmat((1:ModObjecct.ModulationOrder).',1,ModObjecct.ModulationOrder),1,ModObjecct.ModulationOrder^2);
ML_Mapping = ModObjecct.SymbolMapping([ML_MapIndex1;ML_MapIndex2]);
%% Pilot Matrix: 0=DataSymbol, 1=PilotSymbol, -1=ZeroSymbol;
PilotMatrixBlockAntenna1 = zeros(L,K/2);
PilotMatrixBlockAntenna1(2:6:end,1:8:end)=1; 
PilotMatrixBlockAntenna1(5:6:end,5:8:end)=1;
PilotMatrixBlockAntenna1(2:6:end,2:8:end)=-1; 
PilotMatrixBlockAntenna1(5:6:end,6:8:end)=-1;
PilotMatrixBlockAntenna2 = PilotMatrixBlockAntenna1*(-1);
PilotMatrixBlockAntenna(:,:,1) = PilotMatrixBlockAntenna1;
PilotMatrixBlockAntenna(:,:,2) = PilotMatrixBlockAntenna2;

NrPilots = sum(PilotMatrixBlockAntenna1(:)==1);
NrDataSymbols = sum(PilotMatrixBlockAntenna1(:)==0);
NrTransmittedSymbols = length(PilotMatrixBlockAntenna1(:));
%% Channel Estimation Objects
ChannelEstimation_OFDM = PilotSymbolAidedChannelEstimation(...
    'Diamond',...                           % Pilot pattern
    [...                                    % Matrix that represents the pilot pattern parameters
    L,...                 % Number of subcarriers
    6; ...                                  % Pilot spacing in the frequency domain
    K/2*3,...                   % Number of OFDM Symbols
    4 ...                                   % Pilot spacing in the time domain
    ],...                                   
    'linear'...                             % Interpolation(Extrapolation) method 'linear','spline','FullAverage,'MovingBlockAverage'
    );
%-------------------------------------------------%
tic;
%% Simulate Over Different Channel Realizations
for i_rep = 1:100:NrRepetitions
    
%% Generate Data and Pilots
% Pilot Symbols: The pilot symbol power is by a factor of two higher because (pilots for the other antenna) are zero!
x_PilotAntenna1 = ModObjecct.SymbolMapping(randi(ModObjecct.ModulationOrder,NrPilots,1));
x_PilotAntenna2 = ModObjecct.SymbolMapping(randi(ModObjecct.ModulationOrder,NrPilots,1));
x_PilotAntenna1 = x_PilotAntenna1./abs(x_PilotAntenna1)*sqrt(2);
x_PilotAntenna2 = x_PilotAntenna2./abs(x_PilotAntenna2)*sqrt(2);

% Binary Data Stream
BinaryDataStream_Alamouti = randi([0 1],NrDataSymbols*log2(ModObjecct.ModulationOrder),1);
BinaryDataStream_SMAntenna1 = randi([0 1],NrDataSymbols*log2(ModObjecct.ModulationOrder),1); %Spatial Multiplexing
BinaryDataStream_SMAntenna2 = randi([0 1],NrDataSymbols*log2(ModObjecct.ModulationOrder),1); %Spatial Multiplexing
 %figure(1);
 %stem(BinaryDataStream_Alamouti);
 %title('Input Signal');
% Transmitted Alamouti Symbols
x_Alamouti = nan(L,K/2);
x_Alamouti(PilotMatrixBlockAntenna1==0) = ModObjecct.Bit2Symbol(BinaryDataStream_Alamouti);
x_Alamouti_Coded = Alamouti.Encoder(x_Alamouti);
x_Alamouti_Coded(PilotMatrixBlockAntenna==1)=[x_PilotAntenna1;x_PilotAntenna2];
x_Alamouti_Coded(PilotMatrixBlockAntenna==-1)=0;
x_Alamouti_Coded_Antenna1 = x_Alamouti_Coded(:,:,1);
x_Alamouti_Coded_Antenna2 = x_Alamouti_Coded(:,:,2);

%figure(2);
% subplot(121);surf(abs(x_Alamouti_Coded_Antenna1));
% title('Equalized Signal1');
% subplot(122);surf(abs(x_Alamouti_Coded_Antenna2));
% title('Equalized Signal2');
 
% Transmitted Spatial Multiplexed Symbols
x_SM_Antenna1 = nan(L,K/2);
x_SM_Antenna1(PilotMatrixBlockAntenna1==0) = ModObjecct.Bit2Symbol(BinaryDataStream_SMAntenna1);
x_SM_Antenna1(PilotMatrixBlockAntenna1==1) = x_PilotAntenna1;
x_SM_Antenna1(PilotMatrixBlockAntenna1==-1) = 0;
x_SM_Antenna2 = nan(L,K/2);
x_SM_Antenna2(PilotMatrixBlockAntenna2==0) = ModObjecct.Bit2Symbol(BinaryDataStream_SMAntenna2);
x_SM_Antenna2(PilotMatrixBlockAntenna2==1) = x_PilotAntenna2;
x_SM_Antenna2(PilotMatrixBlockAntenna2==-1) = 0;

% Data symbols of the the first and second block (chosen randomly to keep it simple)
x_Block1_Antenna1 = ModObjecct.SymbolMapping(randi(ModObjecct.ModulationOrder,L,K/2,1));
x_Block1_Antenna2 = ModObjecct.SymbolMapping(randi(ModObjecct.ModulationOrder,L,K/2,1));
x_Block2_Antenna1 = ModObjecct.SymbolMapping(randi(ModObjecct.ModulationOrder,L,K/2,1));
x_Block2_Antenna2 = ModObjecct.SymbolMapping(randi(ModObjecct.ModulationOrder,L,K/2,1));

%figure(3);
% subplot(221);surf(abs(x_Block1_Antenna1));
% title('Antenna 1 SM block 1');
% subplot(222);surf(abs(x_Block2_Antenna1));
% title('Antenna 2 SM block 2');
% subplot(223);surf(abs(x_Block1_Antenna2));
% title('Antenna 2 SM block 1');
% subplot(224);surf(abs(x_Block2_Antenna2));
% title('Antenna 2 SM block 2');
 
% Account for pilots in block 1 and 3
x_Block1_Antenna1(PilotMatrixBlockAntenna1==1) = x_Block1_Antenna1(PilotMatrixBlockAntenna1==1)./abs(x_Block1_Antenna1(PilotMatrixBlockAntenna1==1))*sqrt(2);
x_Block1_Antenna2(PilotMatrixBlockAntenna2==1) = x_Block1_Antenna2(PilotMatrixBlockAntenna2==1)./abs(x_Block1_Antenna2(PilotMatrixBlockAntenna2==1))*sqrt(2);
x_Block1_Antenna1(PilotMatrixBlockAntenna1==-1) = 0;
x_Block1_Antenna2(PilotMatrixBlockAntenna2==-1) = 0;
x_Block2_Antenna1(PilotMatrixBlockAntenna1==1) = x_Block2_Antenna1(PilotMatrixBlockAntenna1==1)./abs(x_Block2_Antenna1(PilotMatrixBlockAntenna1==1))*sqrt(2);
x_Block2_Antenna2(PilotMatrixBlockAntenna2==1) = x_Block2_Antenna2(PilotMatrixBlockAntenna2==1)./abs(x_Block2_Antenna2(PilotMatrixBlockAntenna2==1))*sqrt(2);
x_Block2_Antenna1(PilotMatrixBlockAntenna1==-1) = 0;
x_Block2_Antenna2(PilotMatrixBlockAntenna2==-1) = 0;

%figure(4);
% subplot(221);surf(abs(x_Block1_Antenna1));
% title('Antenna-1 Blk 1 Modulated(Pilots)');
% subplot(222);surf(abs(x_Block2_Antenna1));
% title('Antenna 1 Blk 2 Modulated(Pilots)');
% subplot(223);surf(abs(x_Block1_Antenna2));
% title('Antenna 2 Blk 1 Modulated(Pilots)');
% subplot(224);surf(abs(x_Block2_Antenna2));
% title('Antenna 2 Blk 2 Modulated(Pilots)');
%% Transmitted Signal
TransmittedSymbols_OFDM_Alamouti_Antenna1 = [x_Block1_Antenna1 x_Alamouti_Coded_Antenna1 x_Block2_Antenna1];
TransmittedSymbols_OFDM_Alamouti_Antenna2 = [x_Block1_Antenna2 x_Alamouti_Coded_Antenna2 x_Block2_Antenna2];
s_OFDM_Alamouti_Antenna1 = 1/sqrt(2)*OFDM.Modulation(TransmittedSymbols_OFDM_Alamouti_Antenna1);
s_OFDM_Alamouti_Antenna2 = 1/sqrt(2)*OFDM.Modulation(TransmittedSymbols_OFDM_Alamouti_Antenna2);

TransmittedSymbols_OFDM_SM_Antenna1 = [x_Block1_Antenna1 x_SM_Antenna1 x_Block2_Antenna1];
TransmittedSymbols_OFDM_SM_Antenna2 = [x_Block1_Antenna2 x_SM_Antenna2 x_Block2_Antenna2];
s_OFDM_SM_Antenna1 = 1/sqrt(2)*OFDM.Modulation(TransmittedSymbols_OFDM_SM_Antenna1);
s_OFDM_SM_Antenna2 = 1/sqrt(2)*OFDM.Modulation(TransmittedSymbols_OFDM_SM_Antenna2);

%figure(5);
% subplot(221);stem(s_OFDM_Alamouti_Antenna1);
% title('Antenna 1 IFFT-UB');
% subplot(222);stem(s_OFDM_Alamouti_Antenna2);
% title('Antenna 2 IFFT-UB');
% subplot(223);stem(s_OFDM_SM_Antenna1);
% title('Antenna 1 IFFT-Gauss ');
% subplot(224);stem(s_OFDM_SM_Antenna2);
% title('Antenna 2 IFFT-Gauss)');

%% Visible Light Channel 

phi = 6;%irradiant angle
space = 2.5;%space between lights

phi =(phi/180 * pi);
%-------------------------------------------------------
%---------------------ENTER PARAMETERS------------------------%
% Distance between tx and rx ( Meter )
heightLED = 1.48;

% Speed of Light
c = 300E6; 
% Time
t = 0:0.01:4;
% Radius of light cone
radius = heightLED * tan(phi);
%---------------------END OF PARAMETERS-----------------------%
%%%%%%%%%%%%%%%%%% Meshgrid X-axis and Y-axis %%%%%%%%%%%%%%%
% Incidence angles of receiver according to X-Y axis % 
[X,Y] = meshgrid(-4:0.2:4); 
xydist = sqrt((X-space).^2 + (Y-space).^2);
hdist = sqrt(xydist.^2 + heightLED.^2);
incidence = atand(xydist.* heightLED ^(-1));
%%%%%%%%%%%%%%%%%% End Mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%Plot light coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);hold on
[X1,Y1]=circle([space/2,space/2],radius,150);
plot(X1,Y1,'+')
% axis([-4 4 -4 4])
xlabel('Length of surface [m]')
ylabel('Breadth of surface [m]')
title('Top View of LEDs');hold on
[X1,Y1]=circle([space/2,-space/2],radius,150);
plot(X1,Y1,'+')
hold on
[X1,Y1]=circle([-space/2,space/2],radius,150);
plot(X1,Y1,'+')
hold on
[X1,Y1]=circle([-space/2,-space/2],radius,150);
plot(X1,Y1,'+')
hold off
%%%%%%%%%%%%%%%%%%%%% End light coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Received Power P1 2 3 4 in mW at each X-Y location%%%%%%
%%%%%%%%%%%P1
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location %
P1=P;
Z1=Z;
xydist = sqrt((X+space).^2 + (Y+space).^2);
hdist = sqrt(xydist.^2 + heightLED.^2); 
incidence = atand(xydist.* heightLED ^(-1));
%%%%%%%%%%%P2
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location % 
P2=P;
Z2=Z;
%%%%%%%%%%%P3
xydist = sqrt((X-space).^2 + (Y+space).^2);
hdist = sqrt(xydist.^2 + heightLED.^2); 
incidence = atand(xydist.* heightLED ^(-1));
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location % 
P3=P; 
Z3=Z;
%%%%%%%%%%%P4
xydist = sqrt((X+space).^2 + (Y-space).^2);
hdist = sqrt(xydist.^2 + heightLED.^2); 
incidence = atand(xydist.* heightLED ^(-1));
[P, PO, Z, impulset, impulsetd, impulsef, impulsefd]=RxSNR(incidence,hdist,t,phi); % SNR in dB at each X-Y location % 
P4=P;
Z4=Z;
%%%%%%%%%%%P TOTAL
P=P1+P2+P3+P4;
Z=Z1+Z2+Z3+Z4;
%%%%%%%%%%%%%%%%%%%%%%%END RxSNR TOTAL RECEIVED POWER%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%PLot 3D diagram for DATARATE  %%%%%%%%%%%%%%%%%%%%%%%%%%
C= 50*10^6*log2(1+Z);
figure(2);
[X1,Y1] = meshgrid(-5:0.25:5); 
surf(X1,Y1,C)
colormap(hsv);
shading interp
xlabel('Length of room [m]')
ylabel('Width of room [m]')
zlabel('Datarate in (bits/sec)')
title('DataRate Distribution')
hold off

%%%%%%%%%%%%%%%%%%%%%Plot 3D diagram for Received Power%%%%%%%%%%%%%%%%%
figure(3);
surf(X,Y,P)
colormap(hsv);
shading interp
axis([-5 5 -5 5 0 1.5e-4])
xlabel('Length of room [m]')
ylabel('Width of room [m]')
zlabel('Received Power in (W)')
title('Receivered Power Distribution')
hold off
R_space = 0.5;
xlocation =0.7;
ylocation =0.5;
% underwater channel propagation
signal=[real(s_OFDM_SM_Antenna1) real(s_OFDM_SM_Antenna2)];
% UnderwaterOC_propagation(signal)

%%%%%%%%%%%%%%%%%%%%%%%%%% GAIN VALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h11=channelgain(xlocation-R_space-space,ylocation+R_space+space,phi,space);
h12=channelgain(xlocation+R_space-space,ylocation+R_space+space,phi,space);
h13=channelgain(xlocation+R_space-space,ylocation-R_space+space,phi,space);
h14=channelgain(xlocation-R_space-space,ylocation-R_space+space,phi,space);

h21=channelgain(xlocation-R_space+space,ylocation+R_space+space,phi,space);
h22=channelgain(xlocation+R_space+space,ylocation+R_space+space,phi,space);
h23=channelgain(xlocation+R_space+space,ylocation-R_space+space,phi,space);
h24=channelgain(xlocation-R_space+space,ylocation-R_space+space,phi,space);

h31=channelgain(xlocation-R_space+space,ylocation+R_space-space,phi,space);
h32=channelgain(xlocation+R_space+space,ylocation+R_space-space,phi,space);
h33=channelgain(xlocation+R_space+space,ylocation-R_space-space,phi,space);
h34=channelgain(xlocation-R_space+space,ylocation-R_space-space,phi,space);

h41=channelgain(xlocation-R_space-space,ylocation+R_space-space,phi,space);
h42=channelgain(xlocation+R_space-space,ylocation+R_space-space,phi,space);
h43=channelgain(xlocation+R_space-space,ylocation-R_space-space,phi,space);
h44=channelgain(xlocation-R_space-space,ylocation-R_space-space,phi,space);
%%%%%%%%%%%%%%%%%%%%%%%%%% END GAIN VALUE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=[h11,h12,h13,h14;h11,h22,h23,h24;h31,h32,h33,h34;h41,h42,h43,h44];

figure(4);
subplot(211);surf(h)
colormap(jet);
shading interp
xlabel('Number of Transmit antenna')
ylabel('Number of Receive antenna')
zlabel('Channel Gain matrix')
title('VLC channel Gain matrix')
hold off
subplot(212);
bar(reshape(h,1,16));
xlabel('Time(s)')
ylabel('Normalized Intensity')
title('Impulse response of fading channel')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ptx =0.1;

Prx = Ptx * h;
% (supportet by our measurements for a small number of subcarriers and multicarrier symbols)
h11 = sqrt(1/2)*(randn+1j*randn)*mean(Prx(1,:))*0.2e6;
h12 = sqrt(1/2)*(randn+1j*randn)*mean(Prx(1,:))*0.2e6;   
h21 = sqrt(1/2)*(randn+1j*randn)*mean(Prx(1,:))*0.2e6;
h22 = sqrt(1/2)*(randn+1j*randn)*mean(Prx(1,:))*0.2e6;   
powerPbit=length(M_SNR_dB_OFDM);
%% Receiver Part
for i_S = 1:powerPbit
SNR_dB_OFDM = M_SNR_dB_OFDM(i_S); 

Pn = SamplingRate/(SubcarrierSpacing*L)*10^(-SNR_dB_OFDM/10); 
noise_Antenna1 = sqrt(Pn/2)*(randn(size(s_OFDM_Alamouti_Antenna1))+1j*randn(size(s_OFDM_Alamouti_Antenna1)));
noise_Antenna2 = sqrt(Pn/2)*(randn(size(s_OFDM_Alamouti_Antenna1))+1j*randn(size(s_OFDM_Alamouti_Antenna1)));


%% Received Signal
r_OFDM_Alamouti_Antenna1 = h11*s_OFDM_Alamouti_Antenna1+h12*s_OFDM_Alamouti_Antenna2+noise_Antenna1;

r_OFDM_SM_Antenna1 = h11*s_OFDM_SM_Antenna1+h12*s_OFDM_SM_Antenna2+noise_Antenna1;
r_OFDM_SM_Antenna2 = h21*s_OFDM_SM_Antenna1+h22*s_OFDM_SM_Antenna2+noise_Antenna2;

%% demodulation
 y_OFDM_Alamouti_3Blocks = OFDM.Demodulation(r_OFDM_Alamouti_Antenna1);
y_OFDM_Alamouti = y_OFDM_Alamouti_3Blocks(:,(1:K/2)+K/2);

 y_OFDM_SM_Antenna1_3Blocks = OFDM.Demodulation(r_OFDM_SM_Antenna1);
y_OFDM_SM_Antenna1 = y_OFDM_SM_Antenna1_3Blocks(:,(1:K/2)+K/2);
y_OFDM_SM_Antenna2_3Blocks = OFDM.Demodulation(r_OFDM_SM_Antenna2);
y_OFDM_SM_Antenna2 = y_OFDM_SM_Antenna2_3Blocks(:,(1:K/2)+K/2);
%% Noise and Power Estimation
noise_OFDM_Antenna1_3Blocks = OFDM.Demodulation(noise_Antenna1);
noise_OFDM_Antenna1 = noise_OFDM_Antenna1_3Blocks(:,(1:K/2)+K/2);
noise_OFDM_Antenna2_3Blocks = OFDM.Demodulation(noise_Antenna2);
noise_OFDM_Antenna2 = noise_OFDM_Antenna2_3Blocks(:,(1:K/2)+K/2);

Pn_OFDM_Antenna1(i_S,i_rep) = mean(abs(noise_OFDM_Antenna1(:)).^2);
Pn_OFDM_Antenna2(i_S,i_rep) = mean(abs(noise_OFDM_Antenna2(:)).^2);

PSignalPlusNoise_OFDM_Alamouti(i_S,i_rep) = mean(abs(y_OFDM_Alamouti(:)).^2);

PSignalPlusNoise_OFDM_SM_Antenna1(i_S,i_rep) = mean(abs(y_OFDM_SM_Antenna1(:)).^2);
PSignalPlusNoise_OFDM_SM_Antenna2(i_S,i_rep) = mean(abs(y_OFDM_SM_Antenna2(:)).^2);

%% Channel Estimation
EstimatedChannel_OFDM_Alamouti(:,:,1) = mean(y_OFDM_Alamouti(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1)*ones(L,K/2);
EstimatedChannel_OFDM_Alamouti(:,:,2) = mean(y_OFDM_Alamouti(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2)*ones(L,K/2);

% OFDM  :  PowerNormalization (Note that the noise power is also increasd=> The same SNR (ignoring zero guard))

H11_Est_OFDM = mean(y_OFDM_SM_Antenna1(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
H21_Est_OFDM = mean(y_OFDM_SM_Antenna2(PilotMatrixBlockAntenna1==1)./x_PilotAntenna1);
H12_Est_OFDM = mean(y_OFDM_SM_Antenna1(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);
H22_Est_OFDM = mean(y_OFDM_SM_Antenna2(PilotMatrixBlockAntenna2==1)./x_PilotAntenna2);


if PerfectChannelKnowledge
    % The power is split between TX Antennta 1 and 2, which shows in the channel   
    EstimatedChannel_OFDM_Alamouti(:,:,1) = h11/sqrt(2)*ones(L,K/2);
    EstimatedChannel_OFDM_Alamouti(:,:,2) = h12/sqrt(2)*ones(L,K/2);
    
    H11_Est_OFDM = h11/sqrt(2);
    H21_Est_OFDM = h21/sqrt(2);
    H12_Est_OFDM = h12/sqrt(2);
    H22_Est_OFDM = h22/sqrt(2);

end

H_Est_OFDM =[H11_Est_OFDM,H12_Est_OFDM;H21_Est_OFDM,H22_Est_OFDM];

%% Data Detection
%Estimated Symbols

x_est_OFDM_Alamouti = Alamouti.Decoder(y_OFDM_Alamouti,EstimatedChannel_OFDM_Alamouti*sqrt(2));


% ML Detection
y_OFDM_SM_Temp = repmat(reshape([y_OFDM_SM_Antenna1(:).';y_OFDM_SM_Antenna2(:).'],2,1,[]),1,ModObjecct.ModulationOrder^2);
[~,indexMin] = min(sum(abs(y_OFDM_SM_Temp-reshape(repmat(H_Est_OFDM*ML_Mapping,1,L*K/2),size(y_OFDM_SM_Temp))).^2),[],2);
x_est_OFDM_SM_ML = reshape(ML_Mapping(:,indexMin(:)).',L,K/2,2);


% Symbols To Bit
DetectedBitStream_OFDM_Alamouti = ModObjecct.Symbol2Bit(x_est_OFDM_Alamouti(PilotMatrixBlockAntenna1==0));
DetectedBitStream_OFDM_SM_ML = ModObjecct.Symbol2Bit(x_est_OFDM_SM_ML(PilotMatrixBlockAntenna==0));

% Bit Error Ratio
BER_OFDM_Alamouti(i_S,i_rep) = mean(BinaryDataStream_Alamouti~=DetectedBitStream_OFDM_Alamouti);
BER_OFDM_SM_ML(i_S,i_rep) = mean([BinaryDataStream_SMAntenna1;BinaryDataStream_SMAntenna2]~=DetectedBitStream_OFDM_SM_ML);


 
end

TimePassed = toc;
if mod(i_rep,100)==0
disp(['Realization ' int2str(i_rep) ' of ' int2str(NrRepetitions) '. Time left: ' int2str(TimePassed/i_rep*(NrRepetitions-i_rep)/60) 'minutes']);
end
end
UnderwaterOC_propagation(signal);


%figure();
%ChannelEstimation_OFDM.PlotPilotPattern;
%title('Symbols of MIMO');

[Power_OFDM,t_OFDM] = OFDM.PlotTransmitPower;
%figure();
%plot(t_OFDM,Power_OFDM,'b-s','linewidth',2);
%ylabel('Transmit Power');
%xlabel('Time(s)');
%title('SM-MIMO signal transmission power');

%% Calculate Power Spectral Density
[PSD_OFDM,t_OFDM] = OFDM.PlotPowerSpectralDensity;
%figure();
%plot(t_OFDM,fftshift(10*log10((PSD_OFDM))),'red ');
%ylabel('Power Spectral Density (dB)');
%xlabel('Frequency (Hz)');
%title('Spectrum of UWOC signal');
Disin=M_SNR_dB_OFDM;
load Paper_result
figure();
semilogy(Disin,mean(BER_OFDM_SM_ML,2),'-s black','linewidth',2);
hold on;
semilogy(Disin,mean(BER_OFDM_Alamouti,2),'-o green','linewidth',2);
hold on;
semilogy(Disin,D1,'-s red','linewidth',2);
hold on;
semilogy(Disin,D2,'-o blue','linewidth',2);
grid on;
legend('MIMO 2x2 with EGC-CPM','MIMO 2x2 with OC-CPM','MIMO 2x2 with EGC-OOK',...
    'MIMO 2x2 with OC-OOK','location','best');
xlabel('Average transmitted power per bit[dBm]');
ylabel('Bit Error Ratio');
title('BER of 25 m coastal water link');

%{
figure();
semilogy(Disin,BER_OFDM_SM_ML(:,1:100:800),'-s','linewidth',2);
grid on;
legend('MIMO 2x2-UB with approx','MISO 3x1-UB with approx',...
    'MIMO 2x2-UB without approx','MISO 3x1-UB without approx',...
    'MIMO 2x2 -Gaussian approx','MISO 3x1-Gaussian approx ',...
    'MIMO 2x2 -Saddle point approx ','MISO 3x1-Saddle point approx',...
    'location','best');
xlabel('Average transmitted power per bit[dBm]');
ylabel('Bit Error Ratio');
title('BER of different approximations with \sigma_x=0.3');
%}

figure;semilogy(M_SNR_dB_OFDM,BER_OFDM_SM_ML(:,1:100:500),'-d','linewidth',2);
grid on;
legend('MIMO 2x2-Sim','MISO 2x1-Sim','SISO-Analytic','MIMO 2x2-Analytic','MISO 2x1-Analytic','location','best');
xlabel('Average transmitted power per bit[dBm]');
ylabel('Bit Error Ratio');
title('BER of 25 m coastal water link with \sigma_x=0.1');
figure;semilogy(M_SNR_dB_OFDM,BER_OFDM_SM_ML(:,501:100:end),'-p','linewidth',2);
grid on;
legend('MIMO 2x2-Sim','MISO 2x1-Sim','SISO-Analytic','MIMO 2x2-Analytic','MISO 2x1-Analytic','location','best');
xlabel('Average transmitted power per bit[dBm]');
ylabel('Bit Error Ratio');
title('BER of 25 m coastal water link with \sigma_x=0.4');


figure;semilogy(M_SNR_dB_OFDM,BER_OFDM_Alamouti(:,1:100:500),'-d','linewidth',2);
grid on;
legend('MIMO 2x2','MISO 2x1','SISO','MISO 3x1','SIMO 1x2','location','best');
xlabel('Average transmitted power per bit[dBm]');
ylabel('Bit Error Ratio');
title('BER of 25 m coastal water link with R_b=0.5Gbps');
figure;semilogy(M_SNR_dB_OFDM,BER_OFDM_Alamouti(:,501:100:end),'-p','linewidth',2);
grid on;
legend('MISO 2x1','MISO 3x1','MIMO 2x2','SISO','SIMO 1x2','location','best');
xlabel('Average transmitted power per bit[dBm]');
ylabel('Bit Error Ratio');
title('BER of 25 m coastal water link with R_b=50Gbps');

%figure;
%[C1 h1]=contourf(incidence(31:36,5:12)*1e-1);
%clabel (C1,h1);colormap(bone);
%colorbar
%xlabel('Horizontal Displacemnt');
%ylabel('Vertical Displacement');
%title('BER at various points on a visible light plane');