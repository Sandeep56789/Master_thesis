%
% Generate the impulse response as a time serie (a signal) from the individual
% amplitudes, phases and arrival times calculated by Bellhop.
%
% Proceed in the frequency domain
%

% produce a real impulse response, with a symetric spectrum %

% produce the baseband complex impulse responses (without energy correction) %
%============================================================================%
vf = fc + [-ntau/2:ntau/2-1]/ntau*fs_x;
vf = vf(:);
for k = 1:nhyd
    %indexes = ((travel_time(k,:) > 0.0)&(~isnan(amplitude(k,:))));
    indexes = (travel_time(k,:) > 0.0);
    taum = travel_time(k,indexes);
    am = ones(ntau,1)*amplitude(k,indexes);
    phim = ones(ntau,1)*phi(k,indexes);
    Hsingle = sum(am.*exp(sqrt(-1)*phim).*exp(-sqrt(-1)*2*pi*vf*(taum-tstart)),2);
    Hsingle = [Hsingle(ntau/2+1:end);Hsingle(1:ntau/2)];
    hit(tt,:,k) = ifft(Hsingle);
    clear indexes taum am phim Hsingle
end
