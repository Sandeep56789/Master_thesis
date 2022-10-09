function y = TVReceiverReconstruct(x,h,fs_x,fs_h)

%
% implement a time-varying Reconstruction
%
% inputs :  - x : input signal, length Nx vector
%           - h : the initial TV impulse response, matrix Nh*M (time*delay)
%           - fs_x : input signal's sampling frequencies along time
%           - fs_h : IR's sampling frequencies along time
%
% output :  - y : the output, length Nx+M-1 vector
%
% uses cubic splines to interpolate the time-varying impulse response h(t,tau)
% along the time t.
% only one hydrophone is considered
%


% parameters %
%============%
x = x(:);
Nx = length(x);                     % length of the input signal
[Nh,M] = size(h);                   % dimension of the initial IR %
y = zeros(Nx+M-1,1);                % output signal
x = [zeros(M-1,1);x;zeros(M-1,1)];  % zero-pad the input signal %

interpsolution = 2;     % interpolation choice :    1 => amplitude and phase
                        %                           2 => real and imaginary parts

if(interpsolution == 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation in amplitude and phase %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spline polynomial coefficients evaluation %
%===========================================%
ppamp = spline([0:Nh-1]/fs_h,abs(shiftdim(h,1)));
ppphi = spline([0:Nh-1]/fs_h,unwrap(angle(shiftdim(h,1))));

% perform time varying filtering %
%================================%
for t = 1:Nx+M-1
    hcurrent = flipud(ppval(ppamp,t/fs_x).*exp(sqrt(-1)*ppval(ppphi,t/fs_x)));
    y(t) = sum(x(t:t+M-1).*hcurrent,1)/fs_x; 
end

elseif(interpsolution == 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation in real and imaginary parts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spline polynomial coefficients evaluation %
%===========================================%
pp = spline([0:Nh-1]/fs_h,shiftdim(h,1));

% perform time varying filtering %
%================================%
for t = 1:Nx+M-1
    hcurrent = flipud(ppval(pp,t/fs_x));
    y(t) = sum(x(t:t+M-1).*hcurrent,1)/fs_x; 
end
y=y(1:Nx)*5e2;
end
