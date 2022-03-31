clear all;
m = 51;
qm = floor(0.1*m);
rm = m - 10*qm;
BL = 5+3*qm + 11*rm;
BH = BL + 25;
transbw = 3*10^3;

% Band Edge specifications
freq_p1 = BL*10^3-transbw;
freq_s1 = BL*10^3;
freq_s2 = BH*10^3;
freq_p2 = BH*10^3+transbw;
freqsamp = 400e3;
delta = 0.15;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(freq_s1+freq_p1)/2 (freq_s2+freq_p2)/2]*2/freqsamp;        %average value of the two paramters
N_min = 1 + ceil((A-7.95) / (2.285*(6/400)*pi));       %empirical formula for N_min
disp(N_min);
%Window length for Kaiser Window
n=N_min+15;
disp(N_min);
disp(n);

%Ideal bandstop impulse response of length "n"

bsideal =  ideal_lp(pi,n) -ideal_lp(((freq_p2+freq_s2)/freqsamp)*pi,n) + ideal_lp(((freq_s1+freq_p2)/freqsamp)*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_window = (kaiser(n,beta))';

FIR_BandStop = bsideal .* kaiser_window;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, freqsamp);
plot(f,abs(H))
grid