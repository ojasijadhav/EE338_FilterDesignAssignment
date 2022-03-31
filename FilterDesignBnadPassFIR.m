clear all;
m = 51;
transbw = 3*10^3;
qm = floor(0.1*m-0.0001);
rm = m - 10*qm;
BL = 10+5*qm + 13*rm;
BH = BL + 45;

% Band Edge specifications
delta = 0.15;
freqsamp = 540e3;
freq_s2 = BH*10^3+transbw;
freq_s1 = BL*10^3-transbw;
freq_p1 = BL*10^3;
freq_p2 = BH*10^3;
W_c1 = freq_p1*2*pi/freqsamp;
W_c2  = freq_p2*2*pi/freqsamp;

%Kaiser paramters calculations
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

% formula for N_min
Nmin = 1 + ceil((A-7.95) / (2.285*(6/540)*pi));          

%Window length for Kaiser Window
n=Nmin ;
disp(n)

%Ideal bandpass impulse response of length "n"

bpideal = ideal_lp(((freq_p2+freq_s2)/freqsamp)*pi,n) - ideal_lp(((freq_s1+freq_p2)/freqsamp)*pi,n);


%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';
FIR_BandPass = bpideal .* kaiser_win;
%frequency response
fvtool(FIR_BandPass);        

%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, freqsamp);
plot(f,abs(H))
grid
