clear all;
m = 51;
transbw = 3*10^3;
qm = floor(0.1*m-0.0001);
rm = m - 10*qm;
BL = 10+5*qm + 13*rm;
BH = BL + 45;

% Band Edge specifications
freqsamp = 540e3;
delta = 0.15;
freq_s1 = BL*10^3-transbw;
freq_s2 = BH*10^3+transbw;
freq_p1 = BL*10^3;
freq_p2 = BH*10^3;


%Bilinear Transformation
w_s1 = tan(freq_s1/freqsamp*pi); 
w_s2 = tan(freq_s2/freqsamp*pi);
w_p1 = tan(freq_p1/freqsamp*pi);
w_p2 = tan(freq_p2/freqsamp*pi);


%Parameters for Bandpass Transformation
W0 = sqrt(w_p1*w_p2);
B = w_p2-w_p1;
wl_s1 = ((w_s1)^2 - (W0)^2)/(B*w_s1);
wl_s2 = ((w_s2)^2 - (W0)^2)/(B*w_s2);
wl_s = min(abs(wl_s1),abs(wl_s2));


%Butterworth LPF parameters
D1 = 1/((1-delta)^2)-1;
D2 = 1/(delta)^2 - 1;

%Butterworth approximation
N = log(sqrt(D2/D1))/log(wl_s); 
N = ceil(N);         %order
w_c = ((1/D1^(1/(2*N))) + (wl_s/(D2^(1/(2*N)))))/2;

syms x;
lowpf2 = (1+(x/(1i*w_c))^(2*N));
rootslowpf2 = double(solve(lowpf2));

reallowpf2 = real(rootslowpf2);
imglowpf2  = imag(rootslowpf2);
rootslowpf = zeros([N 1]);
j=1;
for i=1:2*N
    if reallowpf2(i)<0 
        rootslowpf(j) = rootslowpf2(i);
        j=j+1;
    end
    
end

%disp(roots_lpf)

s = tf('s');
z = tf('z');
den_lowpf = 1;
den_bandpf = 1;
den_filter = 1;
for i=1:N
    den_lowpf = den_lowpf*(s-rootslowpf(i));
    den_bandpf = den_bandpf*(((s^2+W0^2)/(B*s))-rootslowpf(i));
    den_filter = den_filter*(((((z-1)/(z+1))^2+W0^2)/(B*((z-1)/(z+1))))-rootslowpf(i));
end
lowpf = (w_c^N)/den_lowpf;
bandpf = (w_c^N)/den_bandpf;
des_filter = (w_c^N)/den_filter;
%[num,den]  = tfdata(lowpf,'v');
%freqs(num,den);
%[num,den]  = tfdata(bandpf,'v');
%freqs(num,den);
%
%freqs(num,den);
%discrete
[num,den]  = tfdata(des_filter,'v');
fvtool(num,den);