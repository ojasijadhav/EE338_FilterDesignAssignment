clear all;
m = 51;
transbw = 3*10^3;
qm = floor(0.1*m);
rm = m - 10*qm;
BL = 5+3*qm + 11*rm;
BH = BL + 25;

% Band Edge specifications
freqsamp = 400e3;
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


%Parameters for Bandpstop Transformation
W0 = sqrt(w_p1*w_p2);
B = w_p2-w_p1;
wl_s1 = ((w_s1)^2 - (W0)^2)/(B*w_s1);
wl_s2 = ((w_s2)^2 - (W0)^2)/(B*w_s2);
wl_s = min(abs(wl_s1),abs(wl_s2));


%LPF parameters
D1 = 1/((1-delta)^2)-1;
D2 = 1/(delta)^2 - 1;
%Chebyschev approximation
N = acosh(sqrt(D2/D1))/acosh(wl_s); 
N = ceil(N);
w_c = ((1/D1^(1/(2*N))) + (wl_s/(D2^(1/(2*N)))))/2;

syms x;
lowpf2 = (1+(D1*(chebyshevT(N,x/1i))^2));
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
den_bandsf = 1;
den_filter = 1;
norm = 1;
for i=1:N
    den_lowpf = den_lowpf*(s-rootslowpf(i));
    den_bandsf = den_bandsf*(((B*s)/(s^2+W0^2))-rootslowpf(i));
    den_filter = den_filter*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-rootslowpf(i));
    norm = norm*rootslowpf(i);
end
lowpf = (norm)/(den_lowpf);
bandsf = (norm)/(den_bandsf);
des_filter = (norm)/(den_filter);
%[num,den]  = tfdata(lowpf,'v');
%freqs(num,den);
%[num,den]  = tfdata(bandpf,'v');
%freqs(num,den);
%
%freqs(num,den);
%discrete
[num,den]  = tfdata(des_filter,'v');
fvtool(num,den);