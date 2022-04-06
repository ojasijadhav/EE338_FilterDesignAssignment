clear all;
m = 51;
qm = floor(0.1*m-0.0001);
rm = m - 10*qm;
BL = 10+5*qm + 13*rm;
BH = BL + 45;
transbw = 3*10^3;

% Band Edge specifications
freq_s1 = BL*10^3-transbw;
freq_p1 = BL*10^3;
freq_p2 = BH*10^3;
freq_s2 = BH*10^3+transbw;
freqsamp = 540e3;
delta = 0.15;
fn_s1 = freq_s1/freqsamp*2;          
fn_p1 = freq_p1/freqsamp*2;
fn_p2 = freq_p2/freqsamp*2;
fn_s2 = freq_s2/freqsamp*2;
twn = transbw/freqsamp*2;

%Transformed specs using Bilinear Transformation
w_s1 = tan(freq_s1/freqsamp*pi);          
w_p1 = tan(freq_p1/freqsamp*pi);
w_p2 = tan(freq_p2/freqsamp*pi);
w_s2 = tan(freq_s2/freqsamp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(w_p1*w_p2);
B = w_p2-w_p1;
wl_s1 = ((w_s1)^2 - (W0)^2)/(B*w_s1);
wl_s2 = ((w_s2)^2 - (W0)^2)/(B*w_s2);
wl_s = min(abs(wl_s1),abs(wl_s2));


G_p = 0.85; G_s =0.15 ; % filter specifications
W_p = 1; W_s = wl_s;
f_p = W_p/(2*pi);
f_s = W_s/(2*pi);
ep = sqrt(1/G_p^2 - 1); es = sqrt(1/G_s^2 - 1); % ripples εp =

k = W_p/W_s; % k = 
k1 = ep/es; % k1 = 

[K,Kp] = ellipk(k); % elliptic integrals K = , K' = 
[K1,K1p] = ellipk(k1); % elliptic integrals K1 = , K'1 = 

Nexact = (K1p/K1)/(Kp/K); N = ceil(Nexact); % Nexact = , N = 
disp(N)
k = ellipdeg(N,k1); % recalculated k = 

fs_new = f_p/k; % new stopband fs = 

L = floor(N/2); r = mod(N,2); q = (1:L); % L = , r = , i = [1; ]
u = (2*q-1)/N; zeta_i = cde(u,k); % ui = [; ], ζi = [;]
disp(r)
za = W_p * 1j./(k*zeta_i); % filter zeros
disp('za')
disp(za)
v0 = -1j*asne(1j/ep, k1)/N; % v0 = 

pa = W_p * 1j*cde(u-1j*v0, k); % filter poles
disp('pa')
disp(pa)
pa0 = W_p * 1j*sne(1j*v0, k);

s= tf('s');
z=tf('z');
lowpf = ((s-za(1))*(s-za(2))*(s-conj(za(1)))*(s-conj(za(2))))/((s-pa(1))*(s-pa(2))*(s-conj(pa(1)))*(s-conj(pa(2))));
norm = (za(1)*za(2)*conj(za(1))*conj(za(2)))/((conj(pa(2))*(conj(pa(1)))*(pa(2))*(pa(1))));
dcgain  = G_p^(1-r);
lowpf = dcgain*lowpf/norm;
[numl,denl]  = tfdata(lowpf,'v');
% freqs(num,den);
[h,w] =freqs(numl,denl);
plot(w,abs(h));
figure;
plot(w,angle(h));

bandpf = ((((s^2+W0^2)/(B*s))-za(1))*(((s^2+W0^2)/(B*s))-za(2))*(((s^2+W0^2)/(B*s))-conj(za(1)))*(((s^2+W0^2)/(B*s))-conj(za(2))))/((((s^2+W0^2)/(B*s))-pa(1))*(((s^2+W0^2)/(B*s))-pa(2))*(((s^2+W0^2)/(B*s))-conj(pa(1)))*(((s^2+W0^2)/(B*s))-conj(pa(2))));
bandpf = bandpf*dcgain/(norm);

[numb,denb] = tfdata(bandpf,'v');
[hb,wb]= freqs(numb,denb);
figure;
plot(wb,abs(hb));
figure;
plot(wb,angle(hb));

bandpf_dis = (((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-za(1))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-za(2))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(za(1)))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(za(2))))/(((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-pa(1))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-pa(2))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(pa(1)))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(pa(2))));
bandpf_dis = bandpf_dis*dcgain/(norm);

[num2,den2] = tfdata(bandpf_dis,'v');
[h2,w2]= freqz(num2,den2);
figure;
plot(w2,abs(h2));
figure;
plot(w2,angle(h2));
fvtool(num2,den2);
