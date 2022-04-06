clear all;
m = 51;
qm = floor(0.1*m);
rm = m - 10*qm;
BL = 5+3*qm + 11*rm;
BH = BL + 25;
transbw = 3*10^3;

%Freqeuncy specifications
freq_p1 = BL*10^3-transbw;
freq_s1 = BL*10^3;
freq_s2 = BH*10^3;
freq_p2 = BH*10^3+transbw;
freqsamp = 400e3;
delta = 0.15;

%Bilinear Transformation
w_s1 = tan(freq_s1/freqsamp*pi);          
w_p1 = tan(freq_p1/freqsamp*pi);
w_p2 = tan(freq_p2/freqsamp*pi);
w_s2 = tan(freq_s2/freqsamp*pi);

%Bandpass Transformation
W0 = sqrt(w_p1*w_p2);
B = w_p2-w_p1;
wl_s1 = (B*w_s1)/((W0)^2 - (w_s1)^2);
wl_s2 = (B*w_s2)/((W0)^2 - (w_s2)^2);
wl_s = min(abs(wl_s1),abs(wl_s2));

% filter specifications
G_p = 0.85; G_s =0.15 ;
W_p = 1; W_s = wl_s;
freqp = W_p/(2*pi);
freqs = W_s/(2*pi);
ep = sqrt(1/G_p^2 - 1); es = sqrt(1/G_s^2 - 1); %ripple parameters

k = W_p/W_s; % k = 
k1 = ep/es; % k1 = 

[K,Kp] = ellipk(k); 
[K1,K1p] = ellipk(k1); 

Nexact = (K1p/K1)/(Kp/K); N = ceil(Nexact); 
disp(N)
k = ellipdeg(N,k1);  

fsnew = freqp/k; 

L = floor(N/2); r = mod(N,2); q = (1:L); 
u = (2*q-1)/N; zeta_i = cde(u,k); 

% filter zeros
za = W_p * 1j./(k*zeta_i); 
disp('za')
disp(za)
v0 = -1j*asne(1j/ep, k1)/N; 

% filter poles
pa = W_p * 1j*cde(u-1j*v0, k);
pa0 = W_p * 1j*sne(1j*v0, k);
disp('pa')
disp(pa)


s= tf('s');
z=tf('z');
lowpf = ((s-za(1))*(s-conj(za(1))))/(((s-pa0)^r)*(s-pa(1))*(s-conj(pa(1))));
norm = (za(1)*conj(za(1)))/((conj(pa(1)))*(pa(1))*(pa0^r));
dcgain  = G_p^(1-r);
lowpf = dcgain*lowpf/norm;
[numl,denl]  = tfdata(lowpf,'v');
%freqs(numl,denl);
%[h,w] = freqs(numl,denl);
%plot(w,abs(h));
%figure;
%plot(w,angle(h));

bandsf =  ((((B*s)/(s^2+W0^2))-za(1))*(((B*s)/(s^2+W0^2))-conj(za(1))))/((((B*s)/(s^2+W0^2))-pa0)*(((B*s)/(s^2+W0^2))-pa(1))*(((B*s)/(s^2+W0^2))-conj(pa(1))));
bandsf = bandsf*dcgain/(norm);
 
[numb,denb] = tfdata(bandsf,'v');
%[hb,wb]= freqs(numb,denb);
%figure;
%plot(wb,abs(hb));
%figure;
%plot(wb,angle(hb));
 
 
 
 bandsf_dis = ((((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-za(1))*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-conj(za(1))))/((((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-pa0)*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-pa(1))*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-conj(pa(1))));
 bandsf_dis = bandsf_dis*dcgain/(norm);
 [num2,den2] = tfdata(bandsf_dis,'v');
fvtool(num2,den2);
disp(dcgain/(norm))
disp('num2')
disp(num2)
disp('den2')
disp(den2)
pzmap(lowpf)