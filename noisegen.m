function Y= noisegen(SIGNAL,EbN0,T0,fs)
% 添加高斯白噪声，SNR单位dB
fb=1/T0;
NOISE=randn(size(SIGNAL));
Eb =sum(SIGNAL.*SIGNAL)/(length(SIGNAL)*fb);
ebN0=10^((EbN0)/10);
N0=Eb/ebN0;
noise_variance =fs/2*N0;
Y=SIGNAL+sqrt(noise_variance)*NOISE;