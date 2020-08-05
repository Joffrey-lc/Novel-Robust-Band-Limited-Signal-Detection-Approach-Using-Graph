clear;
T0=1;%符号周期
MPSK=2;%BPSK
r=10;%quantization level
K=300;%simple size
M=300;%FFT size
Gamma=50;%over-simpling factor
fs=Gamma/T0;%采样率
t=0:1/fs:T0-1/fs;%时间向量
fc=2/T0; %载波频率  也是一个码元多少个正弦周期
c=sqrt(2)*exp(1i*2*pi*fc*t);%载波信号
nsymbol=K/Gamma;%每种信噪比下的发送符号数
theta=1;%menta carlo number
r_set=1/r:1/r:1;%量化区间
snr_dB=[-3 -7 -11 -1000];
msg=randi([0 MPSK-1],1,nsymbol); %生成基带数据       
msgmod=pskmod(msg,MPSK).'; %基带B-PSK调制
tx=real(msgmod*c);%载波调制
tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'的每一列是一个码元代表的采样点,现展开为一行       
for indx=1:length(snr_dB)
    lamda0=zeros(1,theta);
    for indx2=1:1:theta
        rx=noisegen_err(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rxy=abs(fft(rx,300));%fft
        figure(1)%显示fft图像
        subplot(4,1,indx)
        plot(rxy);
        Ux=zeros(1,length(rxy));
        r_level=zeros(1,length(rxy));
    %%%Normalized
        theta_min=min(rxy);
        theta_max=max(rxy);
        for m=1:1:length(rxy)
            Ux(m)=(rxy(m)- sita_min)/(sita_max-sita_min);  
        end
    %%%quantization
        for mm=1:1:length(Ux)
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%找到量化等级
        end
    figure(2);%显示量化等级图像 同论文中的Fig.1
    subplot(4,1,indx)
    plot(r_level);
    title(['SNR=',num2str(snr_dB(indx)),'dB']);
    xlabel('M','position',[320 -20]);
    ylabel('Qx(m)');
        Lx=get_LaplacianMatrix(r,r_level);%得到laplacian 矩阵
        [~,lamda]=eig(Lx);%计算特征值
        [not_sort,~]=max(lamda);%提取特征值
        lamda_sort=sort(not_sort);%特征值排序
        lamda0(indx2)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave=sum(lamda0)/length(lamda0);
end
