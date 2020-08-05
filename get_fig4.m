clear;
T0=1;%符号周期
MPSK=2;%BPSK
K=160;%simple size
M=300;%FFT size
Gamma=4;%over-simpling factor
fs=Gamma/T0;%采样率
t=0:1/fs:T0-1/fs;%时间向量
fc=2/T0; %载波频率  也是一个码元多少个正弦周期
c=sqrt(2)*exp(1i*2*pi*fc*t);%载波信号
nsymbol=K/fs;%每种信噪比下的发送符号数
theta=25;%menta carlo number
snr_dB=-20:2:10;
msg=randi([0 MPSK-1],1,nsymbol); %生成基带数据       
msgmod=pskmod(msg,MPSK).'; %基带B-PSK调制
tx=real(msgmod*c);%载波调制
tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'的每一列是一个码元代表的采样点,现展开为一行 

r=10;%quantization level
r_set=1/r:1/r:1;%量化区间
for indx=1:length(snr_dB)
    lamda1=zeros(1,theta);
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rxy=abs(fft(rx,M));%fft
        Ux=zeros(1,length(rxy));
        r_level=zeros(1,length(rxy));
    %%%Normalized
        sita_min=min(rxy);
        sita_max=max(rxy);
        for m=1:1:length(rxy)
            Ux(m)=(rxy(m)- sita_min)/(sita_max-sita_min);  
        end
    %%%quantization
        for mm=1:1:length(Ux)
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%找到量化等级
        end
        Lx=get_LaplacianMatrix(r,r_level);%得到laplacian 矩阵
        [~,lamda]=eig(Lx);%计算特征值
        [not_sort,~]=max(lamda);%提取特征值
        lamda_sort=sort(not_sort);%特征值排序
        lamda1(indx2)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave_BPSK(indx)=sum(lamda1)/length(lamda1);
    
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=rx-tx1;%AWGN
        rxy=abs(fft(rx,M));%fft
        Ux=zeros(1,length(rxy));
        r_level=zeros(1,length(rxy));
    %%%Normalized
        sita_min=min(rxy);
        sita_max=max(rxy);
        for m=1:1:length(rxy)
            Ux(m)=(rxy(m)- sita_min)/(sita_max-sita_min);  
        end
    %%%quantization
        for mm=1:1:length(Ux)
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%找到量化等级
        end
        Lx=get_LaplacianMatrix(r,r_level);%得到laplacian 矩阵
        [~,lamda]=eig(Lx);%计算特征值
        [not_sort,~]=max(lamda);%提取特征值
        lamda_sort=sort(not_sort);%特征值排序
        lamda1(indx2)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave_AWGN(indx)=sum(lamda1)/length(lamda1);
end
figure(1)
h11=plot(snr_dB,lamda_ave_BPSK,'ro-');
hold on;
h12=plot(snr_dB,lamda_ave_AWGN,'ro--');

r=15;%quantization level
r_set=1/r:1/r:1;%量化区间
for indx=1:length(snr_dB)
    lamda1=zeros(1,theta);
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rxy=abs(fft(rx,M));%fft
        Ux=zeros(1,length(rxy));
        r_level=zeros(1,length(rxy));
    %%%Normalized
        sita_min=min(rxy);
        sita_max=max(rxy);
        for m=1:1:length(rxy)
            Ux(m)=(rxy(m)- sita_min)/(sita_max-sita_min);  
        end
    %%%quantization
        for mm=1:1:length(Ux)
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%找到量化等级
        end
        Lx=get_LaplacianMatrix(r,r_level);%得到laplacian 矩阵
        [~,lamda]=eig(Lx);%计算特征值
        [not_sort,~]=max(lamda);%提取特征值
        lamda_sort=sort(not_sort);%特征值排序
        lamda1(indx2)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave_BPSK(indx)=sum(lamda1)/length(lamda1);
    
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=rx-tx1;%AWGN
        rxy=abs(fft(rx,M));%fft
        Ux=zeros(1,length(rxy));
        r_level=zeros(1,length(rxy));
    %%%Normalized
        sita_min=min(rxy);
        sita_max=max(rxy);
        for m=1:1:length(rxy)
            Ux(m)=(rxy(m)- sita_min)/(sita_max-sita_min);  
        end
    %%%quantization
        for mm=1:1:length(Ux)
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%找到量化等级
        end
        Lx=get_LaplacianMatrix(r,r_level);%得到laplacian 矩阵
        [~,lamda]=eig(Lx);%计算特征值
        [not_sort,~]=max(lamda);%提取特征值
        lamda_sort=sort(not_sort);%特征值排序
        lamda1(indx2)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave_AWGN(indx)=sum(lamda1)/length(lamda1);
end
h21=plot(snr_dB,lamda_ave_BPSK,'g+-');
hold on;
h22=plot(snr_dB,lamda_ave_AWGN,'g+--');

r=20;%quantization level
r_set=1/r:1/r:1;%量化区间
for indx=1:length(snr_dB)
    lamda1=zeros(1,theta);
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rxy=abs(fft(rx,M));%fft
        Ux=zeros(1,length(rxy));
        r_level=zeros(1,length(rxy));
    %%%Normalized
        sita_min=min(rxy);
        sita_max=max(rxy);
        for m=1:1:length(rxy)
            Ux(m)=(rxy(m)- sita_min)/(sita_max-sita_min);  
        end
    %%%quantization
        for mm=1:1:length(Ux)
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%找到量化等级
        end
        Lx=get_LaplacianMatrix(r,r_level);%得到laplacian 矩阵
        [~,lamda]=eig(Lx);%计算特征值
        [not_sort,~]=max(lamda);%提取特征值
        lamda_sort=sort(not_sort);%特征值排序
        lamda1(indx2)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave_BPSK(indx)=sum(lamda1)/length(lamda1);
    
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%加入高斯白噪声
        rx=rx-tx1;%AWGN
        rxy=abs(fft(rx,M));%fft
        Ux=zeros(1,length(rxy));
        r_level=zeros(1,length(rxy));
    %%%Normalized
        sita_min=min(rxy);
        sita_max=max(rxy);
        for m=1:1:length(rxy)
            Ux(m)=(rxy(m)- sita_min)/(sita_max-sita_min);  
        end
    %%%quantization
        for mm=1:1:length(Ux)
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%找到量化等级
        end
        Lx=get_LaplacianMatrix(r,r_level);%得到laplacian 矩阵
        [~,lamda]=eig(Lx);%计算特征值
        [not_sort,~]=max(lamda);%提取特征值
        lamda_sort=sort(not_sort);%特征值排序
        lamda1(indx2)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave_AWGN(indx)=sum(lamda1)/length(lamda1);
end
h31=plot(snr_dB,lamda_ave_BPSK,'bx-');
hold on;
h32=plot(snr_dB,lamda_ave_AWGN,'bx--');
hold off;
legend([h11(1),h12(1) h21(1) h22(1) h31(1) h32(1)],'BPSK(r=10)','AWGN(r=10)','BPSK(r=15)','AWGN(r=15)','BPSK(r=20)','AWGN(r=20)','location', 'southwest');
axis([-20 10 -5 20]);
ylabel('lamda1_bar');
xlabel('Signal-to-Noise Ratio (dB)');