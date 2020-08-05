clear;
T0=1;%符号周期
MPSK=2;%BPSK
%r=[10 15 20];%quantization level
K=100:200:20000;%simple size
M=300;%FFT size
Gamma=50;%over-simpling factor
fs=Gamma/T0;%采样率
t=0:1/fs:T0-1/fs;%时间向量
fc=2/T0; %载波频率  也是一个码元多少个正弦周期
c=sqrt(2)*exp(1i*2*pi*fc*t);%载波信号
theta=1000;%menta carlo number
snr_dB=-1000;
r=10;
r_set=1/r:1/r:1;%量化区间
for kk=1:1:length(K)
    lamda1=zeros(1,theta);
    for jj=1:1:theta
        nsymbol=round(K(kk)/Gamma);%每种信噪比下的发送符号数
        msg=randi([0 MPSK-1],1,nsymbol); %生成基带数据       
        msgmod=pskmod(msg,MPSK).'; %基带B-PSK调制
        tx=real(msgmod*c);%载波调制
        tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'的每一列是一个码元代表的采样点,现展开为一行      
        rx=noisegen(tx1,snr_dB,T0,fs);%加入高斯白噪声
        rx=rx-tx1;
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
        lamda1(jj)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave1=sum(lamda1)/length(lamda1);
    h1=plot(K(kk),lamda_ave1,'r.');
    hold on;
end

r=15;
r_set=1/r:1/r:1;%量化区间
for kk=1:1:length(K)
    lamda1=zeros(1,theta);
    for jj=1:1:theta
        nsymbol=round(K(kk)/Gamma);%每种信噪比下的发送符号数
        msg=randi([0 MPSK-1],1,nsymbol); %生成基带数据       
        msgmod=pskmod(msg,MPSK).'; %基带B-PSK调制
        tx=real(msgmod*c);%载波调制
        tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'的每一列是一个码元代表的采样点,现展开为一行      
        rx=noisegen(tx1,snr_dB,T0,fs);%加入高斯白噪声
        rx=rx-tx1;
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
        lamda1(jj)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave2=sum(lamda1)/length(lamda1);
    h2=plot(K(kk),lamda_ave2,'b.');
    hold on;
end

r=20;
r_set=1/r:1/r:1;%量化区间
for kk=1:1:length(K)
    lamda1=zeros(1,theta);
    for jj=1:1:theta
        nsymbol=round(K(kk)/Gamma);%每种信噪比下的发送符号数
        msg=randi([0 MPSK-1],1,nsymbol); %生成基带数据       
        msgmod=pskmod(msg,MPSK).'; %基带B-PSK调制
        tx=real(msgmod*c);%载波调制
        tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'的每一列是一个码元代表的采样点,现展开为一行      
        rx=noisegen(tx1,snr_dB,T0,fs);%加入高斯白噪声
        rx=rx-tx1;
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
        lamda1(jj)=lamda_sort(end-1);%找到第二大特征值
    end
    lamda_ave3=sum(lamda1)/length(lamda1);
    h3=plot(K(kk),lamda_ave3,'g.');
    hold on;
end
hold off;
legend([h1(1),h2(1),h3(1)],'r=10','r=15','r=20','location', 'southeast');
xlabel('Sample Size K');
ylabel('lamda1_bar');
axis([0 20000 0 20])