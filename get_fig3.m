clear;
T0=1;%��������
MPSK=2;%BPSK
%r=[10 15 20];%quantization level
K=100:200:20000;%simple size
M=300;%FFT size
Gamma=50;%over-simpling factor
fs=Gamma/T0;%������
t=0:1/fs:T0-1/fs;%ʱ������
fc=2/T0; %�ز�Ƶ��  Ҳ��һ����Ԫ���ٸ���������
c=sqrt(2)*exp(1i*2*pi*fc*t);%�ز��ź�
theta=1000;%menta carlo number
snr_dB=-1000;
r=10;
r_set=1/r:1/r:1;%��������
for kk=1:1:length(K)
    lamda1=zeros(1,theta);
    for jj=1:1:theta
        nsymbol=round(K(kk)/Gamma);%ÿ��������µķ��ͷ�����
        msg=randi([0 MPSK-1],1,nsymbol); %���ɻ�������       
        msgmod=pskmod(msg,MPSK).'; %����B-PSK����
        tx=real(msgmod*c);%�ز�����
        tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'��ÿһ����һ����Ԫ����Ĳ�����,��չ��Ϊһ��      
        rx=noisegen(tx1,snr_dB,T0,fs);%�����˹������
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
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%�ҵ������ȼ�
        end
        Lx=get_LaplacianMatrix(r,r_level);%�õ�laplacian ����
        [~,lamda]=eig(Lx);%��������ֵ
        [not_sort,~]=max(lamda);%��ȡ����ֵ
        lamda_sort=sort(not_sort);%����ֵ����
        lamda1(jj)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave1=sum(lamda1)/length(lamda1);
    h1=plot(K(kk),lamda_ave1,'r.');
    hold on;
end

r=15;
r_set=1/r:1/r:1;%��������
for kk=1:1:length(K)
    lamda1=zeros(1,theta);
    for jj=1:1:theta
        nsymbol=round(K(kk)/Gamma);%ÿ��������µķ��ͷ�����
        msg=randi([0 MPSK-1],1,nsymbol); %���ɻ�������       
        msgmod=pskmod(msg,MPSK).'; %����B-PSK����
        tx=real(msgmod*c);%�ز�����
        tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'��ÿһ����һ����Ԫ����Ĳ�����,��չ��Ϊһ��      
        rx=noisegen(tx1,snr_dB,T0,fs);%�����˹������
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
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%�ҵ������ȼ�
        end
        Lx=get_LaplacianMatrix(r,r_level);%�õ�laplacian ����
        [~,lamda]=eig(Lx);%��������ֵ
        [not_sort,~]=max(lamda);%��ȡ����ֵ
        lamda_sort=sort(not_sort);%����ֵ����
        lamda1(jj)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave2=sum(lamda1)/length(lamda1);
    h2=plot(K(kk),lamda_ave2,'b.');
    hold on;
end

r=20;
r_set=1/r:1/r:1;%��������
for kk=1:1:length(K)
    lamda1=zeros(1,theta);
    for jj=1:1:theta
        nsymbol=round(K(kk)/Gamma);%ÿ��������µķ��ͷ�����
        msg=randi([0 MPSK-1],1,nsymbol); %���ɻ�������       
        msgmod=pskmod(msg,MPSK).'; %����B-PSK����
        tx=real(msgmod*c);%�ز�����
        tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'��ÿһ����һ����Ԫ����Ĳ�����,��չ��Ϊһ��      
        rx=noisegen(tx1,snr_dB,T0,fs);%�����˹������
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
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%�ҵ������ȼ�
        end
        Lx=get_LaplacianMatrix(r,r_level);%�õ�laplacian ����
        [~,lamda]=eig(Lx);%��������ֵ
        [not_sort,~]=max(lamda);%��ȡ����ֵ
        lamda_sort=sort(not_sort);%����ֵ����
        lamda1(jj)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
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