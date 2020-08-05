clear;
T0=1;%��������
MPSK=2;%BPSK
K=160;%simple size
M=300;%FFT size
Gamma=4;%over-simpling factor
fs=Gamma/T0;%������
t=0:1/fs:T0-1/fs;%ʱ������
fc=2/T0; %�ز�Ƶ��  Ҳ��һ����Ԫ���ٸ���������
c=sqrt(2)*exp(1i*2*pi*fc*t);%�ز��ź�
nsymbol=K/fs;%ÿ��������µķ��ͷ�����
theta=25;%menta carlo number
snr_dB=-20:2:10;
msg=randi([0 MPSK-1],1,nsymbol); %���ɻ�������       
msgmod=pskmod(msg,MPSK).'; %����B-PSK����
tx=real(msgmod*c);%�ز�����
tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'��ÿһ����һ����Ԫ����Ĳ�����,��չ��Ϊһ�� 

r=10;%quantization level
r_set=1/r:1/r:1;%��������
for indx=1:length(snr_dB)
    lamda1=zeros(1,theta);
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
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
        lamda1(indx2)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave_BPSK(indx)=sum(lamda1)/length(lamda1);
    
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
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
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%�ҵ������ȼ�
        end
        Lx=get_LaplacianMatrix(r,r_level);%�õ�laplacian ����
        [~,lamda]=eig(Lx);%��������ֵ
        [not_sort,~]=max(lamda);%��ȡ����ֵ
        lamda_sort=sort(not_sort);%����ֵ����
        lamda1(indx2)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave_AWGN(indx)=sum(lamda1)/length(lamda1);
end
figure(1)
h11=plot(snr_dB,lamda_ave_BPSK,'ro-');
hold on;
h12=plot(snr_dB,lamda_ave_AWGN,'ro--');

r=15;%quantization level
r_set=1/r:1/r:1;%��������
for indx=1:length(snr_dB)
    lamda1=zeros(1,theta);
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
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
        lamda1(indx2)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave_BPSK(indx)=sum(lamda1)/length(lamda1);
    
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
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
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%�ҵ������ȼ�
        end
        Lx=get_LaplacianMatrix(r,r_level);%�õ�laplacian ����
        [~,lamda]=eig(Lx);%��������ֵ
        [not_sort,~]=max(lamda);%��ȡ����ֵ
        lamda_sort=sort(not_sort);%����ֵ����
        lamda1(indx2)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave_AWGN(indx)=sum(lamda1)/length(lamda1);
end
h21=plot(snr_dB,lamda_ave_BPSK,'g+-');
hold on;
h22=plot(snr_dB,lamda_ave_AWGN,'g+--');

r=20;%quantization level
r_set=1/r:1/r:1;%��������
for indx=1:length(snr_dB)
    lamda1=zeros(1,theta);
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
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
        lamda1(indx2)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave_BPSK(indx)=sum(lamda1)/length(lamda1);
    
    for indx2=1:1:theta
        %rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
        rx=noisegen(tx1,snr_dB(indx),T0,fs);%�����˹������
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
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%�ҵ������ȼ�
        end
        Lx=get_LaplacianMatrix(r,r_level);%�õ�laplacian ����
        [~,lamda]=eig(Lx);%��������ֵ
        [not_sort,~]=max(lamda);%��ȡ����ֵ
        lamda_sort=sort(not_sort);%����ֵ����
        lamda1(indx2)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
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