clear;
T0=1;%��������
MPSK=2;%BPSK
r=10;%quantization level
K=300;%simple size
M=300;%FFT size
Gamma=50;%over-simpling factor
fs=Gamma/T0;%������
t=0:1/fs:T0-1/fs;%ʱ������
fc=2/T0; %�ز�Ƶ��  Ҳ��һ����Ԫ���ٸ���������
c=sqrt(2)*exp(1i*2*pi*fc*t);%�ز��ź�
nsymbol=K/Gamma;%ÿ��������µķ��ͷ�����
theta=1;%menta carlo number
r_set=1/r:1/r:1;%��������
snr_dB=[-3 -7 -11 -1000];
msg=randi([0 MPSK-1],1,nsymbol); %���ɻ�������       
msgmod=pskmod(msg,MPSK).'; %����B-PSK����
tx=real(msgmod*c);%�ز�����
tx1=reshape(tx.',1,length(msgmod)*length(c));   %tx'��ÿһ����һ����Ԫ����Ĳ�����,��չ��Ϊһ��       
for indx=1:length(snr_dB)
    lamda0=zeros(1,theta);
    for indx2=1:1:theta
        rx=noisegen_err(tx1,snr_dB(indx),T0,fs);%�����˹������
        rxy=abs(fft(rx,300));%fft
        figure(1)%��ʾfftͼ��
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
            [~,r_level(mm)]=min(abs(Ux(mm)-r_set));%�ҵ������ȼ�
        end
    figure(2);%��ʾ�����ȼ�ͼ�� ͬ�����е�Fig.1
    subplot(4,1,indx)
    plot(r_level);
    title(['SNR=',num2str(snr_dB(indx)),'dB']);
    xlabel('M','position',[320 -20]);
    ylabel('Qx(m)');
        Lx=get_LaplacianMatrix(r,r_level);%�õ�laplacian ����
        [~,lamda]=eig(Lx);%��������ֵ
        [not_sort,~]=max(lamda);%��ȡ����ֵ
        lamda_sort=sort(not_sort);%����ֵ����
        lamda0(indx2)=lamda_sort(end-1);%�ҵ��ڶ�������ֵ
    end
    lamda_ave=sum(lamda0)/length(lamda0);
end
