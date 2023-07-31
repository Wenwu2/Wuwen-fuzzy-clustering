% function w_rflicm=rflicm(data_noise,v1,m,n,c,mc,e,ct)
% function[w_rflicm,w_rflicm_Vpc,w_rflicm_SA,w_rflicm_psnr]=rflicm(data_noise,v1,m,n,c,mc,e,ct);
function[I_RFLICM,RFLICM_Vpc,RFLICM_psnr,RFLICM_SA,RFLICM_Acc,RFLICM_Sen,RFLICM_Jaccard,RFLICM_Kappa]=RFLICM(data_noise,v1,m,n,c,mc,e,ct)
data_noise=double(data_noise);
data_noise=data_noise(3:m-2,3:n-2);
[m,n]=size(data_noise);
%��ʼ������
d=zeros(m,n,c);
d1=zeros(m,n,c);
% ��ʼ��������
u=zeros(m,n,c);
t6=clock;
et=0.01;
t=0;
while e>0.0001 && ct<1000 %ѭ������
    v=v1;
    % �����:����data(i,j)����k��ľ���
    for k=1:c
        for  i=1:m
            for j=1:n
                d1(i,j,k)=(data_noise(i,j)-v1(k))^2+0.0001;
            end
        end
    end
    % ��������
    for j=1:n
        for i=1:m
            tp1=0.0;
            for k=1:c
                tp1=tp1+(1/d1(i,j,k))^(1/(mc-1));
            end
            for k=1:c
                u(i,j,k)=(1/d1(i,j,k))^(1/(mc-1))/tp1;
            end
        end
    end
    % ���¾�������
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for j=1:n
            for i=1:m
                tp1=tp1+u(i,j,k)^mc*data_noise(i,j);
                tp2=tp2+u(i,j,k)^mc;
            end
        end
        v1(k)=tp1/tp2;            %��������
    end
    % ��ֹ����
    temp=0.0;
    for k=1:c
        temp=temp+(v(k)-v1(k))^2;
    end
    if   temp<0.0001
        e=0.0001;
    end
    ct=ct+1;
end

G=zeros(m,n,c);
Ck=rand(m,n);
Cbar=rand(m,n);

%��Ck
for i=2:m-1
    for j=2:n-1
        tp1=0.0;
        daa=zeros(1,9);
        r=1;
        for i1=-1:1
            for j1=-1:1
                tp1=tp1+(data_noise(i+i1,j+j1));
                daa(1,r)=data_noise(i+i1,j+j1);
                r=r+1;
            end
        end
        ave=tp1/9;
        Ck(i,j)= var(daa)/(ave^2+0.0001);
    end
end

%����Cbar
for i=2:m-1
    for j=2:n-1
        tp1=0.0;
        for i1=-1:1
            for j1=-1:1
                tp1=tp1+ Ck(i+i1,j+j1);
            end
        end
        Cbar(i,j)=tp1/9;
    end
end

%��Wik
for i=2:m-1
    for j=2:n-1
        for i1=-1:1
            for j1=-1:1
                if Ck(i+i1,j+j1)>Cbar(i+i1,j+j1)
                    W(i+i1,j+j1)=1/(2+min((Ck(i+i1,j+j1)/(Ck(i,j)+0.0001))^2,(Ck(i,j)/(Ck(i+i1,j+j1)+0.0001))^2));
                else
                    W(i+i1,j+j1)=1/(2-min((Ck(i+i1,j+j1)/(Ck(i,j)+0.0001))^2,(Ck(i,j)/(Ck(i+i1,j+j1)+0.0001))^2));
                end
            end
        end
    end
end

while et>0.0001 && t<1000   %ѭ������
    v=v1;
    % ��G��
    for k=1:c
        for i=2:m-1
            for j=2:n-1
                tp1=0.0;
                for i1=-1:1
                    for j1=-1:1
                        if i1~=0||j1~=0
                            tp1=tp1+W(i+i1,j+j1)*((1-u(i+i1,j+j1,k))^mc*(data_noise(i+i1,j+j1)-v1(k))^2)/((i1^2+j1^2)^0.5+1);
                        end
                    end
                end
                G(i,j,k)=tp1;
            end
        end
    end
    %�����
    for k=1:c
        for i=1:m
            for j=1:n
                d(i,j,k)=(data_noise(i,j)-v1(k))^2+G(i,j,k)+0.0001;
            end
        end
    end
    % ��������
    for i=1:m
        for j=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+d(i,j,k)^(-1/(mc-1));
            end
            for k=1:c
                u(i,j,k)=d(i,j,k)^(-1/(mc-1))/tp1;
            end
        end
    end
    % ���¾�������
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for i=1:m
            for j=1:n
                tp1=tp1+u(i,j,k)^mc*data_noise(i,j);
                tp2=tp2+u(i,j,k)^mc;
            end
        end
        v1(k)=tp1/tp2;            %��������
    end
    % ��ֹ����
    temp=0.0;
    for k=1:c
        temp=temp+(v(k)-v1(k))^2;
    end
    if temp<0.0001
        et=0.0001;
    end
%     fprintf('��������: %d ; Du: %f\n',t,temp);
    t=t+1;
end
disp(['RFLICM����ʱ�䣺',num2str(etime(clock,t6))]);
% ����
I_RFLICM=zeros(m,n);
if c==2
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)
                I_RFLICM(i,j)=0;
            else
                I_RFLICM(i,j)=227;
            end
        end
    end
elseif c==3
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)
                I_RFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)
                I_RFLICM(i,j)=127;
            else
                I_RFLICM(i,j)=255;
            end
        end
    end
elseif c==4
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)
                I_RFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)
                I_RFLICM(i,j)=80;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)
                I_RFLICM(i,j)=150;
            else
                I_RFLICM(i,j)=255;
            end
        end
    end
elseif c==5
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)&&u(i,j,1)>u(i,j,5)
                I_RFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)&&u(i,j,2)>u(i,j,5)
                I_RFLICM(i,j)=100;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)&&u(i,j,3)>u(i,j,5)
                I_RFLICM(i,j)=155;
            elseif u(i,j,4)>u(i,j,1)&& u(i,j,4)>u(i,j,2)&&u(i,j,4)>u(i,j,3)&&u(i,j,4)>u(i,j,5)
                I_RFLICM(i,j)=200;
            else
                I_RFLICM(i,j)=255;
            end
        end
    end
end
% I_RFLICM=I_RFLICM(3:m-2,3:n-2);
I=I_RFLICM;
[m,n]=size(I);
% [rflicm_Vpc,rflicm_SA,rflicm_psnr]=huafenzhibiao(u,m,n,I,c);
[RFLICM_Vpc,RFLICM_psnr,RFLICM_SA,RFLICM_Acc,RFLICM_Sen,RFLICM_Jaccard,RFLICM_Kappa]=huafenzhibiao(u,m,n,I,c);