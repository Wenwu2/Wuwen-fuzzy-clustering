% function w_wflicm=wflicm(data_noise,v1,m,n,c,mc,e,ct)
% function[w_wflicm,w_wflicm_Vpc,w_wflicm_SA,w_wflicm_psnr]=wflicm(data_noise,v1,m,n,c,mc,e,ct);
function[I_WFLICM,WFLICM_Vpc,WFLICM_psnr,WFLICM_SA,WFLICM_Acc,WFLICM_Sen,WFLICM_Jaccard,WFLICM_Kappa]=WFLICM(data_noise,v1,m,n,c,mc,e,ct)
data_noise=double(data_noise);
data_noise=data_noise(3:m-2,3:n-2);
[m,n]=size(data_noise);
%初始化距离
d=zeros(m,n,c);
d1=zeros(m,n,c);
% 初始化隶属度
u=zeros(m,n,c);
t8=clock;
et=0.01;
t=0;
while e>0.0001 && ct<1000 %循环条件
    v=v1;
    % 算距离:样本data(i,j)到第k类的距离
    for k=1:c
        for  i=1:m
            for j=1:n
                d1(i,j,k)=(data_noise(i,j)-v1(k))^2+0.0001;
            end
        end
    end    
    % 算隶属度
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
    % 更新聚类中心
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for j=1:n
            for i=1:m
                tp1=tp1+u(i,j,k)^mc*data_noise(i,j);
                tp2=tp2+u(i,j,k)^mc;
            end
        end
        v1(k)=tp1/tp2; %聚类中心
    end
    % 终止条件
    temp=0.0;
    for k=1:c
        temp=temp+(v(k)-v1(k))^2;
    end
    if temp<0.0001
        e=0.0001;
    end
    ct=ct+1;
end

G=zeros(m,n,c);
Vj=rand(m,n);
Vbar=rand(m,n);
chi=rand(m,n); 

%算Vj
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
        Vj(i,j)= var(daa)/(ave^2+0.0001);
    end
end

%计算Vbar
for i=2:m-1
    for j=2:n-1
        tp1=0.0;
        for i1=-1:1
            for j1=-1:1
                tp1=tp1+ Vj(i+i1,j+j1);
            end
        end
        Vbar(i,j)=tp1/9;
    end
end

%求chi
for i=2:m-1
    for j=2:n-1
        tp1=0.0;
        for i1=-1:1
            for j1=-1:1
                tp1=tp1+exp(-( Vj(i+i1,j+j1)-Vbar(i+i1,j+j1)));
            end
        end
        chi(i,j)=exp(-( Vj(i,j)-Vbar(i,j)))/tp1;
    end
end

%求wgc
for i=2:m-1
    for j=2:n-1
        for i1=-1:1
            for j1=-1:1                
                if  Vj(i+i1,j+j1)<Vbar(i+i1,j+j1)
                    Wgc(i+i1,j+j1)=2+chi(i+i1,j+j1);
                else
                    Wgc(i+i1,j+j1)=2-chi(i+i1,j+j1);
                end
            end
        end
    end
end

%求Wij
for i=2:m-1
    for j=2:n-1
        for i1=-1:1
            for j1=-1:1
                Wij(i+i1,j+j1)=(1/((i1^2+j1^2)^0.5+1))*Wgc(i+i1,j+j1);
            end
        end
    end
end

while et>0.0001 && t<1000   %循环条件
    v=v1;
    % 算G
    for k=1:c
        for i=2:m-1
            for j=2:n-1
                tp1=0.0;
                for i1=-1:1
                    for j1=-1:1
                        if i1~=0||j1~=0
                            tp1=tp1+Wij(i+i1,j+j1)*((1-u(i+i1,j+j1,k))^mc*(data_noise(i+i1,j+j1)-v1(k))^2);
                        end
                    end
                end
                G(i,j,k)=tp1;
            end
        end
    end
    %算距离
    for k=1:c
        for i=1:m
            for j=1:n
                d(i,j,k)=(data_noise(i,j)-v1(k))^2+G(i,j,k)+0.0001;
            end
        end
    end    
    % 算隶属度
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
    % 更新聚类中心
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for i=1:m
            for j=1:n
                tp1=tp1+u(i,j,k)^mc*data_noise(i,j);
                tp2=tp2+u(i,j,k)^mc;
            end
        end
        v1(k)=tp1/tp2; %聚类中心
    end    
    % 终止条件
    temp=0.0;
    for k=1:c
        temp=temp+(v(k)-v1(k))^2;
    end
    if temp<0.0001
        et=0.0001;
    end
%     fprintf('迭代次数: %d ; Du: %f\n',t,temp);
    t=t+1;
end
disp(['WFLICM运行时间:',num2str(etime(clock,t8))]);
% 聚类
I_WFLICM=zeros(m,n);
if c==2
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)
                I_WFLICM(i,j)=0;
            else
                I_WFLICM(i,j)=227;
            end
        end
    end
elseif c==3
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)
                I_WFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)
                I_WFLICM(i,j)=127;
            else
                I_WFLICM(i,j)=255;
            end
        end
    end
elseif c==4
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)
                I_WFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)
                I_WFLICM(i,j)=80;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)
                I_WFLICM(i,j)=150;
            else
                I_WFLICM(i,j)=255;
            end
        end
    end
elseif c==5
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)&&u(i,j,1)>u(i,j,5)
                I_WFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)&&u(i,j,2)>u(i,j,5)
                I_WFLICM(i,j)=100;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)&&u(i,j,3)>u(i,j,5)
                I_WFLICM(i,j)=155;
            elseif u(i,j,4)>u(i,j,1)&& u(i,j,4)>u(i,j,2)&&u(i,j,4)>u(i,j,3)&&u(i,j,4)>u(i,j,5)
                I_WFLICM(i,j)=200;
            else
                I_WFLICM(i,j)=255;
            end
        end
    end
end
% I_WFLICM=I_WFLICM(3:m-2,3:n-2);
I=I_WFLICM;
[m,n]=size(I);
% [wflicm_Vpc,wflicm_SA,wflicm_psnr]=huafenzhibiao(u,m,n,I,c);
[WFLICM_Vpc,WFLICM_psnr,WFLICM_SA,WFLICM_Acc,WFLICM_Sen,WFLICM_Jaccard,WFLICM_Kappa]=huafenzhibiao(u,m,n,I,c);