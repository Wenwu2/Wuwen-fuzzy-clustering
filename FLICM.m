% function w_flicm=flicm(data_noise,v1,m,n,c,mc,e,ct)
% function[w_flicm,fcm_flicm_Vpc,fcm_flicm_SA,fcm_flicm_psnr]=flicm(data_noise,v1,m,n,c,mc,e,ct);
function[I_FLICM,FLICM_Vpc,FLICM_psnr,FLICM_SA,FLICM_Acc,FLICM_Sen,FLICM_Jaccard,FLICM_Kappa]=FLICM(data_noise,v1,m,n,c,mc,e,ct)

% [uu]=fcm(data_noise,m,n,c,v1);
% u1=uu;
data_noise=double(data_noise);
data_noise=data_noise(3:m-2,3:n-2);
[m,n]=size(data_noise);
et=0.1;
t=0; 
d=zeros(m,n,c);
u1=zeros(m,n,c);
while e>0.0001 && ct<500
    v=v1;
    for k=1:c
        for x=1:m 
            for y=1:n
                d(x,y,k)=(data_noise(x,y)-v(k))^2+0.0001; 
            end
        end
    end
    for x=1:m
        for y=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+d(x,y,k)^(-1/(mc-1));
            end
            for k=1:c
                u1(x,y,k)=d(x,y,k)^(-1/(mc-1))/tp1;
            end
        end
    end
    %这里的for是用来计算聚类中心的
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for x=1:m
            for y=1:n
                tp1=tp1+u1(x,y,k)^mc*data_noise(x,y);
                tp2=tp2+u1(x,y,k);
            end
        end
        v1(k)=tp1/tp2;
    end
    temp=0.0;
    for k=1:c
        temp=temp+(v(k)-v1(k))^2;
    end
    if temp<0.0001
        e=0.0001;
    end
    ct=ct+1;
end

d1=zeros(m,n,c);
G=zeros(m,n,c);
t5=clock;
while et>0.0001 && t<1000
    v=v1;
    u=u1;
    %计算距离
    for k=1:c
        for x=1:m
            for y=1:n
                tp=0.0;
                for i=x-1:x+1
                    for j=y-1:y+1
                        if (i>=1&&j>=1&&i<=m&&j<=n&&i~=x&&j~=y)
                            tp=tp+((1-u(i,j,k))^mc)*(data_noise(i,j)-v(k))^2/(((x-i)^2+(y-j)^2)^0.5+1);                           
                        end
                    end
                end
                G(x,y,k)=tp;
            end
        end
    end
%     for k=1:c
%         for x=2:m-1
%             for y=2:n-1
%                 tp=0.0;
%                 for x1=-1:1
%                     for y1=-1:1
%                         if(x1~=0||y1~=0)
%                             tp=tp+((1-u(x+x1,y+y1,k))^mc)*(data_noise(x+x1,y+y1)-v(k))^2/(1+sqrt(x1^2+y1^2));
%                         end
%                     end
%                 end
%                 G(x,y,k)=tp;
%             end
%         end
%     end
    for k=1:c
        for x=1:m
            for y=1:n
                d1(x,y,k)=(data_noise(x,y)-v(k))^2+G(x,y,k)+0.0001;
            end
        end
    end
        
    %更新隶属度矩阵
    for x=1:m
        for y=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+d1(x,y,k)^(-1/(mc-1));
            end
            for k=1:c
                u1(x,y,k)=d1(x,y,k)^(-1/(mc-1))/tp1;%隶属度
            end
        end
    end
    %更新聚类中心
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for x=2:m-1
            for y=2:n-1
                tp1=tp1+u1(x,y,k)^mc*data_noise(x,y);
                tp2=tp2+u1(x,y,k)^mc;
            end
        end
        v1(k)=tp1/(tp2+0.0001);%聚类中心
    end
    %终止条件
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
disp(['FLICM运行时间：',num2str(etime(clock,t5))]);
u=u1;
I_FLICM=zeros(m,n);
for x=1:m
    for y=1:n
        if c==4
            if u(x,y,1)>u(x,y,2)&&u(x,y,1)>u(x,y,3)&&u(x,y,1)>u(x,y,4)
                I_FLICM(x,y)=255;
            end
            if u(x,y,2)>u(x,y,1)&&u(x,y,2)>u(x,y,3)&&u(x,y,2)>u(x,y,4)
                I_FLICM(x,y)=150;
            end
            if u(x,y,3)>u(x,y,1)&&u(x,y,3)>u(x,y,2)&&u(x,y,3)>u(x,y,4)
                I_FLICM(x,y)=80;
            end
            if u(x,y,4)>u(x,y,1)&&u(x,y,4)>u(x,y,2)&&u(x,y,4)>u(x,y,3)
                I_FLICM(x,y)=0;
            end
         end
         if c==3
             if u(x,y,1)>u(x,y,2)&&u(x,y,1)>u(x,y,3)
                 I_FLICM(x,y)=0;
             end
             if u(x,y,2)>u(x,y,1)&&u(x,y,2)>u(x,y,3) 
                 I_FLICM(x,y)=127;
             end
             if u(x,y,3)>u(x,y,1)&&u(x,y,3)>u(x,y,2) 
                 I_FLICM(x,y)=255;
             end
         end
         if c==2
             if u(x,y,1)>u(x,y,2)
                 I_FLICM(x,y)=0;
             else
                 I_FLICM(x,y)=227;
             end
         end
    end
end
% I_FLICM=I_FLICM(3:m-2,3:n-2);
I=I_FLICM;
[m,n]=size(I);
% [flicm_Vpc,flicm_SA,flicm_psnr]=huafenzhibiao(u,m,n,I,c);
[FLICM_Vpc,FLICM_psnr,FLICM_SA,FLICM_Acc,FLICM_Sen,FLICM_Jaccard,FLICM_Kappa]=huafenzhibiao(u,m,n,I,c);
