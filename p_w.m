function [I_p_w,p_w_Vpc,p_w_psnr,p_w_SA,p_w_Acc,p_w_Sen,p_w_Jaccard,p_w_Kappa]=p_w(data_noise,v1,m,n,c,e,ct)

K=1;
p=2;q=2;
a=6;b=1;
r1=zeros(c,1);
r1(1)=1;r1(2)=1;r1(3)=1;r1(4)=1;
tt=clock;
d1=zeros(m,n,c);
G1=zeros(m,n,c);
G2=zeros(m,n,c);
G3=zeros(m,n,c);
G4=zeros(m,n,c);
[u1]=fcm(m,n,c,v1,data_noise,p);%初始化u
[t1]=pcm(m,n,c,v1,data_noise,q,u1);%初始化t
% t1=u1;
while e>0.0001 && ct<1000
    v=v1;
    u=u1;
    t=t1;
    r=r1;
    for k=1:c
        for x=1:m
            for y=1:n
                tp1=0.0;
                tp2=0.0;
                tp3=0.0;
                for i=x-1:x+1
                    for j=y-1:y+1
                        if (i>=1&&j>=1&&i<=m&&j<=n&&i~=x&&j~=y)
                            tp1=tp1+(a*(1-u(i,j,k))^p+b*(1-t(i,j,k))^q)/(((x-i)^2+(y-j)^2)^0.5+1);
                            tp2=tp2+(data_noise(i,j)-v(k))^2*(a*(1-u(i,j,k))^p+b*(1-t(i,j,k))^q)/(((x-i)^2+(y-j)^2)^0.5+1);
                            tp3=tp3+data_noise(i,j)*(a*(1-u(i,j,k))^p+b*(1-t(i,j,k))^q)/(((x-i)^2+(y-j)^2)^0.5+1);   
                        end
                    end
                end  
                G1(x,y,k)=tp1;
                G2(x,y,k)=tp2;
                G3(x,y,k)=tp3;
            end
        end
    end
    
     for k=1:c
        for x=1:m
            for y=1:n
                tp4=0.0;
                for i=x-1:x+1
                    for j=y-1:y+1
                        if (i>=1&&j>=1&&i<=m&&j<=n&&i~=x&&j~=y)
                            tp4=tp4+(data_noise(i,j)-v(k))^2+G2(i,j,k);
                        end
                    end
                end
                G4(x,y,k)=tp4;
            end
        end
    end
     for x=1:m
          for y=1:n
              for k=1:c
                   d1(x,y,k)=(data_noise(x,y)-v(k))^2+G1(x,y,k)*G4(x,y,k)+0.0001;
              end
          end
     end
     for x=1:m
        for y=1:n
            tep=0.0;
            for k=1:c
                tep=tep+d1(x,y,k)^(-1/(p-1));
            end
            for k=1:c
                u1(x,y,k)=d1(x,y,k)^(-1/(p-1))/(tep);                             
            end
        end
     end
     
     %计算可能t
     for k=1:c
         for x=1:m
%              tpp=0.0;
             for y=1:n
%                  tpp=tpp+b*q*(1-t(x,y,k))*((data_noise(x,y)-v(k))^2+G1(x,y,k)*G4(x,y,k));
 %                  t1(x,y,k)=(1+(tpp/(r(k)*(q-1)))^(1/(q-1)))^(-1);%q-1
                 t1(x,y,k)=(1+(b*((data_noise(x,y)-v(k))^2+G1(x,y,k)*G4(x,y,k))/r(k))^(1/(q-1)))^(-1);%q
             end
         end
     end
    
    %更新聚类中心
     for k=1:c
         tep1=0.0;
         tep2=0.0;
         for x=1:m
             for y=1:n
                 tep1=tep1+(a*(u1(x,y,k))^p+b*(t1(x,y,k))^q)*(data_noise(x,y)+G3(x,y,k));
                 tep2=tep2+(a*(u1(x,y,k))^p+b*(t1(x,y,k))^q)*(1+G1(x,y,k));
             end
         end
         v1(k)=tep1/(tep2+0.0001);%聚类中心
     end

    for k=1:c
         tep3=0.0;
         tep4=0.0;
         for x=1:m
             for y=1:n
                 tep3=tep3+u1(x,y,k)^p*(data_noise(x,y)-v1(k))^2;
                 tep4=tep4+u1(x,y,k)^p;
             end
         end
         r1(k)=K*tep3/tep4;
     end

    %终止条件
    temp=0.0;
    for k=1:c
        temp=temp+(v(k)-v1(k))^2;
    end
    if temp<0.0001
        e=0.0001;
    end
%     fprintf('迭代次数: %d ; Du: %f\n',ct,temp);
    ct=ct+1;
end
disp(['p_w运行时间:',num2str(etime(clock,tt))]);
u=u1;
I_p_w=zeros(m,n);
for x=1:m
    for y=1:n
        if c==4
            if u(x,y,1)>u(x,y,2)&&u(x,y,1)>u(x,y,3)&&u(x,y,1)>u(x,y,4)
                I_p_w(x,y)=150;
            end
            if u(x,y,2)>u(x,y,1)&&u(x,y,2)>u(x,y,3)&&u(x,y,2)>u(x,y,4)
                I_p_w(x,y)=0;
            end
            if u(x,y,3)>u(x,y,1)&&u(x,y,3)>u(x,y,2)&&u(x,y,3)>u(x,y,4)
                I_p_w(x,y)=80;
            end
            if u(x,y,4)>u(x,y,1)&&u(x,y,4)>u(x,y,2)&&u(x,y,4)>u(x,y,3)
                I_p_w(x,y)=255;
            end
        end
        if c==3
            if u(x,y,1)>u(x,y,2)&&u(x,y,1)>u(x,y,3)
                I_p_w(x,y)=0;
            end
            if u(x,y,2)>u(x,y,1)&&u(x,y,2)>u(x,y,3)
                I_p_w(x,y)=127;
            end
            if u(x,y,3)>u(x,y,1)&&u(x,y,3)>u(x,y,2)
                I_p_w(x,y)=255;
            end 
        end
        if c==2
            if u(x,y,1)>u(x,y,2)
                I_p_w(x,y)=0;
            else
                I_p_w(x,y)=227;
            end
         end         
    end
end
I_p_w=I_p_w(3:m-2,3:n-2);
I=I_p_w;
[m,n]=size(I);
[p_w_Vpc,p_w_psnr,p_w_SA,p_w_Acc,p_w_Sen,p_w_Jaccard,p_w_Kappa]=huafenzhibiao(u,m,n,I,c);
