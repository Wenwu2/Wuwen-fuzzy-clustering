function[I_PFLICM,PFLICM_Vpc,PFLICM_psnr,PFLICM_SA,PFLICM_Acc,PFLICM_Sen,PFLICM_Jaccard,PFLICM_Kappa]=PFLICM(data_noise,v1,m,n,c,e,ct)

p=2;q=2;
a=6;b=1;
d=zeros(m,n,c);
[u1,v1]=fcm(m,n,c,v1,data_noise,p);
u=u1;
tt=clock;
% u=zeros(m,n,c);
t=zeros(m,n,c);
G=zeros(m,n,c);
R=zeros(m,n,c);
for k=1:c
    tp1=0.0;
	tp2=0.0;
    for x=1:m
        for y=1:n
            tp1=tp1+(u(x,y,k).^p)*(v1(k)-data_noise(x,y))^2;
            tp2=tp2+u(x,y,k)^p;
        end
    end
    r(k)=2*tp1/(tp2+0.00001);
end
    
while e>0.0001 && ct<200
%     disp(['itteration ' num2str(ct)]);%打印迭代次数
    v=v1;
   for k=1:c
      for x=1:m
          for y=1:n
              d(x,y,k)=(data_noise(x,y)-v(k))^2+0.00001;
          end
      end
   end
   
   for k=1:c
    for x=1:m
        for y=1:n
            tp1=0.0;
            tp2=0.0;
            for i=x-1:x+1
                for j=y-1:y+1
                  if(i>=1&&j>=1&&i<=m&&j<=n)
                    tp1=tp1+((1-u(i,j,k))^p*d(i,j,k)/(((x-i)^2+(y-j)^2)^0.5+1)+0.000001);
%                     tp2=tp2+u(i,j,k)^2/(((x-i)^2+(y-j)^2)^0.5+1);
                   end
                end
            end
            G(x,y,k)=tp1; 
%             R(x,y,k)=tp2;
        end
    end
    
    
   end
    
    
    for x=1:m
        for y=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+(d(x,y,k)+G(x,y,k))^(-1/(p-1));
            end
            for k=1:c
                u(x,y,k)=(d(x,y,k)+G(x,y,k))^(-1/(p-1))/(tp1);%隶属度矩阵
            end
        
        end
    end
      
  for k=1:c
       for x=1:m
           for y=1:n
               t(x,y,k)=1/(1+(b*d(x,y,k)/(r(k)))^(1/(q-1)));
           end
       end
  end 
    
   for k=1:c
        tp1=0.0;
        tp2=0.0;
        for x=1:m
            for y=1:n
                tp1=tp1+(a*u(x,y,k)^p+b*t(x,y,k)^q)*data_noise(x,y);
                tp2=tp2+(a*u(x,y,k)^p+b*t(x,y,k)^q);
            end
        end
        v1(k)=tp1/(tp2);%中心进行迭代 找到最优中心
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
disp(['PFLICM运行时间:',num2str(etime(clock,tt))]);
I_PFLICM=zeros(m,n);
for x=1:m
    for y=1:n
        if c==4
            if u(x,y,1)>u(x,y,2) && u(x,y,1)>u(x,y,3) && u(x,y,1)> u(x,y,4)
                I_PFLICM(x,y)=80;%
            end
            if u(x,y,2)>u(x,y,1) && u(x,y,2)>u(x,y,3) && u(x,y,2)> u(x,y,4)
                I_PFLICM(x,y)=150;%
            end
            if u(x,y,3)>u(x,y,1) && u(x,y,3)>u(x,y,2) && u(x,y,3)>u(x,y,4)
                I_PFLICM(x,y)=255;%
            end
            if u(x,y,4)>u(x,y,1)&& u(x,y,4)>u(x,y,2) && u(x,y,4)>u(x,y,3)
                I_PFLICM(x,y)=0;
            end
        end
        if c==3
            if u(x,y,1)>u(x,y,2) && u(x,y,1)>u(x,y,3)
                I_PFLICM(x,y)=0;
            end
            if u(x,y,2)>u(x,y,1) && u(x,y,2)>u(x,y,3)
                I_PFLICM(x,y)=127;
            end
            if u(x,y,3)>u(x,y,1) && u(x,y,3)>u(x,y,2)
                I_PFLICM(x,y)=255;
            end
        end
        if c==2
            if u(x,y,1)>u(x,y,2)
                I_PFLICM(x,y)=0;
            else
                I_PFLICM(x,y)=227;
            end
        end
    end
end
I_PFLICM=I_PFLICM(3:m-2,3:n-2);
I=I_PFLICM;
[m,n]=size(I);
% [wflicm_Vpc,wflicm_SA,wflicm_psnr]=huafenzhibiao(u,m,n,I,c);
[PFLICM_Vpc,PFLICM_psnr,PFLICM_SA,PFLICM_Acc,PFLICM_Sen,PFLICM_Jaccard,PFLICM_Kappa]=huafenzhibiao(u,m,n,I,c);
     