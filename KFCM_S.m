function [I_KFCM_S,KFCM_S_Vpc,KFCM_S_psnr,KFCM_S_SA,KFCM_S_Acc,KFCM_S_Sen,KFCM_S_Jaccard,KFCM_S_Kappa]=KFCM_S(data_noise,v1,m,n,c,mc,e,ct)
% function I_KFCM_S=KFCM_S(data,v1,m,n,c,mc,e,ct)
data_noise=double(data_noise);
data_noise=data_noise(3:m-2,3:n-2);
[m,n]=size(data_noise);
alpha=3.8;

% 初始化隶属度
u=zeros(m,n,c);
K1=rand(m,n,c);
K2=rand(m,n,c);
t1=clock;
%算Pbar
  tp1=0.0;
  for i=1:m
      for j=1:n
          tp1=tp1+data_noise(i,j);
      end
  end
  Pbar=tp1/(m*n);
  %算da
  for i=1:m
      for j=1:n
          da(i,j)=abs(data_noise(i,j)-Pbar)^2;
      end
  end
  %算dbar
  tp1=0.0;
  for i=1:m
      for j=1:n
          tp1=tp1+da(i,j);
      end
  end
  dbar=tp1/(m*n);
  
  %算sigma
 tp1=0.0;
  for i=1:m
      for j=1:n
           tp1=tp1+(da(i,j)- dbar)^2;
      end
  end
  sigma=sqrt(tp1/(m*n-1));


while e>0.0001 && ct<1000     %终止条件
    v=v1;
    
 %算K
  for k=1:c
      for i=1:m
          for j=1:n
              K(i,j,k)=exp(-(data_noise(i,j)-v1(k))^2/sigma);
          end
      end
  end
    
  
    % 算隶属度
    for k=1:c
        for j=2:n-1
            for i=2:m-1
                tp1=0.0;
                for i1=-1:1
                    for j1=-1:1
                        tp1=tp1+(1-K(i+i1,j+j1,k))^mc;
                    end
                end
                K1(i,j,k)=alpha/9*tp1;
            end
        end
    end
    
    for i=1:m
        for j=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+((1-K(i,j,k))+K1(i,j,k)+0.0001)^(-1/(mc-1));
            end
            for k=1:c
                u(i,j,k)=((1-K(i,j,k))+K1(i,j,k)+0.0001)^(-1/(mc-1))/tp1;
            end
        end
    end
               
   
    
    % 更新聚类中心
    
     for k=1:c
        for j=2:n-1
            for i=2:m-1
                tp1=0.0;
                for i1=-1:1
                    for j1=-1:1
                        tp1=tp1+K(i+i1,j+j1,k);
                    end
                end
                K2(i,j,k)=alpha/9*tp1;
            end
        end
    end
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for i=1:m
            for j=1:n
                tp1=tp1+u(i,j,k)^mc*(K(i,j,k)+K2(i,j,k))*data_noise(i,j);
                tp2=tp2+u(i,j,k)^mc*(K(i,j,k)+K2(i,j,k));
            end
        end
        v1(k)=tp1/tp2;            %聚类中心
    end
    
    % 终止条件
   temp=0.0;
   for k=1:c
         temp=temp+(v(k)-v1(k))^2;
   end
   if   temp<0.0001
        e=0.0001;
   end
ct=ct+1;
end
disp(['KFCM_S运行时间',num2str(etime(clock,t1))]);
% 聚类
I_KFCM_S=zeros(m,n);
if c==2
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)
                I_KFCM_S(i,j)=0;
            else
                I_KFCM_S(i,j)=227;
            end
        end
    end
elseif c==3
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)
                I_KFCM_S(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)
                I_KFCM_S(i,j)=127;
            else
                I_KFCM_S(i,j)=255;
            end
        end
    end
elseif c==4
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)
                I_KFCM_S(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)
                I_KFCM_S(i,j)=80;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)
                I_KFCM_S(i,j)=150;
            else
                I_KFCM_S(i,j)=255;
            end
        end
    end
elseif c==5
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)&&u(i,j,1)>u(i,j,5)
                I_KFCM_S(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)&&u(i,j,2)>u(i,j,5)
                I_KFCM_S(i,j)=100;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)&&u(i,j,3)>u(i,j,5)
                I_KFCM_S(i,j)=155;
            elseif u(i,j,4)>u(i,j,1)&& u(i,j,4)>u(i,j,2)&&u(i,j,4)>u(i,j,3)&&u(i,j,4)>u(i,j,5)
                I_KFCM_S(i,j)=200;
            else
                I_KFCM_S(i,j)=255;
            end
        end
    end
end
I=I_KFCM_S;

[KFCM_S_Vpc,KFCM_S_psnr,KFCM_S_SA,KFCM_S_Acc,KFCM_S_Sen,KFCM_S_Jaccard,KFCM_S_Kappa]=huafenzhibiao(u,m,n,I,c);

