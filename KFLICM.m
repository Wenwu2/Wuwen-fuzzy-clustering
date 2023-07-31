 function [I_KFLICM,KFLICM_Vpc,KFLICM_psnr,KFLICM_SA,KFLICM_Acc,KFLICM_Sen,KFLICM_Jaccard,KFLICM_Kappa]=KFLICM(data_noise,v1,m,n,c,mc,e,ct)
%  function [I_KFLICM,KFLICM_Vpc,KFLICM_SA]=KFLICM(data,v1,m,n,c,mc,e,ct)
% function I_KFLICM=KFLICM(data,v1,m,n,c,mc,e,ct)
data_noise=double(data_noise);
data_noise=data_noise(3:m-2,3:n-2);
[m,n]=size(data_noise);
et=0.01;
t=0;

%��ʼ������
d=zeros(m,n,c);
% ��ʼ��������
u=zeros(m,n,c);

 
t1=clock;
while e>0.0001 && ct<1000     %ѭ������
    v=v1;
    
   % �����:����data(i,j)����k��ľ���
    for k=1:c
        for  i=1:m
            for j=1:n
                d(i,j,k)=(data_noise(i,j)-v1(k))^2+0.0001;
            end
        end
    end
    
    % ��������
    for j=1:n
        for i=1:m
            tp1=0.0;
            for k=1:c
                tp1=tp1+(1/d(i,j,k))^(1/(mc-1));
            end
            for k=1:c
                u(i,j,k)=(1/d(i,j,k))^(1/(mc-1))/tp1;
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
K=zeros(m,n,c);


 %��Pbar
  tp1=0.0;
  for i=1:m
      for j=1:n
          tp1=tp1+data_noise(i,j);
      end
  end
  Pbar=tp1/(m*n);
  %��da
  for i=1:m
      for j=1:n
          da(i,j)=abs(data_noise(i,j)-Pbar)^2;
      end
  end
  %��dbar
  tp1=0.0;
  for i=1:m
      for j=1:n
          tp1=tp1+da(i,j);
      end
  end
  dbar=tp1/(m*n);
  
  %��sigma
 tp1=0.0;
  for i=1:m
      for j=1:n
           tp1=tp1+(da(i,j)- dbar)^2;
      end
  end
  sigma=sqrt(tp1/(m*n-1));

while et>0.0001 && t<1000 %ѭ������
    v=v1;

  %��K
  for k=1:c
      for i=1:m
          for j=1:n
              
              K(i,j,k)=exp(-(data_noise(i,j)-v1(k))^2/sigma)+0.0001;
          end
      end
  end
  
 % ��G��

 for k=1:c
    for i=2:m-1
        for j=2:n-1
            tp1=0.0;
            for i1=-1:1
                for j1=-1:1
                   if i1~=0 || j1~=0
                    tp1=tp1+(1-u(i+i1,j+j1,k))^mc*(1-K(i+i1,j+j1,k))/((i1^2+j1^2)^0.5+1);
                   end
                end
            end
            G(i,j,k)=tp1;
       end
    end
 end
 
 % ��������
    for i=1:m
        for j=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+(1-K(i,j,k)+G(i,j,k)+0.0001)^(-1/(mc-1));
            end
            for k=1:c
                u(i,j,k)=(1-K(i,j,k)+G(i,j,k)+0.0001)^(-1/(mc-1))/tp1;
            end
        end
    end
 
  % ���¾�������
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for i=1:m
            for j=1:n
                tp1=tp1+(u(i,j,k)^mc)*K(i,j,k)*data_noise(i,j);
                tp2=tp2+(u(i,j,k)^mc)*K(i,j,k);
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
        et=0.0001;
   end
t=t+1;
end
 disp(['KFLICM����ʱ��',num2str(etime(clock,t1))]);
% [Iout]=qubian(data1,m,n);
% [m,n]=size(Iout);
% ����
% ����
I_KFLICM=zeros(m,n);
if c==2
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)
                I_KFLICM(i,j)=0;
            else
                I_KFLICM(i,j)=227;
            end
        end
    end
elseif c==3
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)
                I_KFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)
                I_KFLICM(i,j)=127;
            else
                I_KFLICM(i,j)=255;
            end
        end
    end
elseif c==4
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)
                I_KFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)
                I_KFLICM(i,j)=80;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)
                I_KFLICM(i,j)=150;
            else
                I_KFLICM(i,j)=255;
            end
        end
    end
elseif c==5
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)&&u(i,j,1)>u(i,j,5)
                I_KFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)&&u(i,j,2)>u(i,j,5)
                I_KFLICM(i,j)=100;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)&&u(i,j,3)>u(i,j,5)
                I_KFLICM(i,j)=155;
            elseif u(i,j,4)>u(i,j,1)&& u(i,j,4)>u(i,j,2)&&u(i,j,4)>u(i,j,3)&&u(i,j,4)>u(i,j,5)
                I_KFLICM(i,j)=200;
            else
                I_KFLICM(i,j)=255;
            end
        end
    end
end
I=I_KFLICM;
% [KFLICM_Vpc,KFLICM_SA]=huafenzhibiao3(u,m,n,I,c);
[KFLICM_Vpc,KFLICM_psnr,KFLICM_SA,KFLICM_Acc,KFLICM_Sen,KFLICM_Jaccard,KFLICM_Kappa]=huafenzhibiao(u,m,n,I,c);
