function [Vpc,psnr,SA,Acc,Sen,Jaccard,Kappa]=huafenzhibiao(u,m,n,I,c)
% I1=imread('E:\小论文\论文3\对比\理想\医学图\eFCM.png');
% I1=imread('E:\小论文\论文3\对比\理想\iFCM.png');
I1=imread('E:\小论文\论文3\对比\椒盐\a原图.png');
I=double(I);
I1=double(I1);

%Vpc的计算
tp1=0.0;
for i=1:m
    for j=1:n
        for k=1:c
            tp1=tp1+u(i,j,k)^2;
        end
    end
end
Vpc=tp1/(m*n);

%Vpe的计算
% tp1=0.0;
% for i=1:m
%     for j=1:n
%         for k=1:c
%             tp1=tp1+u(i,j,k)*log(u(i,j,k));
%         end
%     end
% end
% Vpe=tp1/(-m*n);

%Vxb的计算
% tp1=0.0;
% for i=1:m
%     for j=1:n
%         for k=1:c
%             tp1=tp1+u(i,j,k)^2*(data(i,j)-v1(k))^2;
%         end
%     end
% end
% tp2=tp1/(m*n);
% i=1;
% for k=1:c
%     for k1=1:c
%         if k1~=k
%             V(i)=(v1(k1)-v1(k))^2;
%             i=i+1;
%         end
%     end
% end
%             
% Vxb=tp2/min(V);

%SC的计算
% for k=1:c
%     tp1=0.0;
%     for i=1:m
%         for j=1:n
%             tp1=tp1+u(i,j,k)^2*(data(i,j)-v1(k))^2;
%         end
%     end
%     S(k)=tp1;
% end
% 
% tp2=0.0;
% for k1=1:c
%     for k2=1:c
%         if k1~=k2
%             tp2=tp2+(v1(k1)-v1(k2))^2;
%         end
%     end
% end
% tp3=0.0;
% for k=1:c
%     tp3=tp3+S(k)/(tp2*m*n);
% end
% SC=tp3;

%峰值信噪比(psnr)的计算
tp1=0;
for i=1:m
    for j=1:n
        tp1=tp1+(I1(i,j)-I(i,j))^2;
    end
end
mse=tp1/(m*n);
psnr=10*log10(double(255^2/mse));

%SA的计算
geshu=0;
for i=1:m
    for j=1:n
        if abs(I(i,j)-I1(i,j))==0
            geshu=geshu+1;
        end
    end
end
SA=geshu/(m*n)*100;

%ME的计算
% geshu=0;
% for i=1:m
%     for j=1:n
%         if abs(I(i,j)-I1(i,j))==0
%             geshu=geshu+1;
%         end
%     end
% end
% ME=(1-geshu/(m*n))*100;

% Acc、Sen、Spe的计算
% 分割指标
Lab=zeros(m,n,c);
Lab1=zeros(m,n,c);
label=zeros(m,n);
Label=zeros(m,n);
for i=1:m   %实际标准图像的标签
    for j=1:n
           if c==5
            if I1(i,j)==0
                label(i,j)=1;
            end
            if I1(i,j)==65
                label(i,j)=2;
            end
            if I1(i,j)==130
                label(i,j)=3;
            end
             if I1(i,j)==200
                label(i,j)=4;
             end
              if I1(i,j)==255
                label(i,j)=5;
              end
          end
        if c==4
            if I1(i,j)==0
                label(i,j)=1;
            end
            if I1(i,j)==50
                label(i,j)=2;
            end
            if I1(i,j)==150
                label(i,j)=3;
            end
             if I1(i,j)==255
                label(i,j)=4;
             end
        end
        
        if c==2
            if I1(i,j)==0
                label(i,j)=1;
            end
            if I1(i,j)==255
                label(i,j)=2;
            end
        end
        if c==3
            if I1(i,j)==0
                label(i,j)=1;
            end
            if I1(i,j)==150
                label(i,j)=2;
            end
            if I1(i,j)==255
                label(i,j)=3;
            end
        end
        for k=1:c
            if  label(i,j)==k
                Lab(i,j,k)=1.0;
            else
                Lab(i,j,k)=0;
            end
        end
    end
end
 %分割后所得的估计标签
for i=1:m 
    for j=1:n
        if c==4
            if I(i,j)==0
                Label(i,j)=1;
            end
            if I(i,j)==50
                Label(i,j)=2;
            end
            if I(i,j)==150
                Label(i,j)=3;
            end
             if I(i,j)==255
                Label(i,j)=4;
             end
        end
        
        if c==2
            if I(i,j)==0
                Label(i,j)=1;
            end
            if I(i,j)==255
                Label(i,j)=2;
            end
        end
        if c==3
            if I(i,j)==0
                Label(i,j)=1;
            end
            if I(i,j)==150
                Label(i,j)=2;
            end
            if I(i,j)==255
                Label(i,j)=3;
            end
        end
        for k=1:c
            if  Label(i,j)==k
                Lab1(i,j,k)=1.0;
            else
                Lab1(i,j,k)=0;
            end
        end
    end
end

% 计算TP，FP，TN，FN
tp11=0.0;
tp12=0.0;
tp13=0.0;
tp14=0.0;
for i=1:m
    for j=1:n
        for k=1:c
            tp11=tp11+Lab1(i,j,k)*Lab(i,j,k);
            tp12=tp12+Lab1(i,j,k)*(1-Lab(i,j,k));
            tp13=tp13+(1-Lab1(i,j,k))*(1-Lab(i,j,k));
            tp14=tp14+(1-Lab1(i,j,k))*Lab(i,j,k);
        end
    end
end
TP=tp11;
FP=tp12;
TN=tp13;
FN=tp14;
% 计算准确度，灵敏度，特异度
Acc=(TP+TN)/(TP+TN+FP+FN);
Sen=TP/(TP+FN);
Jaccard=TP/(TP+FP+FN);
% Jaccard = double(sum(uint8(image2(:) & data0(:)))) / double(sum(uint8(image2(:) | data0(:))));

L1=find(Label==1);N1=size(L1,1);
L2=find(Label==2);N2=size(L2,1);
L3=find(Label==3);N3=size(L3,1);
L4=find(Label==4);N4=size(L4,1);

l1=find(label==1);n1=size(l1,1);
l2=find(label==2);n2=size(l2,1);
l3=find(label==3);n3=size(l3,1);
l4=find(label==4);n4=size(l4,1);

pe=(N1*n1+N2*n2+N3*n3+N4*n4)/(n*m*n*m);

i=find(Label==label);
po=size(i,1)/(n*m); %等同于SA
Kappa=(abs(po-pe))/(1-pe);

% Acc、Sen、Spe的计算
% 分割指标
% label=zeros(m,n);
% Lab=zeros(m,n,c);
% Lab1=zeros(m,n,c);
% if c==2
%     for i=1:m   %实际标准图像的标签
%         for j=1:n
%             if I1(i,j)==0
%                 label(i,j)=2;
%             else
%                 label(i,j)=1;
%             end
%             for k=1:c
%                 if abs(label(i,j)-k)<0.01
%                     Lab(i,j,k)=1.0;
%                 else
%                     Lab(i,j,k)=0;
%                 end
%             end
%         end
%     end
% elseif c==3
%     for i=1:m   %实际标准图像的标签
%         for j=1:n
%             if I1(i,j)==0
%                 label(i,j)=1;
%             elseif   I1(i,j)==128
%                 label(i,j)=2;
%             else
%                 label(i,j)=3;
%             end
%             for k=1:c
%                 if  abs(label(i,j)-k)<0.01
%                     Lab(i,j,k)=1.0;
%                 else
%                     Lab(i,j,k)=0;
%                 end
%             end
%         end
%     end
% elseif c==4
%     for i=1:m   %实际标准图像的标签
%         for j=1:n
%             if I1(i,j)==0;
%                 label(i,j)=1;
%             elseif  I1(i,j)==100
%                 label(i,j)=2;
%             elseif I1(i,j)==155
%                 label(i,j)=3;
%             else
%                 label(i,j)=4;
%             end
%             for k=1:c
%                 if  abs(label(i,j)-k)<0.01
%                     Lab(i,j,k)=1.0;
%                 else
%                     Lab(i,j,k)=0;
%                 end
%             end
%         end
%     end
% elseif c==5
%     for i=1:m   %实际标准图像的标签
%         for j=1:n
%             if I1(i,j)==0
%                 label(i,j)=1;
%             elseif I1(i,j)==100
%                 label(i,j)=2;
%             elseif I1(i,j)==155
%                 label(i,j)=3;
%             elseif I1(i,j)==200
%                 label(i,j)=4;
%             else
%                 label(i,j)=5;
%             end
%             for k=1:c
%                 if  abs(label(i,j)-k)<0.01
%                     Lab(i,j,k)=1.0;
%                 else
%                     Lab(i,j,k)=0;
%                 end
%             end
%         end
%     end
% end
% 
% for i=1:m  %分割后所得的估计标签
%     for j=1:n
%         for k=1:c
%             tp=u(i,j,k);
%             t1=k;
%             for k1=1:c
%                 if tp<u(i,j,k1)
%                     tp=u(i,j,k1);
%                     t1=k1;
%                 end
%             end
%             if k==t1
%                 Lab1(i,j,k)=1.0;
%             else
%                 Lab1(i,j,k)=0;
%             end
%         end
%     end
% end
% 
% %计算TP，FP，TN，FN
% tp11=0.0;
% tp12=0.0;
% tp13=0.0;
% tp14=0.0;
% for i=1:m
%     for j=1:n
%         for k=1:c
%             tp11=tp11+Lab1(i,j,k)*Lab(i,j,k);
%             tp12=tp12+Lab1(i,j,k)*(1-Lab(i,j,k));
%             tp13=tp13+(1-Lab1(i,j,k))*(1-Lab(i,j,k));
%             tp14=tp14+(1-Lab1(i,j,k))*Lab(i,j,k);
%         end
%     end
% end
% TP=tp11;
% FP=tp12;
% TN=tp13;
% FN=tp14;
% %计算准确度，灵敏度，特异度
% Acc=(TP+TN)/(TP+TN+FP+FN)*100;
% Sen=TP/(TP+FN)*100;
% % Spe=TN/(TN+FP)*100;
% Jaccard=TP/(TP+FP+FN);
% 
% 
% L1=find(Lab1==1);N1=size(L1,1);
% L2=find(Labl==2);N2=size(L2,1);
% L3=find(Labl==3);N3=size(L3,1);
% L4=find(Labl==4);N4=size(L4,1);
% 
% l1=find(labl==1);n1=size(l1,1);
% l2=find(labl==2);n2=size(l2,1);
% l3=find(labl==3);n3=size(l3,1);
% l4=find(labl==4);n4=size(l4,1);
% 
% pe=(N1*n1+N2*n2+N3*n3+N4*n4)/(n*m*n*m);
% i=find(Labl==label);
% po=size(i,1)/(n*m); %等同于SA
% Kappa=(abs(po-pe))/(1-pe);


