clc;clear all;close all;
data=imread('06.png');
data=rgb2gray(data);
figure(1);imshow(uint8(data)); 
imwrite(data,'E:\小论文\论文3\对比\椒盐\a原图.png');
title('原图');

% data=padarray(data,[2,2]);
data=padarray(data,[2,2],'replicate');

% data_noise=imnoise(data,'gaussian',0,0.4);
% figure(2);imshow(data_noise);
% imwrite(data_noise,'E:\小论文\论文3\对比\标准数据集\a噪声图_0.4.png');
% title('高斯噪声图');
 
data_noise=imnoise(data,'salt & pepper',0.5);
figure(2);imshow(data_noise);
imwrite(data_noise,'E:\小论文\论文3\对比\椒盐\a噪声图_0.5.png');
title('椒盐噪声图');

% data_noise=imnoise(data,'speckle',0.1);
% figure(2);imshow(data_noise);
% imwrite(data_noise,'E:\小论文\论文3\对比\标准数据集\d噪声图_0.1.png');
% title('噪声图');

% data=double(data);
% ss=70;
% noise_real=ss*randn(size(data))+data;
% noise_image=ss*randn(size(data))+data;
% data_noise=round(sqrt(noise_real.^2+noise_image.^2));
% data_noise=double(data_noise);
% figure(2);
% imshow(uint8(data_noise));
% f=getimage(figure(2)); 
% imwrite(f,'E:\小论文\论文3\对比\莱斯\b噪声图_70.png');
% title('莱斯噪声图像')

data_noise=double(data_noise);
[m,n]=size(data_noise);
mc=2;
e=0.1;
ct=0;
c=3;
if c==2
    v1(1)=0;v1(2)=255; 
elseif c==3
    v1(1)=0;v1(2)=125;v1(3)=255; 
elseif c==4
    v1(1)=0;v1(2)=80;v1(3)=150;v1(4)=255;
end
 
[I_KFCM_S,KFCM_S_Vpc,KFCM_S_psnr,KFCM_S_SA,KFCM_S_Acc,KFCM_S_Sen,KFCM_S_Jaccard,KFCM_S_Kappa]=KFCM_S(data_noise,v1,m,n,c,mc,e,ct);
% [I_KFCM_S1,KFCM_S1_Vpc,KFCM_S1_psnr,KFCM_S1_SA,KFCM_S1_Acc,KFCM_S1_Sen,KFCM_S1_Jaccard,KFCM_S1_Kappa]=KFCM_S1(data_noise,v1,m,n,c,mc,e,ct);
% [I_KFCM_S2,KFCM_S2_Vpc,KFCM_S2_psnr,KFCM_S2_SA,KFCM_S2_Acc,KFCM_S2_Sen,KFCM_S2_Jaccard,KFCM_S2_Kappa]=KFCM_S2(data_noise,v1,m,n,c,mc,e,ct);
[I_FLICM,FLICM_Vpc,FLICM_psnr,FLICM_SA,FLICM_Acc,FLICM_Sen,FLICM_Jaccard,FLICM_Kappa]=FLICM(data_noise,v1,m,n,c,mc,e,ct);
[I_RFLICM,RFLICM_Vpc,RFLICM_psnr,RFLICM_SA,RFLICM_Acc,RFLICM_Sen,RFLICM_Jaccard,RFLICM_Kappa]=RFLICM(data_noise,v1,m,n,c,mc,e,ct);
[I_KFLICM,KFLICM_Vpc,KFLICM_psnr,KFLICM_SA,KFLICM_Acc,KFLICM_Sen,KFLICM_Jaccard,KFLICM_Kappa]=KFLICM(data_noise,v1,m,n,c,mc,e,ct);
[I_WFLICM,WFLICM_Vpc,WFLICM_psnr,WFLICM_SA,WFLICM_Acc,WFLICM_Sen,WFLICM_Jaccard,WFLICM_Kappa]=WFLICM(data_noise,v1,m,n,c,mc,e,ct);
[I_KWFLICM,KWFLICM_Vpc,KWFLICM_psnr,KWFLICM_SA,KWFLICM_Acc,KWFLICM_Sen,KWFLICM_Jaccard,KWFLICM_Kappa]=KWFLICM(data_noise,v1,m,n,c,mc,e,ct);
% [I_PFCM,PFCM_Vpc,PFCM_psnr,PFCM_SA,PFCM_Acc,PFCM_Sen,PFCM_Jaccard,PFCM_Kappa]=PFCM(data_noise,v1,m,n,c,e,ct);
[I_PFLICM,PFLICM_Vpc,PFLICM_psnr,PFLICM_SA,PFLICM_Acc,PFLICM_Sen,PFLICM_Jaccard,PFLICM_Kappa]=PFLICM(data_noise,v1,m,n,c,e,ct);
[I_p_w,p_w_Vpc,p_w_psnr,p_w_SA,p_w_Acc,p_w_Sen,p_w_Jaccard,p_w_Kappa]=p_w(data_noise,v1,m,n,c,e,ct);
% [I_p_w1,p_w1_Vpc,p_w1_psnr,p_w1_SA,p_w1_Acc,p_w1_Sen,p_w1_Jaccard,p_w1_Kappa]=p_w1(data_noise,v1,m,n,c,e,ct);
% [I_w,w_Vpc,w_psnr,w_SA,w_Acc,w_Sen,w_Jaccard,w_Kappa]=w(data_noise,v1,m,n,c,mc,e,ct);

figure(3)
imshow(uint8(I_KFCM_S));
f=getimage(figure(3));
imwrite(f,'E:\小论文\论文3\对比\椒盐\aKFCM_S_0.5.png');
title('KFCM_S')
fprintf('KFCM_S:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',KFCM_S_Vpc,KFCM_S_psnr,KFCM_S_SA,KFCM_S_Acc,KFCM_S_Sen,KFCM_S_Jaccard,KFCM_S_Kappa);

% figure(4);
% imshow(uint8(I_KFCM_S1));
% f=getimage(figure(4));
% imwrite(f,'E:\小论文\论文3\对比\莱斯\dKFCM_S1_0.5.png');
% title('KFCM_S1')
% fprintf('KFCM_S1:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',KFCM_S1_Vpc,KFCM_S1_psnr,KFCM_S1_SA,KFCM_S1_Acc,KFCM_S1_Jaccard,KFCM_S1_Sen,KFCM_S1_Kappa);
% 
% figure(5)
% imshow(uint8(I_KFCM_S2));
% f=getimage(figure(5));
% imwrite(f,'E:\小论文\论文3\对比\莱斯\dKFCM_S2_0.5.png');
% title('KFCM_S2')
% fprintf('KFCM_S2:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',KFCM_S2_Vpc,KFCM_S2_psnr,KFCM_S2_SA,KFCM_S2_Acc,KFCM_S2_Sen,KFCM_S2_Jaccard,KFCM_S2_Kappa);

figure(6)
imshow(uint8(I_FLICM));
f=getimage(figure(6));
imwrite(f,'E:\小论文\论文3\对比\椒盐\aFLICM_0.5.png');
title('FLICM')
fprintf('FLICM:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',FLICM_Vpc,FLICM_psnr,FLICM_SA,FLICM_Acc,FLICM_Sen,FLICM_Jaccard,FLICM_Kappa);

figure(7)
imshow(uint8(I_RFLICM));
f=getimage(figure(7));
imwrite(f,'E:\小论文\论文3\对比\椒盐\bRFLICM_0.5.png');
title('RFLICM')
fprintf('RFLICM:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',RFLICM_Vpc,RFLICM_psnr,RFLICM_SA,RFLICM_Acc,RFLICM_Sen,RFLICM_Jaccard,RFLICM_Kappa);

figure(8)
imshow(uint8(I_KFLICM));  
f=getimage(figure(8));
imwrite(f,'E:\小论文\论文3\对比\椒盐\aKFLICM_0.5.png');
title('KFLICM')
fprintf('KFLICM:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',KFLICM_Vpc,KFLICM_psnr,KFLICM_SA,KFLICM_Acc,KFLICM_Sen,KFLICM_Jaccard,KFLICM_Kappa);

figure(9)
imshow(uint8(I_WFLICM));
f=getimage(figure(9));
imwrite(f,'E:\小论文\论文3\对比\椒盐\aWFLICM_0.5.png');
title('WFLICM')
fprintf('WFLICM:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',WFLICM_Vpc,WFLICM_psnr,WFLICM_SA,WFLICM_Acc,WFLICM_Sen,WFLICM_Jaccard,WFLICM_Kappa);

figure(10)
imshow(uint8(I_KWFLICM));
f=getimage(figure(10));
imwrite(f,'E:\小论文\论文3\对比\椒盐\aKWFLICM_0.5.png');
title('KWFLICM')
fprintf('KWFLICM:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',KWFLICM_Vpc,KWFLICM_psnr,KWFLICM_SA,KWFLICM_Acc,KWFLICM_Sen,KWFLICM_Jaccard,KWFLICM_Kappa);

% figure(11)
% imshow(uint8(I_PFCM));
% f=getimage(figure(11));
% imwrite(f,'E:\小论文\论文3\对比\椒盐\aPFCM_0.5.png');
% title('PFCM')
% fprintf('PFCM:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',PFCM_Vpc,PFCM_psnr,PFCM_SA,PFCM_Acc,PFCM_Sen,PFCM_Jaccard,PFCM_Kappa);

figure(12)
imshow(uint8(I_PFLICM));
f=getimage(figure(12));
imwrite(f,'E:\小论文\论文3\对比\椒盐\aPFLICM_0.5.png');
title('PFLICM')
fprintf('PFLICM:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',PFLICM_Vpc,PFLICM_psnr,PFLICM_SA,PFLICM_Acc,PFLICM_Sen,PFLICM_Jaccard,PFLICM_Kappa);

figure(13)
imshow(uint8(I_p_w));
f=getimage(figure(13));
imwrite(f,'E:\小论文\论文3\对比\椒盐\ap_w_0.5.png');
title('p_w')
fprintf('p_w:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',p_w_Vpc,p_w_psnr,p_w_SA,p_w_Acc,p_w_Sen,p_w_Jaccard,p_w_Kappa);

% figure(14)
% imshow(uint8(I_p_w1));
% f=getimage(figure(14));
% imwrite(f,'E:\小论文\论文3\对比\高斯\fp_w1_0.05.png');
% title('p_w1')
% fprintf('p_w1:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',p_w1_Vpc,p_w1_psnr,p_w1_SA,p_w1_Acc,p_w1_Sen,p_w1_Jaccard,p_w1_Kappa);
% 
% figure(15);
% imshow(uint8(I_w)); 
% f=getimage(figure(15));
% imwrite(f,'E:\小论文\论文3\对比\高斯\fw_0.1.png');
% title('w')
% fprintf('w:\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n%.4f\n',w_Vpc,w_psnr,w_SA,w_Acc,w_Sen,w_Jaccard,w_Kappa);
