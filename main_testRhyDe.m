%
clear;clc;close all;
addpath('RhyDe');

load Simulatedtestdata\simulatedData_noisy_sigma04_Noniid_20p31;
 %load Simulatedtestdata\simulatedData_noisy_sigma02_Noniid_25p65; 
% load Simulatedtestdata\simulatedData_noisy_sigma013_Noniid_30p84; 
% load Simulatedtestdata\simulatedData_noisy_sigma007_Noniid_35p40; 
% load Simulatedtestdata\simulatedData_noisy_sigma004_Noniid_40p19; 


img_noisy = y_noisy; clear y_noisy;
[r, c, b] = size(img_noisy);
n  = r*c;
Y_noisy = reshape(img_noisy, r*c, b)';
Y_clean = y_clean; clear y_clean;
 img_clean = reshape(Y_clean', r,c,b);



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   RhyDe(Robust hyperspectral Denoising) 

 
pa=sqrt(chi2inv(0.99,b) ); % parameter of the regularizer promoting sparse matrix
pa = 10;
p =10;  %dimension of subspace 
%RhyDe is extremely robust to the subspace dimension overestimation. 
%For real images, without prior knowledge, you may set a larger value for p.
draw = 1;

t1=clock;
[EZ_est,S_est, Rw_save] =  RhyDe(img_noisy,pa,p,draw);
t2=clock;
time_rhyde = etime(t2,t1);
Y_rhyde = EZ_est + S_est;


result_rhyde = [ MSNR(Y_rhyde,Y_clean);MSSIM(Y_rhyde,Y_clean, r, c);time_rhyde];


img_rhyde = reshape(Y_rhyde',r,c,b);
S_est_img = reshape( sum(S_est.^2), r,c);


%% NAILRMA
addpath('NAIRSVD_1');
M= 20;
stepsize=8;
noise_type='additive';
t1 = clock;
[img_nailrma,rsvdPSNR,rsvdSSIM] =NAIRLMA_denosing(img_noisy,img_clean,M,stepsize,1);
time_nailrma = etime(clock,t1)
Y_nailrma = reshape(img_nailrma, [], b)';
result_NAIRLMA = [MSNR(Y_nailrma,Y_clean);MSSIM(Y_nailrma,Y_clean, r, c);time_nailrma];

%---------------------------------------------
figure(9); 
idx = 50;
 tmp =  img_clean(:,:,idx) ;
    tmp = tmp(:);
    tmp = sort(tmp);
    min_d = tmp(fix(n*0.02),1);
    max_d = tmp(fix(n*0.98),1);
subplot(1,4,1);
imshow(img_clean(:,:,idx),[min_d,max_d]),
title(['Clean band #',num2str(idx)]);
subplot(1,4,2);
imshow(img_noisy(:,:,idx),[min_d,max_d]),
title(['Noisy band #',num2str(idx)]);
subplot(1,4,3);
imshow(img_nailrma(:,:,idx),[min_d,max_d]),
title({['NAIRLMA band #',num2str(idx)]; ['MSNR=',num2str(result_NAIRLMA(1)),...
    'dB,  MSSIM=',num2str(result_NAIRLMA(2))]});
subplot(1,4,4);
imshow(img_rhyde(:,:,idx),[min_d,max_d]),
 title({['RhyDe band #',num2str(idx)]; ['MSNR=',num2str(result_rhyde(1)),...
    'dB,  MSSIM=',num2str(result_rhyde(2))]}); 


figure(99);
subplot(1,2,1);
 
loc = [  23,   49208  ];
name_title = {'(a) a normal pixel', '(b) a anomaly pixel'};
for i=1:2
    subplot(1,2,i );
    plot( [Y_clean(:,loc(i)),Y_noisy(:,loc(i)),Y_nailrma(:,loc(i)),Y_rhyde(:,loc(i))]);
    legend({'Clean', 'Noisy','NAILRMA','RhyDe'});
    title(name_title{i});
end
% Note that the noise in anomaly pixel is not removed completelyby RhyDe 
% since our main objective w.r.t.  anomalies is to keep them rather than 
% to denoise them and our output result is \hat{Z}+\hat{S}.
 
