 
clear all;
clc;
% addpath 'demo_HySime'
% %% simulated data
% 
% load 'Pavia_80'
% load 'noisePavia80';
% % 
% [m,n,p] = size(OriData3);
% %%%%%%-----experiment with different noise variance----%%%%%%%%%
% sigma = noisePavia80/255;
% %   sigma = 2*noiselevel/255;
% % sigma = noisetoy;
% oriData3_noise = zeros(m,n,p);
% for i =1:p
%     oriData3_noise(:,:,i) = OriData3(:,:,i)+sigma(i)*randn(m,n);
% end
% 
% %%%%%%-----experiment with constant noise variance----%%%%%%%%%
% % sigma = 10/255;
% % oriData3_noise = OriData3 + sigma *randn(m,n,p);
% 
% %%
% % parameter 
%  M=20;
%  stepsize=4;%8;
% %  lambda = 1;
% %%%%%%%%%%%%%%%%%%%%
% tic;
% [output_image,rsvdPSNR,rsvdSSIM] =NAIRLMA_denosing(oriData3_noise2,OriData3,M,stepsize,1);  %NAIRSVD
% % % [output_image] = NAIRLMA_denosing(oriData3_noise,OriData3,blocksize,stepsize,power);
% toc;



%% test indian pine
clear;clc;
addpath('C:\1Denoising\FastHyDe');
load Indian_pines;
oriData3_noise = indian_pines;

for i=1:220
    otmp = oriData3_noise(:,:,i);
    
    min_x(i) = min(otmp(:));
    max_x(i) = max(otmp(:));
    otmp = otmp - min_x(i);
    scale(i) = max_x(i)-min_x(i);
    %scale to [0,1]
    oriData3_noise2(:,:,i) = otmp/scale(i);
end

OriData3 = oriData3_noise2;
%%
% parameter 
 M=20;
 stepsize=4;%8;
%  lambda = 1;
%%%%%%%%%%%%%%%%%%%%
tic;
[output_image,rsvdPSNR,rsvdSSIM] =NAIRLMA_denosing(oriData3_noise2,OriData3,M,stepsize,1);  %NAIRSVD
% % [output_image] = NAIRLMA_denosing(oriData3_noise,OriData3,blocksize,stepsize,power);
toc;
figure;imagesc(output_image(:,:,1));


for i=1:220
   
    output_image2(:,:,i) = output_image(:,:,i)*scale(i)+min_x(i);
   
end
figure;imagesc(output_image2(:,:,1));