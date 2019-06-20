%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demo file: Signal subspace estimation with
%            HySime method
%
addpath 'D:\»•‘Î\SVT∑÷Ω‚'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, clc, close all
verbose='on';

p = 15;              % number of endmembers
N = 10^4;           % number of pixels
SNR = 25;           % in dB

% eta=0;              % white noise
%or
eta = 1/18;         % color_noise;

%noise_type = 'poisson';
%or
noise_type = 'additive';

[y,x,M,s,n,Cn,wl]=generate_hyperspectral_data(p,N,SNR,eta,noise_type,verbose);
L=length(wl);  % number of bands

% % denoising
% g = 20;
% power = 1;
% 
% y = y';
% ystand = y;
% 
% L_new = BRP(y,g,power);
% for i =1:10
%    y = L_new + 0.5*(ystand-L_new); 
%    L_new = BRP(y,g,power);
% end
% y = L_new';
% % denoising

[w Rn] = estNoise(y,noise_type,verbose);

% sel_pixel = round((N-1)*rand(1));
% figure(1);
%   set(gca,'FontSize',12,'FontName','times new roman');
%   subplot(221);
%      plot(wl,x(:,sel_pixel),'r--',wl,y(:,sel_pixel),'b-','Linewidth',1);
%      axis([min(wl) max(wl) 0 1.2*max(y(:,sel_pixel))])
%      xlabel('\lambda(\mum)')  
%      title('pixel signature');
%      legend('x','x+n')
%   subplot(222);
%      plot(wl,n(:,sel_pixel),'r--',wl,w(:,sel_pixel),'b-','Linewidth',1);
%      axis([min(wl) max(wl) min([n(:,sel_pixel); w(:,sel_pixel)]) max([n(:,sel_pixel); w(:,sel_pixel)])])
%      xlabel('\lambda(\mum)')  
%      title('noise and its estimates');
%      legend('n','w')
%   subplot(223);
%      plot(wl,x(:,sel_pixel),'r--',wl,y(:,sel_pixel)-w(:,sel_pixel),'b-','Linewidth',1);
%      axis([min(wl) max(wl) 0 1.2*max(y(:,sel_pixel))])
%      xlabel('\lambda(\mum)')  
%      title('pixel signature and its estimates');
%      legend('x','y-w')
%   subplot(224);
%      plot(wl,diag(Cn),'r--',wl,diag(Rn),'b-','Linewidth',1);
%      axis([min(wl) max(wl) 0 1.2*max([Cn(:); Rn(:)])])
%      xlabel('\lambda(\mum)')  
%      title('diagonal of the noise covariance matrix');
%      legend('true noise','noise estimates')

[kf Ek]=hysime(y,w,Rn,verbose);    % estimate the p
if ~strcmp(verbose,'on'),fprintf(1,'The signal subspace dimension is: k = %d\n',kf);end

% sigma = mean(mean(w));
% [Image,iter ] = SVT_soft(y,sigma);
% 
% 
% 
%  [U,Sig,V] = svd(y-w,'econ');
% nnn = zeros(1,224);
% for i=1:224
% nnn(1,i) = U(:,i)'*w*V(:,i);
% end
% 
%  sig = diag(Sig);
% nn = zeros(1,224);
% for i=1:224
% nn(1,i) = sig(i)/nnn(i);
% end
%  return
