function [EZ_est,S_est,Rw] =  RhyDe(img_noisy,pa,p,draw)

%Input:
% img_noisy          noisy 3D data of size row*column*band
% pa                 pa is corresponding to regularization parameter lambda_2 in
%                    our objective funcion, which need to be tuned carefully.
%                    pa=sqrt(chi2inv(0.99,band) ), when noise is Gaussian N(0,1) with std 1.
% p                  subsapce dimension 
% draw               draw=1, meaning display results, otherwise, working silently.

% Output:
% EZ_est             Denoised 3D image
% S_est              3D image with outlier component
% ---------------------------- -------------------------------------------
% Given model Y=EZ+S+N, this script estimates the matrix Z and the sparse matrix S representing
% outliers by solving the optimization:
% {\hat{Z}, \hat{S}} \in \min_{Z,S}  1/2 || EZ+S-Y||_F^2 + \lambda_1
% \phi(Z) + \lambda_2 ||S||_{2,1},
% 
% See more details in papers:
%   [1] L. Zhuang and J. M. Bioucas-Dias, 
%       "Hyperspectral image denoising and anomaly detection based on low-rank 
%        and sparse representations" in Image and Signal Processing for Remote 
%        Sensing XXIII. Vol. 10427. International Society for Optics and Photonics, 2017.
%
%        URL: http://www.lx.it.pt/~bioucas/files/spie_2017_rhyde.pdf
%
%% -------------------------------------------------------------------------
%
% Copyright (May, 2019):        
%             Lina Zhuang (lina.zhuang@lx.it.pt)
%             &
%             Jose Bioucas-Dias (bioucas@lx.it.pt)
%            
%
% RhyDe is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------




 [row, column, band] = size(img_noisy);
N=row*column;
Y_noisy = reshape(img_noisy, [N band])';


[w Rw_correlation] = estNoise(Y_noisy, 'additive');
Rw_correlation=Rw_correlation;%  +eye(band)*0.000001;
%whiten data: Each band is divided by its noise standard deviation, leading to 
%standard Gaussian noise N(0,1).
Y_noisy = sqrt(inv(Rw_correlation))*Y_noisy; 
 



%% subspace estimation: (using HySime or SVD)
%SVD
[E,D]=svd(Y_noisy*Y_noisy');
E(:,p+1:end) = [];
Z = E'*Y_noisy;

 

%% initializtion:
k = 0;
I_nb = eye(band);
I_p = eye(p);
zero_p_nb = zeros(p, band);
mu1 =  1;
mu2 =  1;
mu3 =  1;
S = zeros(size(Y_noisy));
A = [Z;S];
V1 = [E, I_nb]*A;
V2 = [I_p, zero_p_nb]*A;
V3 = [zero_p_nb',I_nb]*A;
D1 = zeros(size(V1));
D2 = zeros(size(V2));
D3 = zeros(size(V3));
Y = Y_noisy;



lambda2 =pa*mu3;
%% iteration
noise_std_ini = [];
while 1
    
    k = k+1;
    %% update A
    A_tmp = mu1*[E,I_nb]'*[E,I_nb] + mu2*[I_p, zero_p_nb]'*[I_p, zero_p_nb] + mu3*[zero_p_nb',I_nb]'*[zero_p_nb',I_nb];
    A = inv(A_tmp)* (mu1*[E,I_nb]'*(V1-D1)+mu2*[I_p, zero_p_nb]'*(V2-D2)+mu3*[zero_p_nb',I_nb]'*(V3-D3) );
    
    
    %% update V1
    V1 = 1/(1+mu1)*(Y+mu1*([E,I_nb]*A+D1));
    
    %% update V2
    V2_tmp = [I_p, zero_p_nb]*A+D2;
    V2_tmp_img = reshape(V2_tmp', row, column ,[]);
    
    
    %%%%%%%%%%%%% estimate noise:
    if k==1
        V2_tmp_img_init=V2_tmp_img;     
            noise_std_ini = ones(p,1); %noise std is initialized as 1, because the image has been whiten.
        noise_std = noise_std_ini;
        
        
    else
        
        w=V2_tmp_img-V2_tmp_img_init;
        w = reshape(w,N,p)';
         Rw=diag(diag(w*w'/N));
        for cp=1:p
            noise_std_tmp(cp,1) = sqrt( Rw(cp,cp));
            if noise_std_ini(cp,1) > noise_std_tmp(cp,1)
                noise_std(cp,1)= noise_std_ini(cp,1)-noise_std_tmp(cp,1);
            else
                noise_std(cp,1) = 0;
            end
        end
        
        
        
    end
 
    if sum(noise_std)==0
        break;
    end
    
    
    
    M=ones(size(V2_tmp_img));
    
    
    
    [V2_fasthyde] = FastHyDe_outliers(V2_tmp_img,noise_std);
    V2 = reshape(V2_fasthyde, N, [])';
    
    %% update V3
    V3_tmp = [zero_p_nb',I_nb]*A+D3;
    V3 = L21norm(V3_tmp, lambda2/mu3);
    sqrt(sum((V3_tmp.^2)));
    if draw
       figure(111); histogram(ans); title({['Histogram of elements in sparse matrix S'];['lambda2= ',num2str(lambda2/mu3), ',   ite=',num2str(k)]});
    end
    %% update D
    D1 = D1 - (V1-[E,I_nb]*A);
    D2 = D2 - (V2-[I_p,zero_p_nb]*A);
    D3 = D3 - (V3-[zero_p_nb',I_nb]*A);
    
    Y_rec = E*A(1:p,:);
 
    if k==10
        break;
    end
    
 
end

EZ_est = E*A(1:p,:);
S_est = A(p+1:end,:);


  EZ_est = sqrt(Rw_correlation)*EZ_est;
  S_est = sqrt(Rw_correlation)*S_est;



