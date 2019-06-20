function k = MSNR(Y, Y_ref)
% input size : Bands*observation

[B,n] = size(Y);
Err = Y-Y_ref;
for i=1:size(Y,1)
k_tmp(i)  = snr(Y_ref(i,:),Y(i,:)-Y_ref(i,:));
end
k=mean(k_tmp);
fprintf('\n The Mean of SNR value is %0.2f', k);