function [image_fasthyde] = FastHyDe_outliers(image_ori,SIGMA)

[Lines, Columns, B] = size(image_ori);
N=Lines*Columns;
for i=1:size(image_ori,3)
    Y(i,1:N)= reshape(image_ori(:,:,i),[1,N]);
end

%% denoising:



% [U_ss,D]=svd(Y*Y');
% U_ss(:,5:B) = [];
% Y_pro = U_ss*U_ss'*Y;
% 
% Rw=diag(diag((Y_pro-Y)*(Y_pro-Y)'/N));
% sqrt(diag(Rw))



y_est_bm3d = zeros(Lines,Columns,B);

eigen_Y_bm3d=[];
 
for i=1:B
    % produce eigen-image
    eigen_im = Y(i,:);
    min_x = min(eigen_im);
    max_x = max(eigen_im);
    eigen_im = eigen_im - min_x;
    scale = max_x-min_x;
    
    %scale to [0,1]
    eigen_im = reshape(eigen_im, Lines, Columns)/scale;
    %     figure(1);
    %     imagesc(eigen_im);
    
    %estimate noise from Rw
   
      sigma =  SIGMA(i)/scale;
      
     filt_eigen_im = eigen_im;
%     if (i==1||2)
%              filt_eigen_im= medfilt2(eigen_im, [3 3]); %if [3 3], we lost some details 
%     else %if i>2
%         % denoise  with BM3D
%         [dummy, filt_eigen_im] = BM3D(1,eigen_im, sigma*255);
%     end
    if sigma>0
        addpath('BM3D');
 [dummy, filt_eigen_im] = BM3D(1,eigen_im, sigma*255);

    end
 


    eigen_Y_bm3d(i,:) = reshape(filt_eigen_im*scale + min_x, 1,N);
    
    if 0
      hfig1 =figure(111);
       set(hfig1, 'unit', 'normalized', 'position', [0,0,1,1]);
     subplot(1,3,1);
     imagesc(eigen_im);
     title(['band=',num2str(i),'Noisy']);
     subplot(1,3,2);
     imagesc(filt_eigen_im);
     title('Denoised');
     if i>1
     subplot(1,3,3);
%      imagesc(abs(filt_eigen_im - eigen_im));
     scatter(eigen_im(:),eigen_im_pre(:)); hold on;
    scatter(filt_eigen_im(:), filt_eigen_im_pre(:),'.'); hold off;
%      xlabel('Error');
     end
      pause(1);
    end
    
      eigen_im_pre = eigen_im;
      filt_eigen_im_pre =filt_eigen_im;
end

% reconstruct data using denoising engin images
Y_reconst_bm3d = eigen_Y_bm3d;


image_fasthyde=[];
for i=1:B
    image_fasthyde(1:Lines,1:Columns,i) = reshape(Y_reconst_bm3d(i,:),[Lines,Columns]);
end

t2=clock;

end