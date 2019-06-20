% BRP

function [L_new] = BRP(D,r,power)

% input: D Noisy image
%        n the column of image(<= the number raw of image)
%        r the rank of iamge
%        power 计算次数，一般设置为0，可以提高精度
%   output: L_new the output image

%       需要m>n
if size(D,1)<size(D,2)
    error('the band should be smaller than the square of blocksize!');
end

L=D;
[m,n] = size(D);

  %%RSVD_BRP
    %Update of L
    Y2=randn(min(m,n),r);     % 改n为m
   for i=1:power+1
        Y1=L*Y2;
        Y2=L'*Y1;
   end
    [Q,R]=qr(Y2,0);
    L_new=(L*Q)*Q';
    
    
% %% RSVD
%    %Update of L
%     G=randn(min(m,n),r);     % 改n为m
%     H=L*G;
%     if power~=0;
%         for j=1:power
%           H = L*L'*H;
%         end
%     else
%     end
%     [Q,R] = qr(H,0);
%     L_new = Q*Q'*L;
        
     
%     L_new=Q*Q'*L;