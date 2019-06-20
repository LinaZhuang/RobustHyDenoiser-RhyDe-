%%  The code and data herein distributed reproduce the results published in
%  the paper 
%
%  Lina Zhuang and Jose M. Bioucas-Dias, "Hyperspectral image denoising and 
%  anomaly detection based on low-rank and sparse representations" 
%  in Image and Signal Processing for Remote Sensing XXIII. Vol. 10427. 
%  International Society for Optics and Photonics, 2017.
%
%  
%  URL: http://www.lx.it.pt/~bioucas/files/spie_2017_rhyde.pdf
%
%%  Notes:
%
%   1) Package instalation: unzip the files to a directory and run the
%   scripts of "main_testRhyDe.m", which reproduces the 
%    results of RhyDe reported in Tab.1 of the above paper.
%
%
%   2) The script FastHyDe.m uses the package BM3D 
%      (v2 (30 January 2014))to implement the denoising algorithm BM3D 
%      introduced in 
%
%      K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image denoising by
%      sparse 3D transform-domain collaborative filtering," IEEE Trans. 
%      Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.
%
%      The BM3D package  is available at the 
%      BM3D web page:  http://www.cs.tut.fi/~foi/GCF-BM3D
%
%      Download this package and install it is the folder /BM3D
%
%
%   3) Function NAIRLMA_denosing.m used in main_testRhyDe.m as a competitor
%   is published in the following paper:
%
%      He, W., Zhang, H., Zhang, L., and Shen, H., “Hyperspectral image denoising 
%      via noise-adjusted iterative low-rank matrix approximation,”IEEE Journal 
%      of Selected Topics in Applied Earth Observations and RemoteSensing, 
%      3050–3061 (Jun. 2015).
%
%      NAIRLMA_denosing.m  is available at the author's homepage: 
%      https://sites.google.com/site/rshewei/home
%       
%      Download the code and install it in the folder /NAIRSVD_1
%
%
%   3) If you want to use RhyDe for your own data, go to folder
%   'RhyDe' and find functions RhyDe.m.
%      
%
%   
%% ACKNOWLEDGMENTS
%
% The authors acknowledge the following individuals and organisations:
%
%
%   - Prof. Paolo Gamba from Pavia university, 
%     for making the Pavia Centre data set available to the community.
%
%   - Authors of NAIRLMA method (Dr. Wei He, Prof. Hongyan Zhang, Prof. 
%     Liangpei Zhang, Prof. Huangfeng Shen) from Wuhan University, for 
%     making the NAIRLMA demo available to the community.
%
%   - Authors of BM3D method (Dr. Kostadin Dabov, Prof. Alessandro Foi, 
%     Prof. Vladimir Katkovnik, and Prof. Karen Egiazarian) from Tampere 
%     University of Technology, for making the BM3D demo software available 
%     to the community.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Lina Zhuang and Jose M. Bioucas Dias, May 2019
%

