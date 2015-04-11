
%
%% script to illustrate the use, and evaluate the efficacy of the Robust FCM implementation
%     also performs a comparison of FCM with the robust counterpart.
%


clear ; clc;

% %=================================
%
% toy example 1 (uncomment to use it)
% TestImg = zeros(100,100);
% TestImg(10:59,10:59) = 25;
% TestImg(70:90,60:90) = 25;
% TestImg(60:90,10:59) = 50;
% TestImg(60:70,60:90) = 50;
% TestImg(10:59,60:90) = 75;

% % toy example 2 (uncomment to use it)
TestImg = zeros(60,60);
TestImg( 5:25, 5:39) = 25;
TestImg(26:39, 5:39) = 50;
TestImg( 5:39,40:55) = 75;
TestImg(40:55, 5:32) = 50;
TestImg(40:55,33:55) = 25;


% adding some gaussian noise
NoiseVar = 9 ;
TestImg = TestImg + NoiseVar*randn(size(TestImg));    
TypeOfNbrhood = '4nbr';
TypeOfExpt = 'ToyImgDzungPham';
NumClusters = 4;
ImgDim = size(TestImg);
figure(gcf); clf;
imagesc(TestImg); 
data = TestImg(:);

% rfcm parameters
% % good params q,beta = 5,0.5;
% % excellent params q,beta,noise_var = 5,0.5,2;
opt.ExpntQ = 1.2;
opt.beta = 300 ;

opt.MaxIter = 20;
opt.tol = 0.01;

opt.TypeOfNbrhood = '4nbr';
opt.visualize = false;
opt.verbose = true;

%% do the FCM clustering
% [ centers_fcm, U_fcm_est, ObjFun_fcm ] = robustFCM(data, ImgDim, NumClusters, opt, TypeOfExpt);
[ centers_fcm, U_fcm_est, ObjFun_fcm ] = fcm(data, NumClusters, [ opt.ExpntQ opt.MaxIter opt.tol opt.verbose ]);
U_fcm_est_ScaledTo255 = U_fcm_est*255;

% do the RFCM clustering
figure;
[ centers_RFCM, U_RFCM_est, ObjFun_RFCM ] = robustFCM(data, ImgDim, NumClusters, opt, TypeOfExpt);
U_RFCM_est_ScaledTo255 = U_RFCM_est*255;

%% visualize the results for toy image

figure;
subplot 221; imagesc(TestImg); colormap gray; axis off;
title('Toy Image With noise');
subplot 222; imagesc(reshape(U_fcm_est_ScaledTo255(1,:),ImgDim)); colormap gray; axis off;
title('FCM: Class 1');
subplot 223; imagesc(reshape(U_fcm_est_ScaledTo255(2,:),ImgDim)); colormap gray; axis off;
title('FCM: Class 2');
subplot 224; imagesc(reshape(U_fcm_est_ScaledTo255(4,:),ImgDim)); colormap gray; axis off;
title('FCM: Class 4');

figure;
subplot 221; imagesc(TestImg); colormap gray; axis off;
title('Toy Image With noise');
subplot 222; imagesc(reshape(U_RFCM_est_ScaledTo255(1,:),ImgDim)); colormap gray; axis off;
title('RFCM: Class 1');
subplot 223; imagesc(reshape(U_RFCM_est_ScaledTo255(2,:),ImgDim)); colormap gray; axis off;
title('RFCM: Class 2');
subplot 224; imagesc(reshape(U_RFCM_est_ScaledTo255(4,:),ImgDim)); colormap gray; axis off;
title('RFCM: Class 4');




%% final seg

% % visualize the results, if requested
[ ~, FCMidx ] = max(U_fcm_est);
[ ~, RFCMidx ] = max(U_RFCM_est);
fcm_seg = reshape(FCMidx,ImgDim);
rfcm_seg = reshape(RFCMidx,ImgDim);
for kk = 1 : NumClusters
    fcm_seg(fcm_seg == kk ) = centers_fcm(kk);
    rfcm_seg(rfcm_seg == kk ) = centers_RFCM(kk);
end

figure; clf;
subplot 151; imagesc(TestImg); title('Toy Image');colormap gray; axis off;
subplot 152; imagesc(fcm_seg); title('FCM Seg');colormap gray; axis off;
subplot 153; imagesc(rfcm_seg); title('RFCM Seg'); colormap gray; axis off;
subplot 154; imagesc(reshape(U_fcm_est_ScaledTo255(1,:),ImgDim)); colormap gray; axis off;
title('FCM: U');
subplot 155; imagesc(reshape(U_RFCM_est_ScaledTo255(2,:),ImgDim)); colormap gray; axis off;
title('RFCM: U');
colormap gray;


