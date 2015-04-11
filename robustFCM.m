function [ centers, U, ObjFun, seg ] = robustFCM(data, ImgDim, NumClusters, opt, TypeOfExpt, varargin)
%{
 Robust Fuzzy C-means clustering
   as presented in Dzung Pham "Spatial Models for Fuzzy Clustering", CVIU, 2001

data        - NumSamples X NumFeaturesPerSample
NumClusters - number of clusters to make
ImgDim      - 2x1 or 3x1 vector specifying the dimensions of the image
                being segmented, to help find the 4- or 6- or 8-nbrhood for
                each voxel
options     - a structure containing the following optional fields
    options.ExpntQ       - exponent for the matrix U     (Default: 2)
    options.beta    - trade-off between spatial smoothing and FCM objective function
    options.MaxIter - max. number of iterations     (Default: 100)
    options.tol     - tolerance/threshold for min improvement per iteration (Default: 1e-5)
    options.verbose - logical variable, whether to display the info per each iteration (Default: 1)

Implemented by Pradeep Reddy Raamana. 
- Comments and bug reports are welcome. 
- No warranty whatsoever, implied or otherwise.

%}

assert(nargin==5 || nargin == 6, 'invalid number of arguments - 5/6 only');
if nargin < 5
    TypeOfExpt = 'MriTissueSeg';
end

default_options.ExpntQ = 2;
default_options.MaxIter = 100;
default_options.beta = 75 ;
default_options.TypeOfNbrhood = '4nbr';
default_options.tol = 1e-5;
default_options.verbose = true;

if nargin < 3
    opt = default_options;
end

ExpntQ = opt.ExpntQ;
beta = opt.beta;

TypeOfNbrhood = opt.TypeOfNbrhood;

MaxIter = opt.MaxIter ;
tol = opt.tol;


nbrhood.given=false;
if nargin > 5
    nbrhood.given=true;
    nbrhood.data=varargin{1};
else
    nbrhood.data = computeNeighbourhood(data,ImgDim,TypeOfNbrhood);
    nbrhood.given=true;
end
nbrhood.TypeOfNbrhood=opt.TypeOfNbrhood;

NumSamples  = size(data,1);
NumFeatures = size(data,2);
ObjFun = nan(MaxIter,1);

U = initRobustFCM(NumClusters,NumSamples);
centers = initCentersRobustFCM(NumClusters, NumFeatures, TypeOfExpt);

for iter =1 : MaxIter
    
    [ U, centers, ObjFun(iter) ] = updateRobustFCM(data, U, centers, NumClusters, ExpntQ, beta, nbrhood.data);
    
    % % detailed info per iteration
    if opt.verbose
        if size(centers,2)==1
            StrCenters = '';
            for cc = 1:length(centers)
                StrCenters = [ StrCenters sprintf('%6.3f   ',centers(cc)) ]; %#ok<AGROW>
            end
            fprintf('\nIter %3d: objective fun = %8.4e; %s',iter, ObjFun(iter),StrCenters);
            
        else
            fprintf('\nIter %3d: objective fun = %8.4e ',iter, ObjFun(iter));
            display(centers);
        end
    end
    
    % % visualize the results, if requested
    if opt.visualize
        [ ~, IdxMax ] = max(U);
        seg = reshape(IdxMax,ImgDim);
        for kk = 1 : NumClusters
            seg(seg == kk ) = centers(kk);
        end
        
        figure(gcf); clf;
        imagesc(seg);
        colormap gray;
        colorbar('YTick',centers);
        title(sprintf('Iter %d',iter));
        drawnow;
    end
    
    
    % convergence check
    if iter > 1 && ( abs(ObjFun(iter)-ObjFun(iter-1)) < tol )
        fprintf('\n Algorithm converged after %d iterations with the following cluster centers \n',iter);
        display(centers);
        break;
    end
    
end

ActNumIter = iter;	% Actual number of iterations
if ActNumIter < MaxIter
    ObjFun(ActNumIter+1:MaxIter) = [];
end

end

function  [ U_new, centers_new, ObjFun_new ] = updateRobustFCM(data, U,centers,NumClusters,ExpntQ, beta, nbrhood)
%
% compute membership functions
% Numerator in equation 8 for U ( page 289 ) in Dzung Pham paper
CalcNumU = @(jj,kk,sncm) realpow(abs(data(jj,:) - centers(kk,:))^2  + beta*sncm, -1/(ExpntQ-1) ) ;
calcSNCMq = @(u,q,othr,nbrd) sum(sum(u(othr,nbrd).^q)) ;

ClassSet = 1 : NumClusters;

U_new = nan(size(U)); % must be of size NumClusters x NumSamples

for jj =1 : size(data,1) % for each sample
    for kk =1 : NumClusters
        
        sncmq = calcSNCMq(U,ExpntQ,ClassSet(ClassSet~=kk),nbrhood(jj,:));
        NumU = CalcNumU(jj,kk,sncmq);
        
        % DenU = sum of NumU over all the classes
        DenU = NumU; % for class kk
        % computing it for remaining classes
        OtherClasses = ClassSet(ClassSet~=kk) ;
        for Othr = OtherClasses
            % notice the different index Othr used to loop thru classes.
            sncmq = calcSNCMq(U,ExpntQ,ClassSet(ClassSet~=Othr),nbrhood(jj,:));
            DenU = DenU + CalcNumU(jj,Othr,sncmq);
        end
        
        U_new(kk,jj) = NumU / DenU;
    end
end

% normalizing the memberships to make sure
% they sum to unity over clases for each sample point
U_new = U_new ./ ( ones(NumClusters, 1)*sum(U_new) );

%
% compute centroids
%
centers_new = nan(size(centers)); % must be of size NumClusters x NumFeaturesPerSample
for kk =1 : NumClusters
    NumCtr = sum( ( U_new(kk,:).^ExpntQ ) * data );
    DenCtr = sum( ( U_new(kk,:).^ExpntQ ) ) ;
    centers_new(kk,:) = NumCtr ./ DenCtr ;
end

%
%
% compute the Robust FCM objective function
%   CostRobustFCM = CostFCM + beta/2*Penalty

CostFCM = 0;
for kk =1 : NumClusters
    CostFCM = CostFCM + sum( (U_new(kk,:)'.^ExpntQ) .* ( (data - repmat(centers_new(kk,:),size(data,1),1)).^2 ) );
end

% Penalty For Spatial Nbrhood
Penalty = 0;
if beta ~= 0
    
    for kk =1 : NumClusters
        for jj = 1 : size(data,1)
            % notice the use of U_new
            sncmq_new = calcSNCMq(U_new,ExpntQ,ClassSet(ClassSet~=kk),nbrhood(jj,:));
            Penalty = Penalty + U_new(kk,jj).^ExpntQ * sncmq_new;
        end
    end
    
end

ObjFun_new = CostFCM + 0.5*beta*Penalty;

end


% function sncm = calcSNCMq(U,ExpntQ,OtherClasses,nbrhood)
% % replaced by an anonymous function as it can be reduced to a single
% % statement to avoid for loop
% % computes Sum Of Neighbouring Classes Memberships raised to the power ExpntQ, denoted by SNCMq
% sncm = 0;
% for mm = OtherClasses
%     sncm = sncm + sum( U(mm,nbrhood).^ExpntQ );
% end
% end

% function sncm = calcSNCMq(U,jj,kk,ExpntQ,ClassSet,ImgDim,TypeOfNbrhood)
% % computes Sum Of Neighbouring Classes Memberships raised to the power ExpntQ, denoted by SNCMq
% sncm = 0;
% Nbrhood     = getNeighborhoodVolumetricSlow(jj,ImgDim,TypeOfNbrhood);
% OtherClasses = ClassSet(ClassSet~=kk) ;
% for mm = OtherClasses
%     sncm = sncm + sum( U(mm,Nbrhood).^ExpntQ );
% end
% end


function centers = initCentersRobustFCM(NumClusters, NumFeatures,TypeOfExpt)

if strcmpi(TypeOfExpt,'MRITISSUESEG')
    if NumClusters == 3 && NumFeatures == 1
        centers = [ 180 ; 110; 5]; % close to ground truth
        % centers = [ 200 ; 80 ; 5];
        display('Smoothed Kernel Estimator is yet to be implemented');
        display('  to automatically find the number of modes in the MRI intensity histogram');
    end
elseif strcmpi(TypeOfExpt,'ToyImgDzungPham')
    if NumClusters == 4 && NumFeatures == 1
        centers = [ 0; 22 ; 53 ; 70 ];
        %         centers = [ 10; 32 ; 55 ; 90 ];
    end
    
    if NumClusters == 3 && NumFeatures == 1
        centers = [ 22 ; 53 ; 77 ];
    end
end

% centers = uint8(centers);

end

function U = initRobustFCM(NumClusters,NumSamples)

U = rand(NumClusters, NumSamples);
col_sum = sum(U);
% normalizing them to make sure the memberships sum to unity for each sample
U = U./col_sum(ones(NumClusters, 1), :);

end
%
