function NbrHood = getNeighborhoodVolumetricSlow(jj,ImgDim,TypeOfNbrhood)
% slower version with for loops and if conditions

assert( ismember(numel(ImgDim), [2 3]), 'Only 2D or 3D neighbourhood queries are allowed at the moment.');

if numel(ImgDim) == 2 && strcmpi(TypeOfNbrhood,'4nbr')
    % subscripts into image
    [ s1 s2 ] = ind2sub(ImgDim,jj);
    
    % find the nbrhood: 4nbr
    SubDim1 = [ s1   s1    s1-1 s1+1 ];
    SubDim2 = [ s2-1 s2+1  s2-1 s2+1 ];

    % checking if these cross the image boundaries
    if s1==1 
        SubDim1( SubDim1 <=0) = 1;
    elseif s1 == ImgDim(1)
        SubDim1( SubDim1 >=ImgDim(1)) = ImgDim(1);
    end
    
    if s2==1 
        SubDim2( SubDim2 <=0) = 1;
    elseif s2 == ImgDim(2)
        SubDim2( SubDim2 >=ImgDim(2)) = ImgDim(2);
    end
    
    NbrHood    = sub2ind(ImgDim,SubDim1,SubDim2);
    
elseif numel(ImgDim) == 3 && strcmpi(TypeOfNbrhood,'6nbr')
    % subscripts into image
    [ s1 s2 s3 ] = ind2sub(ImgDim,jj);
    % find the nbrhood: 4nbr
    SubDim1 = [ s1   s1    s1-1 s1+1 s1   s1   ];
    SubDim2 = [ s2-1 s2+1  s2   s2   s2   s2   ];
    SubDim3 = [ s3   s3    s3   s3   s3-1 s3+1 ];
    % checking if these cross the image boundaries
    SubDim1( SubDim1 <=0) = 1;
    SubDim2( SubDim2 <=0) = 1;
    SubDim3( SubDim3 <=0) = 1;
    
    SubDim1( SubDim1 >=ImgDim(1)) = ImgDim(1);
    SubDim2( SubDim2 >=ImgDim(2)) = ImgDim(2);
    SubDim3( SubDim3 >=ImgDim(3)) = ImgDim(3);
    
    NbrHood = sub2ind(ImgDim,SubDim1,SubDim2,SubDim3);
end


end
