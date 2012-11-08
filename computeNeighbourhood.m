
function nbrhood = computeNeighbourhood(data,ImgDim,TypeOfNbrhood)

if strcmpi(TypeOfNbrhood,'4nbr')
    numnbrs = 4;
elseif strcmpi(TypeOfNbrhood,'6nbr')
    numnbrs = 6;
end

nbrhood = nan(length(data),numnbrs);
for ii =1 :length(data)
    nbrhood(ii,:)=getNeighborhoodVolumetricSlow(ii,ImgDim,TypeOfNbrhood);
end

end
