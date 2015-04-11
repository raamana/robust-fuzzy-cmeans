
function nbrhood = computeNeighbourhood(data,ImgDim,TypeOfNbrhood)

assert(ismember(TypeOfNbrhood,{'4nbr','6nbr'}), 'Invalid choice. Must be ''4nbr'' or ''6nbr''.');

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
