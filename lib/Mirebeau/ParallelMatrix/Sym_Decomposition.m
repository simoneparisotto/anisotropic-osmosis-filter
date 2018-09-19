% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Non-negative decomposition of a positive definite matrix, based on Obtuse superbases

function [offsets,weights] = Sym_Decomposition(m,sb)
    if nargin<2
        sb = Sym_ObtuseSuperbase(m);
    end
    
    Dimension = floor(sqrt(2*size(m,1)));
    nPoints = size(m,2);
    SymDimension = Dimension*(Dimension+1)/2;
    assert(all(size(m,1)==SymDimension));
    assert(all(size(sb)==[Dimension,nPoints,Dimension+1]));
    
    offsets = zeros([Dimension,nPoints,SymDimension]);
    weights = zeros([1,nPoints,SymDimension]);
    klCounter=1;
    for k=1:Dimension
        for l=(k+1):(Dimension+1)
            weights(1,:,klCounter) = -VecSymVec_ScalarProduct(sb(:,:,k),m,sb(:,:,l));
            klExcept = 1:(Dimension+1);
            klExcept = klExcept(klExcept~=k & klExcept~=l);
            offsets(:,:,klCounter) = Vec_MultilinPerp(sb(:,:,klExcept));
            klCounter=klCounter+1;
        end
    end
end