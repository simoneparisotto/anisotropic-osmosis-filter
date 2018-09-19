% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Computes the inverse of a field of positive definite tensors

% !! TODO !! Specializations in dimension 1,2,3 ?

function inv = Sym_Inverse(mat)
    Dimension = floor(sqrt(2*size(mat,1)));
    nPoints = size(mat,2);
    
    SymDimension = Dimension*(Dimension+1)/2;
    assert(size(mat,1)==SymDimension);
    
    inv = zeros([SymDimension,nPoints]);
    klCounter=1;
    for k=1:Dimension
        u = zeros([Dimension,nPoints]);
        u(k,:)=1;
        u=SymVec_Solve(mat,u);
        for l=1:k
            inv(klCounter,:) = u(l,:);
            klCounter=klCounter+1;
        end
    end
end