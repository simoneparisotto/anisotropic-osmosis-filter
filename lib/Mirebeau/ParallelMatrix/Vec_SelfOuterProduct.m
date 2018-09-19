% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Computes the outer product of a vector field with itself, yielding a rank one symmetric matrix field.
% Symmetric matrix format : lower triangular part given serially.



function sym = Vec_SelfOuterProduct(vec)
    Dimension = size(vec,1);
    SymDimension = Dimension*(Dimension+1)/2;
    nPoints = size(vec,2);
    sym = zeros([SymDimension, nPoints]);
    klCounter=1;
    for k=1:Dimension
        for l=1:k
            sym(klCounter,:) = vec(k,:).*vec(l,:);
            klCounter = klCounter+1;
        end
    end
end