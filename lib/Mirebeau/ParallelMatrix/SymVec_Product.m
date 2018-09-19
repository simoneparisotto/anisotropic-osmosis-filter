% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Computes the product of a symmetric matrix field with a vector field, yielding a vector field.
% Symmetric matrix format : lower triangular part given serially.



function v=SymVec_Product(sym,u)
    Dimension = size(u,1);
    nPoints = size(u,2);
    SymDimension = Dimension*(Dimension+1)/2;
    assert(all(size(sym)==[SymDimension,nPoints]));
    v=zeros([Dimension,nPoints]);
    
    klCounter=1;
    for k=1:Dimension
        for l=1:k
            v(k,:) = v(k,:) + sym(klCounter,:).*u(l,:);
            if k~=l
                v(l,:) = v(l,:)+sym(klCounter,:).*u(k,:);
            end
            klCounter = klCounter+1;
        end
    end
end