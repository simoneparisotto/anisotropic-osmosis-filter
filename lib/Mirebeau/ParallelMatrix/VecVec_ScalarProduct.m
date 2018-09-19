% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Computes the scalar product of two vector fields, yielding a scalar field
function scal = VecVec_ScalarProduct(u,v)
    Dimension = size(u,1);
    nPoints = size(u,2);
    assert(all(size(v) == [Dimension,nPoints]));
    
    scal = zeros([1,nPoints]);
    for k=1:Dimension
        scal = scal + u(k,:).*v(k,:);
    end
end