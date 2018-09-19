% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Multiplies a vector field by a scalar field, yielding a vector field

function v=ScalVec_Product(a,u)
    Dimension = size(u,1);
    nPoints = size(u,2);    
    assert(all(size(a)==[1,nPoints]));
    
    v=repmat(a,[Dimension,1]).*u;
end