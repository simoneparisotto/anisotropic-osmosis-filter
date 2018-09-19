% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Solves linear systems Ax=b involving a positive definite symmetric matrix field, and a vector field, yielding a vector field.
% !! Important : matrices must be positive definite !!

% Based on the conjugate gradient algorithm. (Requires positive definiteness of the tensors.)
% Symmetric matrix format : lower triangular part given serially.


function x = SymVec_Solve(A,b,xo)
    Dimension = size(b,1);
    nPoints = size(b,2);
    SymDimension = Dimension*(Dimension+1)/2;
    assert(all(size(A)==[SymDimension,nPoints]));
    
    if nargin<3
        xo=zeros(size(b));
    else
        assert(all(size(xo)==[Dimension,nPoints]));
    end
    
    % We expect Dimension to be small, typically 2 or 3 for applications of interest,
    % so the number of iterations of the CG algorithm is set to the dimension, 
    % in which case ConjugateGradient is an exact solver.
    nIter = Dimension;
    
    % Initial guess
    x=xo;
    
    % Preallocating vector fields p,q,r 
    p=zeros([Dimension,nPoints]);
    q=p;
    r=p;
    
    %Preallocating scalar fields alpha,beta,rr,roro
    alpha = zeros([1,nPoints]);
    beta = alpha;
    rr   = alpha;
    roro = alpha;
    pAp  = alpha;
    
    for k=1:nIter
        if k==1
            r=b-SymVec_Product(A,x);
            p=r;
            rr = VecVec_ScalarProduct(r,r);
        else 
            rr = VecVec_ScalarProduct(r,r);
            roro = VecVec_ScalarProduct(ro,ro);
            pos = roro>0;
            beta(pos)  = rr(pos)./roro(pos);
            beta(~pos) = 0.;
            p=r+ScalVec_Product(beta,p);
        end
        q=SymVec_Product(A,p);
        pAp = VecVec_ScalarProduct(p,q);
        pos = pAp>0;
        alpha(pos) = rr(pos)./pAp(pos);
        alpha(~pos) = 0.;
        x = x+ScalVec_Product(alpha,p);
        ro=r;
        r=r-ScalVec_Product(alpha,q);
    end
end