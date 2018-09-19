% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% This file is a parallel implementation of Selling's algorithm,
% for computing obtuse superbases of 3x3 SDP matrices.
% Key ingredient of the anisotropic diffusion scheme described in

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

function [sb,nIter] = Sym_ObtuseSuperbase(m,maxIter)
    Dimension = floor(sqrt(2*size(m,1)));
    nPoints = size(m,2);
    SymDimension = Dimension*(Dimension+1)/2;
    
    assert(size(m,1)==SymDimension);
    assert(Dimension<=3);
    
    if nargin<2
        maxIter=200;
    end
    
    % Initialize to canonical superbase
    sb = zeros([Dimension,nPoints,Dimension+1]);
    for k=1:Dimension
        sb(k,:,k)=1;
        sb(k,:,Dimension+1)=-1;
    end
    
    % Run Selling's algorithm
    if Dimension==1
        [sb,nIter]=Selling1(m,sb,maxIter);
    elseif Dimension==2
        [sb,nIter]=Selling2(m,sb,maxIter);
    elseif Dimension==3
        [sb,nIter]=Selling3(m,sb,maxIter);
    end
end

function [sb,nIter]=Selling1(~,sb,~)
    nIter=0;
end


function [sb,nIter] = Selling2(m,sb,maxIter)
    cc = [1,2,3; 2,3,1; 3,1,2];
    scalNeg=0;
    nIter=0;
    while scalNeg<3
        c=cc(1+mod(nIter,3),:);
        scal = VecSymVec_ScalarProduct(sb(:,:,c(1)),m,sb(:,:,c(2)));
        scal = scal>0;
        if max(scal(:)) %Superbase update -b1, b2, b1-b2
            sb(:,scal,c(3))= sb(:,scal,c(1))-sb(:,scal,c(2));
            sb(:,scal,c(1))=-sb(:,scal,c(1));
            scalNeg=0;
        else
            scalNeg=scalNeg+1;
        end
        
        nIter=nIter+1;
        if nIter==maxIter;
            disp('Warning : Selling2 algorithm unterminated')
            break
        end
    end
end

function [sb,nIter] = Selling3(m,sb,maxIter)
    cc = [1,2,3,4; 1,3,2,4; 1,4,2,3; 2,3,1,4; 2,4,1,3; 3,4,1,2];
    scalNeg=0;
    nIter=0;
    while scalNeg<6
        c=cc(1+mod(nIter,6),:);
        scal = VecSymVec_ScalarProduct(sb(:,:,c(1)),m,sb(:,:,c(2)));
        scal = scal>0;
        if max(scal(:)) %Superbase update -b1, b2, b3+b1, b4+b1
            sb(:,scal,c(3))= sb(:,scal,c(3))+sb(:,scal,c(1));
            sb(:,scal,c(4))= sb(:,scal,c(4))+sb(:,scal,c(1));
            sb(:,scal,c(1))=-sb(:,scal,c(1));
            scalNeg=0;
        else
            scalNeg=scalNeg+1;
        end
        
        nIter=nIter+1;
        if nIter==maxIter;
            disp('Warning : Selling3 algorithm unterminated')
            break
        end
    end
end