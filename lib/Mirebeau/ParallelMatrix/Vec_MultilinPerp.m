% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% (d-1)-Multilinear antisymmetric product, obeying 
% <u1 x u2 x ... x u_{d-1}, u_d> = det(u_1,...,u_d)

function perp = Vec_MultilinPerp(u)
    Dimension = size(u,1);
    nPoints = size(u,2);
    if Dimension~=2
        assert(all(size(u)==[Dimension,nPoints,Dimension-1]));
    else 
        assert(all(size(u)==[Dimension,nPoints]));
    end
    
    perp = zeros([Dimension,nPoints]);
    switch Dimension
        case 1
            perp(:)=1;
        case 2
            perp(1,:) = -u(2,:,1);
            perp(2,:) =  u(1,:,1);
        case 3
            cc = [1,2,3; 2,3,1; 3,1,2];
            for i=1:3
                c=cc(i,:);
                perp(c(2),:) = u(c(2),:,1).*u(c(3),:,2) - u(c(3),:,1).*u(c(2),:,2);
            end
        otherwise
            disp('Unsupported Vec_MultilinPerp in dimension >= 4 sorry');
            assert(false);
    end


end