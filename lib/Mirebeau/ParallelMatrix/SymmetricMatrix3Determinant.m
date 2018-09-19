% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Matrix format : sdp(:,i) = [xx,xy,yy,xz,yz,zz]

function det = SymmetricMatrix3Determinant(sdp)
    assert(ndims(sdp)==2)
    assert(size(sdp,1)==6);
    det = sdp(1,:).*sdp(3,:).*sdp(6,:)+2*sdp(2,:).*sdp(5,:).*sdp(4,:) ...
    -sdp(1,:).*sdp(5,:).^2-sdp(3,:).*sdp(4,:).^2-sdp(6,:).*sdp(2,:).^2;
end