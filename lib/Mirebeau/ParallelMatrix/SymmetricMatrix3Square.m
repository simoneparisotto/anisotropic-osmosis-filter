% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

%Compute in parallel the square of matrices defined by c.
%Format : c(:,i) represents the i-th 3x3 symmetric matrix, 
%with entries [xx,xy,yy,xz,yz,zz]

function c2 = SymmetricMatrix3Square(c)
    assert(size(c,1)==6);
    assert(ndims(c)==2);
    s=size(c,2);

  c2 = zeros([6,s]);
  c2(1,:)=c(1,:).^2     +c(2,:).^2     +c(4,:).^2;
  c2(2,:)=c(1,:).*c(2,:)+c(2,:).*c(3,:)+c(4,:).*c(5,:);
  c2(3,:)=c(2,:).^2     +c(3,:).^2     +c(5,:).^2;
  c2(4,:)=c(1,:).*c(4,:)+c(2,:).*c(5,:)+c(4,:).*c(6,:);
  c2(5,:)=c(2,:).*c(4,:)+c(3,:).*c(5,:)+c(5,:).*c(6,:);
  c2(6,:)=c(4,:).^2     +c(5,:).^2     +c(6,:).^2;
end