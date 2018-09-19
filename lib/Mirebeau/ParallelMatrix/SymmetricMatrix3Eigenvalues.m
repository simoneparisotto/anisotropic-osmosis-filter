% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Real eigenvalues, sorted by increasing magnitude.
% TO DO : take advantage of the matrix symmetry (?)
% Matrix format : xx,xy,yy,xz,yz,zz

function e = SymmetricMatrix3Eigenvalues(sdp)
assert(ndims(sdp)==2);
assert(size(sdp,1)==6);

  s=size(sdp,2);
  % Using Eig3 (Bruno Luong) for solving multiple eigenvalues.
  % Symmetric tensors
  st = zeros([3,3,s]);
  st(1,1,:)=sdp(1,:);
  st(1,2,:)=sdp(2,:);
  st(2,1,:)=sdp(2,:);
  st(2,2,:)=sdp(3,:);
  st(1,3,:)=sdp(4,:);
  st(3,1,:)=sdp(4,:);
  st(2,3,:)=sdp(5,:);
  st(3,2,:)=sdp(5,:);
  st(3,3,:)=sdp(6,:);
  
  d=eig3(st);
  d=real(d);
  % eigenvalues need to be sorted...
  e = zeros(3,s);
  e(1,:)=min(d,[],1);
  e(3,:)=max(d,[],1);
  e(2,:)=sum(d,1)-e(1,:)-e(3,:);
end