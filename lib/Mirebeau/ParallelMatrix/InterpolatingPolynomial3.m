% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com
  
% Interpolating polynomial P(:,i) from d(:,i) to e(:,i).
% Polynomial degree : 2, interpolates three values.
% TO DO : handle case where several entries are equal

function P = InterpolatingPolynomial3(d,e)
  assert(ndims(d)==2)
  assert(all(size(d)==size(e)));
  assert(size(d,1)==3);
  
  P=zeros(size(d));
  P(3,:) = d(2,:).* e(1,:) - d(3,:).* e(1,:) - d(1,:).*e(2,:) + d(3,:).*e(2,:) + d(1,:).*e(3,:) - d(2,:).*e(3,:);
  P(2,:) = -(d(2,:).^2).*e(1,:) + (d(3,:).^2).* e(1,:) + (d(1,:).^2).*e(2,:) - (d(3,:).^2).*e(2,:) - (d(1,:).^2).*e(3,:) + (d(2,:).^2).* e(3,:);
  P(1,:) = (d(2,:).^2).* d(3,:).* e(1,:) - d(2,:).*( d(3,:).^2).* e(1,:) - (d(1,:).^2).* d(3,:).* e(2,:) + ...
  d(1,:).* (d(3,:).^2) .* e(2,:) + (d(1,:).^2).* d(2,:).* e(3,:) - d(1,:).* (d(2,:).^2).* e(3,:);

  %divide by discriminant
  dprod = (d(1,:)-d(2,:)).*(d(1,:)-d(3,:)).*(d(2,:)-d(3,:));
  for i=1:3 P(i,:)=P(i,:)./dprod; end;

end