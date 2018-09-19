% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

function sdpPow = SymmetricMatrix3Power(sdp,exponent)
assert(ndims(sdp)==2);
assert(size(sdp,1)==6);
  s=size(sdp,2);
  
  d=SymmetricMatrix3Eigenvalues(sdp);
  e=d.^exponent;
  P=InterpolatingPolynomial3(d,e);    
  sdp2=SymmetricMatrix3Square(sdp);
  
  sdpPow = zeros([6,s]);
  for i=1:6
    sdpPow(i,:)=P(1,:)*(i==1||i==3||i==6)+P(2,:).*sdp(i,:)+P(3,:).*sdp2(i,:);
  end
end