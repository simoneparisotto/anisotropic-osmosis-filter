function [V,WW] = compute_W(V,b,umask)

ONE   = ones(size(V,1),size(V,2));
ZERO  = zeros(size(V,1),size(V,2));
MASK2 = repmat(umask(:,:,1),1,1,2)==1;
MASK4 = repmat(umask(:,:,1),1,1,4)==1;

B1 = b{1}.*ONE;
B2 = b{2}.*ONE;

W1 = cat(3,ONE,ZERO,ZERO,ONE);
WW = cat(3, B1.^2.*V(:,:,1).^2         + B2.^2.*V(:,:,2).^2,...
            B1.^2.*V(:,:,1).*V(:,:,2)  - B2.^2.*V(:,:,1).*V(:,:,2),...
            B1.^2.*V(:,:,1).*V(:,:,2)  - B2.^2.*V(:,:,1).*V(:,:,2),...
            B1.^2.*V(:,:,2).^2         + B2.^2.*V(:,:,1).^2);
        
WW(MASK4) = W1(MASK4);
V(MASK2)  = 0;

end