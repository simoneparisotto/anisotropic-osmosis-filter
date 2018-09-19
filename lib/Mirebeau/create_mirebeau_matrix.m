function P3 = create_mirebeau_matrix(umat,umask,W,flag_upwind)

% Create axes
[x,y] = meshgrid(linspace(1,size(umat,1),size(umat,1)),linspace(1,size(umat,2),size(umat,2)));
x=x'; y=y'; % Avoid Matlab's YXZ ordering

%[x,y] = ndgrid(1:size(umat,1),1:size(umat,2));

% Set options
options.upwind    = flag_upwind;
options.dims      = size(umat);
options.gridScale = x(2,1)-x(1,1);

%%% Create tensor field (xx,xy,yy)
m = zeros([3,options.dims]);

%%% W = general (symmetric)
m(1,:,:) = W(:,:,1); 
m(2,:,:) = W(:,:,2); 
m(3,:,:) = W(:,:,4);

%%% Compute image gradient
[gradY,gradX] = gradient(umat,options.gridScale,options.gridScale); % Again, YXZ ordering
gradY = (gradY)./umat;
gradX = (gradX)./umat;

% %%% Compute image gradient
% umatextra = padarray(umat,[1,1],'replicate','pre');
% umat1 = (umatextra(1:end-1,:)+umatextra(2:end,:))./2;
% umat2 = (umatextra(:,1:end-1)+umatextra(:,2:end))./2;
% DU = grad_forward(umat1);
% gradX = -DU(:,1:end-1,1)./umat1(:,1:end-1);
% DU = grad_forward(umat2);
% gradY = DU(1:end-1,:,2)./umat2(1:end-1,:);


%%% SHADOW REMOVAL d==0
gradX(umask==0) = 0;
gradY(umask==0) = 0;

w = [shiftdim(gradX,-1); shiftdim(gradY,-1)];

% Reshape for input
nPoints = prod(options.dims);
m = reshape(m,[3,nPoints]);
w = reshape(w,[2,nPoints]);

%%% SYSTEM
[Rows,Cols,Vals] = ConvectionDiffusionSparseMatrix(m,w,options);
P3 = -sparse(Rows,Cols,Vals,nPoints,nPoints);