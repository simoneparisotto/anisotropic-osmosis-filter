% Copyright Jean-Marie Mirebeau, 2017. jm(dot)mirebeau(at)gmail(dot)com

% This function takes as input 
% - s : positive definite tensor field
% - w : vector field
% - options, with fields specifying : 
%       * dims : the grid dimensions.
%       * gridscale : the pixel scale.

% It assembles the finite differences matrix of 
% - div(s (grad u - w u))

% The second order operator discretization is based on the AD-LBR scheme
% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

% The first order operator discretization is upwind, and uses the same stencils.

function [Rows,Cols,Vals] = ConvectionDiffusionSparseMatrix(s,w,options)

% --- Get meta parameters ---
dims = options.dims;
Dimension = size(dims,2);
SymDimension = Dimension*(Dimension+1)/2;
nPoints = prod(dims);
h = options.gridScale;
upwind = isfield(options,'upwind') && options.upwind;


% --- Construct the stencils ---
[offsets,weights] = Sym_Decomposition(s); % Voronoi's first reduction
nOffsets = size(offsets,3); % Number of scheme offsets at each point (=SymDimension here)
indices = zeros([1,nPoints,nOffsets,2]); % Linear indices of neighbors
indices(:,:,:,1) = GridIndices(dims, offsets); 
indices(:,:,:,2) = GridIndices(dims,-offsets);

% Assemble the sparse matrix
Rows = [];
Cols = [];
Vals = [];

if numel(w)>0
    % Order one : - div( u w), which is the matrix of (v,u) -> <grad(v), u w>
    assert(all(size(w)==[Dimension,nPoints]));

    for iOffset = 1:nOffsets
        coef = h^(-1) * VecVec_ScalarProduct(w,offsets(:,:,iOffset)) .* weights(1,:,iOffset);
        
        if upwind
            % Upwind first order finite differences
            neigh = zeros([1,nPoints]);
            neigh(coef>=0) = indices(1,coef>=0,iOffset,1);
            neigh(coef<0)  = indices(1,coef<0, iOffset,2);
            
            coef = abs(coef);
            pos = 1:nPoints;
            
            % Keep only 'inbox' terms, that fall inside the domain.
            inbox   = (neigh ~= 0);
            pos     = pos(inbox);
            neigh   = neigh(inbox);
            coef    = coef(inbox);
            
            % Assemble the sparse matrix
            Rows = [Rows,   pos,    neigh];
            Cols = [Cols,   pos,    pos];
            Vals = [Vals,   coef,   -coef];
        else
            % Centered second order finite differences
            coef = coef/2;
            % Gather neighbor indices
            neigh1 = indices(1,:,iOffset,1);
            neigh2 = indices(1,:,iOffset,2);
            pos=1:nPoints;
            
            % Fallback to  off-centered, first order, on boundary
            out = (neigh1==0);
            neigh1(out) = pos(out);
            coef(out) = 2*coef(out);
            
            out = (neigh2==0);
            neigh2(out) = pos(out);
            coef(out) = 2*coef(out);
            
            % Assemble the sparse matrix
            Rows = [Rows,   neigh1, neigh2];
            Cols = [Cols,   pos,       pos];
            Vals = [Vals,   -coef,    coef];

        end
    end
end

if numel(s)>0
    %Order two : - div( s grad(u)), which is the matrix of (v,u) -> <grad(v), s grad(u)>
    assert(all(size(s)==[SymDimension,nPoints]));
    for iOffset = 1:nOffsets
        for i=1:2
            coef = h^(-2) * weights(1,:,iOffset);
            neigh = indices(1,:,iOffset,i);
            pos = 1:nPoints;
                        
            % Keep only 'inbox' terms, that fall inside the domain. 
            inbox   = (neigh~=0);
            pos     = pos(inbox);
            neigh   = neigh(inbox);
            coef    = coef(inbox);
            
            % Assemble the sparse matrix
            Rows = [Rows,   pos,    pos,    neigh,  neigh];
            Cols = [Cols,   pos,    neigh,  pos,    neigh];
            coef = coef/2; % Take into account the opposite offset
            Vals = [Vals,   coef,   -coef,  -coef,   coef];
        end
    end
end

inbox   = Rows~=0 & Cols~=0 & Vals~=0; %Out of box index or null matrix coefficient.
Rows = Rows(inbox);
Cols = Cols(inbox);
Vals = Vals(inbox);

end