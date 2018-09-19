% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Indices of neighbors of gridpoints.

function indices = GridIndices(dims,offsets)
    Dimension = size(offsets,1);
    nPoints = size(offsets,2);
    nOffsets = size(offsets,3);
    
    assert(size(dims,2)==Dimension);
    assert(prod(dims)==nPoints);
    
    indices = zeros([1,nPoints,nOffsets]);
    subscripts = zeros([Dimension,nPoints]);
    inbox = ones([1,nPoints]);
    for iOffset=1:nOffsets
        
        % Dimension dependent switch is only due to Matlab's syntax.
        switch Dimension
            case 1
                subscripts(1,:) = ind2sub(dims,1:nPoints);
            case 2
                [subscripts(1,:),subscripts(2,:)] = ind2sub(dims,1:nPoints);
            case 3
                [subscripts(1,:),subscripts(2,:),subscripts(3,:)] = ind2sub(dims,1:nPoints);
            otherwise
                disp('Sorry, Grid indices only supports dimension up to 3'); 
                assert(false);
        end
        
        subscripts = subscripts+offsets(:,:,iOffset);
        inbox(:)=1;
        for k=1:Dimension
            inbox = inbox & 0<subscripts(k,:) & subscripts(k,:)<=dims(k);
        end
        
        switch Dimension
            case 1
                indices(1,inbox,iOffset) = sub2ind([dims,1], subscripts(1,inbox));
            case 2
                indices(1,inbox,iOffset) = sub2ind(dims, subscripts(1,inbox), subscripts(2,inbox));
            case 3
                indices(1,inbox,iOffset) = sub2ind(dims, subscripts(1,inbox), subscripts(2,inbox), subscripts(3,inbox));
            otherwise 
                assert(false);
        end
    end
end