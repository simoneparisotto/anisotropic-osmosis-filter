function norms = SymVecVec_EvalFinsler(m,w,u)
    nPoints=size(w,2);
    Dimension = size(w,1);
    SymDimension = Dimension*(Dimension+1)/2;
    assert(all(size(m)==[SymDimension,nPoints]));
    assert(all(size(w)==[Dimension,nPoints]));
    assert(all(size(w)==[Dimension,nPoints]));
    
    norms = sqrt(VecSymVec_ScalarProduct(u,m,u)) - VecVec_ScalarProduct(w,u);
end