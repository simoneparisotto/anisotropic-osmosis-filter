function scal = VecSymVec_ScalarProduct(u,m,v)
    scal = VecVec_ScalarProduct(u,SymVec_Product(m,v));
end