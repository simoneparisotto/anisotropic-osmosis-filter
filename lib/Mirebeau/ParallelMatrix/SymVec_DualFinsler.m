% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% Computes the dual to a Finsler metric field, of the form v -> sqrt(v.m.v) - w.v

function [m1,w1] = SymVec_DualFinsler(m0,w0)
    s = Sym_Inverse(m0 - Vec_SelfOuterProduct(w0));
    w1 = SymVec_Product(s,w0);
    m1 = ScalVec_Product( 1+VecVec_ScalarProduct(w0,w1), s);    
end