function mInv = inv2(m)

    mInv = m\eye(size(m));
    
    mInv(abs(mInv)<1e-16) = 0;

end