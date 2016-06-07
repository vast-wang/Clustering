function [W, V] = NormalizeWV(K, W, V, NormV, Norm, numC)

if Norm == 2 % normalize to have unit norm L_2
    if NormV % normalize V and scale W
        norms = max(1e-15,sqrt(sum(V.^2,1)))';
        V = V*spdiags(norms.^-1,0,numC,numC);
        W = W*spdiags(norms,0,numC,numC);
    else % normalize W and scale V
        norms = max(1e-15,sqrt(sum(W.*(K*W),1)))';
        W = W*spdiags(norms.^-1,0,numC,numC);
        V = V*spdiags(norms,0,numC,numC);
    end
else % normalize to have unit norm L_1
    if NormV
        norms = max(1e-15,sum(abs(V),1))';
        V = V*spdiags(norms.^-1,0,numC,numC);
        W = W*spdiags(norms,0,numC,numC);
    else
        norms = max(1e-15,sum(W.*(K*W),1))';
        W = W*spdiags(norms.^-1,0,numC,numC);
        V = V*spdiags(norms,0,numC,numC);
    end
end