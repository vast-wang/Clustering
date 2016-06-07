function [ac,nmi_value,f_score,indic] = printResult(X, label, K, clusteringFlag)

if clusteringFlag == 1
    indic = litekmeans(X, K, 'Replicates',20);
else
    [~, indic] = max(X, [] ,2);
end
[ac, nmi_value, f_score, cnt] = CalcMeasures(label, indic); % including bestMap()
disp(sprintf('ac:%0.4f\t%d/%d\tnmi:%0.4f\tFscore:%0.4f\t', ac, cnt, length(label), nmi_value,f_score));