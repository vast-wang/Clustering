%% Calculate objective function
function [obj_NMF, obj_Lap, obj_VVc, obj_Pi] = CalculateObj(data, W, V, L,Vcon, options, pai, view_num)

tempNMF = zeros(1,view_num);
tempLap = zeros(1,view_num);
tempVVc = zeros(1,view_num);
tempPi = zeros(1,view_num);
for i = 1:view_num
    KW = data{i}*W{i};
    KWV = KW*V{i}';
    VL = V{i}'*L{i};
    VLV = VL*V{i};
    tempNMF(i) = tempNMF(i) + norm(data{i}-KWV,'fro')^2;
    tempLap(i) = tempLap(i) + options.alpha*trace(VLV);

    VtVc = V{i} - Vcon;
    tempVVc(i) = tempVVc(i) + options.beta*pai(i)*norm(VtVc,'fro')^2;
    tempPi(i) = tempPi(i) + options.gamma*norm(pai,'fro')^2;
    clear KW KWV VL VLV;
end
obj_NMF = sum(tempNMF);
obj_Lap = sum(tempLap);
obj_VVc = sum(tempVVc);
obj_Pi = sum(tempPi);
end