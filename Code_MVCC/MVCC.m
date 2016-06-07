function [Vall_, Veach_, obj_] = MVCC(data,truelabel,options,view_num,numC)
%function MVCC(data,truelabel,options,view_num,numC)
% Summary of this function goes here
%   Detailed explanation goes here

%% initialize ...
new_error = inf;
differror = options.error;
maxIter = options.maxIter;
minIter = options.minIter;
Rounds = options.Rounds;
alpha = options.alpha;
beta = options.beta;
gamma = options.gamma;
pai = options.pi;
%opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
opts = optimset('Algorithm','interior-point-convex','Display','off');
V_all = zeros();
obj = zeros();

Norm = 2; %normalize to have unit norm L_2; if Norm = 1, to have unit norm L_1 
NormV = 1;

%% kernel matrix X'*X£¬ calculate weight matrix W and Laplace matrix L
for i = 1:view_num
    % nums of dimensions and samples
    [d{i},n{i}] = size(data{i});
    
    % calculate kernel matrix Ker
    Kt = full(data{i}'*data{i});
    Ker{i} = max(Kt,Kt');
    
    % calculate weight matrix W
    TempWt = constructW(data{i}',options);
    if isfield(options,'NormWeight') && strcmpi(options.NormWeight,'NCW')
        D_mhalf = sum(TempWt,2).^-.5;
        D_mhalf = spdiags(D_mhalf,0,n{i},n{i});
        TempWt = D_mhalf*TempWt*D_mhalf;
        clear D_mhalf;
    end
    Wt{i} = TempWt;
    
    % calculate Laplance matrix L
    DCol = full(sum(Wt{i},2));
    D{i} = spdiags(DCol,0,speye(size(Wt{i},1)));
    TempL = D{i} - Wt{i};
    if isfield(options,'NormLap') && options.NormLap % note: now no options.NormLap
        D_mhalf = DCol.^-.5;
        tmpD_mhalf = repmat(D_mhalf,1,nSmp);
        TempL = (tmpD_mhalf.*TempL).*tmpD_mhalf';
        clear D_mhalf tmpD_mhalf;
        TempL = max(TempL, TempL');
    end
    L{i} = TempL;
    
    clear Kt TempWt DCol TempL;
end

iR = 0;
while iR < Rounds
    iR = iR + 1;
    %% initialize W{i}, V{i}, Vcon ...
        W = cell(1, view_num);
        V = cell(1, view_num);
        Vcon = zeros();
    %==================The first method================================
        for i= 1: view_num
            ilabels = zeros(n{i},view_num);
            ilabels(:,i) = kmeans(data{i}',numC,'replicates',20);
            %ilabels(:,i) = litekmeans(data{i}',numC,'Replicates',20);
            G = zeros(n{i},numC);
            for j=1:numC
                G(:,j)=(ilabels(:,i)==j*ones(n{i},1));
            end 
            V{i}=G+0.1*ones(n{i},numC);
            Dw=diag(sum(G,1))^-1;
            W{i}=(G+0.1*ones(n{i},numC))*Dw;
            % normlizing W{i}, V{i}
            [W{i},V{i}] = NormalizeWV(Ker{i}, W{i}, V{i}, NormV, Norm, numC);
            % initialize Vcon
            Vcon = Vcon + pai(i)*V{i};

            clear Kt Dw G tlabel;
        end

    %% Multiplicative optimization
    %     disp(sprintf('updating and optimizing...'));
    iter = 0;
    selectInit = 1;
    count = 0;
    while (iter < maxIter)
        iter = iter+1;
    %% update W{i}, V{i}, Vcon, pi_i
        VC = zeros();
        DVcV = zeros();
        for i = 1:view_num
        %===============  update W{i}  ====================
            KV = Ker{i}*V{i};
            VV = V{i}'*V{i};
            KWVV = Ker{i}*W{i}*VV;
            W{i} = W{i}.*(KV./max(KWVV,1e-9));
            clear KV VV KWVV;
        %===============  update V{i}  ====================
            KW = Ker{i}*W{i};
            WKW = W{i}'*KW;
            VWKW = V{i}*WKW;
            WtV = Wt{i}*V{i};
            DV = D{i}*V{i};
            KWp = KW + beta*pai(i)*Vcon + alpha*WtV;
            VWKWp = VWKW + alpha*DV + beta*pai(i)*V{i};
            V{i} = V{i}.*(KWp./max(VWKWp,1e-9));
            clear KW WKW VWKW WtV DV KWp VWKWp;
         % normlizing W{i}, V{i}
            [W{i},V{i}] = NormalizeWV(Ker{i}, W{i}, V{i}, NormV, Norm, numC);
         %===============  update Vcon  ====================
            VC = VC + pai(i)*V{i};
            % outside of for(), have Vcon = VC;
        end
         %===============  consensus representation  =============
            Vcon = VC;
            clear VC;

        for i = 1:view_num
        %===== calculate the difference between Vcon and each view =====
            DVcV(i) = norm(V{i}-Vcon,'fro')^2;
        end
        %===============  update pi_i  ====================
            %pai = cvxPi(DVcV, options.beta, view_num);
            if options.PiFlag
                H = options.gamma*eye(view_num,view_num);
                Aeq = ones(1,view_num);
                beq = 1;
                lb = zeros(view_num,1);
                pai = quadprog(H,DVcV,[],[],Aeq,beq,lb,[],[],opts);
                %pai = (1/view_num)*ones(view_num,1);
                clear DVcV;
            else
                pai = zeros(view_num,1);
            end

    %% Calculate objective function
        if (selectInit)
            [obj_NMF, obj_Lap, obj_VVc, obj_Pi] = CalculateObj(data, W, V, L,Vcon, options, pai, view_num);
            %new_obj = obj_NMF + obj_Lap + obj_VVc + obj_Pi;
            new_obj = obj_NMF + obj_Lap + obj_VVc;
            objhistory = new_obj;
            selectInit = 0;
        else
            [obj_NMF, obj_Lap, obj_VVc, obj_Pi] = CalculateObj(data, W, V, L,Vcon, options, pai, view_num);
            %new_obj = obj_NMF + obj_Lap + obj_VVc + obj_Pi;
            new_obj = obj_NMF + obj_Lap + obj_VVc;
        end
        if iter == 1
            Vall = Vcon;
            Veach = V;
        end
        if new_obj < objhistory(end)
            Vall = Vcon;
            Veach = V;
        end
        objhistory = [objhistory new_obj];
        new_error = abs(new_obj-objhistory(end-1));
        if (new_error < differror)
            count = count + 1;
            if (count == 5)
                iter = maxIter;
            end
        end
        clear new_obj;
    end
    % restart next round ...
    if iR == 1
        Vall_ = Vall;
        Veach_ = Veach;
        newobj_ = objhistory(end);
        obj_ = objhistory;
        newobj_
    end
    if objhistory(end) < newobj_
        Vall_ = Vall;
        Veach_ = Veach;
        newobj_ = objhistory(end);
        obj_ = objhistory;
        newobj_
        printResult(Vall, truelabel{1}, numC, options.clusteringFlag);
    end
    if iR < Rounds
        disp(sprintf('restart...'));
        clear objhistory;
    else
        disp(sprintf('get the result:'));
    end
end
end