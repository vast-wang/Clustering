clear all
clc
close all
%% set and add the path ...
addpath(genpath('Dataset/'));
addpath(genpath('Results/'));
datadir = 'Dataset/';
resultdir = 'Results/';
dataname = {'3sources','BBC4view_685','BBCSport2view_544','20newsgroups', 'WikipediaArticles'};
numdata = length(dataname);
for cdata = 1:numdata % choose the dataset
%% parameter setting
options = [];
options.WeightMode = 'HeatKernel';
options.NormWeight = 'NCW';
options.k = 5;

options.maxIter = 200;
options.minIter = 50;
options.Rounds = 5;
options.nRepeat = 1;
options.error = 1e-5;
options.clusteringFlag = 0;
options.alpha = 100;
options.beta = 100;
options.gamma = 10;
options.pi = zeros();
options.PiFlag = 1;
%% read dataset
idata = cdata; % choose the dataset
dataf = [datadir, cell2mat(dataname(idata))];
load (dataf);
disp(sprintf('Dataset: %s',cell2mat(dataname(idata))));

%% normalize data matrix
view_num = length (data);
for i = 1:view_num
    data{i} = NormalizeData(data{i},2);
    options.pi(i) = 1/view_num;
end
numC = length(unique(truelabel{1}));

Acc = zeros();
nmi_V = zeros();
f_score = zeros();
Acce = zeros();
nmi_Ve = zeros();
f_each = zeros();
Vcon = cell(1,options.nRepeat);
R_label = cell(1,options.nRepeat);
obj = cell(1,options.nRepeat);

%% Repeat ...
for iter = 1:options.nRepeat

    % multi-view clustering method: MVCC
    disp(sprintf('Repeat %d times, in the %d-th times......',options.nRepeat,iter));
    
    [Vcon{iter}, Veach{iter}, obj{iter}] = MVCC(data,truelabel,options,view_num,numC);
    [Acc(iter), nmi_V(iter), f_score(iter), R_label{iter}]=printResult(Vcon{iter}, truelabel{1}, numC, options.clusteringFlag);
    for i=1:view_num
        [Acce(i,iter),nmi_Ve(i,iter),f_each(i,iter),result_label]=printResult(Veach{iter}{i}, truelabel{i}, numC, options.clusteringFlag);
    end
end
    Result(1,:) = Acc;
    Result(2,:) = nmi_V;
    Result(3,:) = f_score;
    Result(4,1) = mean(Acc);
    Result(4,2) = mean(nmi_V);
    Result(4,3) = mean(f_score);
    Result(5,1) = std(Acc);
    Result(5,2) = std(nmi_V);
    Result(5,3) = std(f_score);
for i=1:view_num
    Result(6,i) = mean(Acce(i,:));
    Result(6,i+view_num) = mean(nmi_Ve(i,:));
    Result(6,i+2*view_num) = mean(f_each(i,:));
    Result(7,i) = std(Acce(i,:));
    Result(7,i+view_num) = std(nmi_Ve(i,:));
    Result(7,i+2*view_num) = std(f_each(i,:));
end
save([resultdir,char(dataname(idata)),'_result.mat'],'Result','Vcon','Veach','obj');
clear Result Vcon Veach obj;
%% Print the objective function value
% plot(obj{1}(1,1:end),'LineWidth',4,'Color','c'); title 'Object';
end
