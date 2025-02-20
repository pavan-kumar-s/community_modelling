clear
changeCobraSolverParams('LP', 'feasTol', 1e-9);
%% Obtaining the file names
p = dir('./models/HDC1_gapseq/');
p = {p(3:end).name}';
p = p(contains(p,'xml'));

%% Obtaining the model names
taxaNames = readtable('./models/HDC1_gapseq/taxa_names.xlsx');
modelNames = {};
for m=1:numel(p)
    temp = taxaNames.Taxa(ismember(taxaNames.FileName,p{m}));
    modelNames =[modelNames;temp];
end

%% Loading all the models
models = cell(numel(p),1);
for m=1:numel(p)
    model = readCbModel(['./models/HDC1_gapseq/',p{m}]);
    models = [models;model];
end

%% getting the biomass ids for all the models
bio_id = [];
for m = 1:numel(models)
    model = models{m};
    bio_id =[bio_id;find(model.c)];
end
%% getting the weights of microbes
wts = readtable('./HDC1_weight.xlsx');


[X_k, V] = comModel(models, modelNames, X_k_exp, bio_id, X0, mu, EX_)