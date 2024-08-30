clear

% Names of the model
modelNames={'Abiotrophia_defectiva_ATCC_49176';'Acinetobacter_baumannii_ABNIH3';...
    'Escherichia_coli_DEC7B';'Herbaspirillum_huttiense_NFYY_53159';...
    'Lactobacillus_crispatus_214_1';'Lactococcus_garvieae_I113';...
    'Mobiluncus_curtisii_ATCC_43063';'Salmonella_enterica_subsp_enterica_serovar_Newport_str_CVM_19536';...
    'Staphylococcus_aureus_M0946';'Ureaplasma_parvum_serovar_3_str_ATCC_700970'};

% loading all the models
models = {};
for m=1:numel(modelNames)
    load(['./models/',modelNames{m}])
    models = [models;model];
end

% giving some random numbers to the experiemntally observed biomass weights in gdw
X_k_exp = rand(numel(modelNames),1)*100;

% getting the biomass ids for all the models
bio_id = [];
for m = 1:numel(modelNames)
    model = models{m};
    bio_id =[bio_id;find(model.c)];
end

% total sum of the biomass
X0 = sum(X_k_exp);

% growth rate
mu = 0.7;

% defining the community export rates. It has to be positive values (ei_c - ui_c)

% here we are randomly choosing 50 reactions as exchange reactions in first model
model = models{1};
exRxns = model.rxns(startsWith(model.rxns,'EX_'));
exRxns = setdiff(exRxns,{'EX_biomass(e)'});
randIDS = sort(randsample(numel(exRxns),20));
exRxns = exRxns(randIDS);
exMets = findMetsFromRxns(model,exRxns);

% giving random numbers as consumption rate
rate = rand(numel(exRxns),1)*10;

EX_.mets = exMets;
EX_.rxns = exRxns;
EX_.rate = rate;

[X_k, V] = comModel(models, modelNames, X_k_exp, bio_id, X0, mu, EX_);
