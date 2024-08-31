function [X_k, V] = comModel(models, modelNames, X_k_exp, bio_id, X0, mu, EX_)
% USAGE:
%   [X_k, V] = comModel(models, X_k_exp, X0, mu, EX_)
%
% INPUTS:
%  models:      Cell array of models
%  modelNames:  Cell array of model names (same order as models)
%  X_k_exp:     Array of gdw (gram dry weight) of the microbes obtained 
%               experimentally. The order of the X_k_exp should match with the
%               respective microbes in 'models'
%  bio_id:      IDs of biomass reactions in the respective models
%  X0:          Non-zero total biomass
%  mu:          Community growth rate (or) the dilution rate
%  EX_:         Community export and import rates as matlab struct
%                .mets: metabolite IDs which are known to be consumed
%                .rxns:   Exchange reactions corresponding to EX_.mets
%                .rate: Consumption rate (S_com - S_gf) as positive values
%
% OUTPUTS:
%  X_k:      Biomass in gram dry weight for all the microbes in same order
%            as models
%  V:        Aggregated flux through all the reactions in the community


% The decision variables are in the below given order
    % X_k (The biomass in gram dry weight), 
    % Z (Decision variables to replace the absolute value in objective function),
    % V (Aggregated reaction flux)

 
% Total no. of microbes in the community
n_models = numel(models);

% Total no. of reactions and metabolites in the community
n_rxns = 0; n_mets = 0; all_rxns = {};
for m =1:n_models
    temp = models{m};
    n_rxns = n_rxns+numel(temp.rxns);
    n_mets = n_mets+numel(temp.mets);
    all_rxns = [all_rxns;temp.rxns];
end

%% Setting the objective function
% the objective function is the summation of all z variables that will be
% equal to the sum of absolute difference between the X_k and X_k_exp
f = [zeros(n_models,1);ones(n_models,1);zeros(n_rxns,1)];

%% Setting the stoichiometric constraints (equality constraints)
model = models{1};
Aeq1 = model.S;
for m = 2:n_models
    model = models{m};
    Aeq1 = blkdiag(Aeq1,model.S);
end
Aeq1 = [sparse(n_mets,2*n_models),Aeq1]; % LHS of the equality constraint
beq1 = zeros(n_mets,1); % RHS of the equality constraint
csenseeq1 = repmat('E',n_mets,1); % 'E' refers to the equality constraints

%% Setting the biomass constraints (equality constraints)
prev_id = 0;
Aeq2 = sparse(n_models,n_rxns);
for m = 1:n_models
    Aeq2(m, prev_id+bio_id(m)) = -1;
    model = models{m};
    prev_id = prev_id+numel(model.rxns);
end
Aeq2 = [mu*speye(n_models,n_models),sparse(n_models,n_models),Aeq2]; % LHS of the equality constraint
beq2 = zeros(n_models,1); % RHS of the equality constraint
csenseeq2 = repmat('E',n_models,1); % 'E' refers to the equality constraints

%% Setting the box constraints: LB (inequality constraints)
Aineq1 = sparse(n_rxns,n_models);
prev_id = 0;
for m = 1:n_models
    model = models{m};
    temp = prev_id+numel(model.rxns);
    Aineq1(prev_id+1:temp,m) = model.lb;
    prev_id = temp;
end
Aineq1 = [Aineq1, sparse(n_rxns,n_models), -1*speye(n_rxns,n_rxns)]; % LHS of the inequality constraint
bineq1 = zeros(n_rxns,1); % RHS of the inequality constraint
csenseineq1 = repmat('L',n_rxns,1); % 'L' refers to the lesser than

%% Setting the box constraints: UB (inequality constraints)
Aineq2 = sparse(n_rxns,n_models);
prev_id = 0;
for m = 1:n_models
    model = models{m};
    temp = prev_id+numel(model.rxns);
    Aineq2(prev_id+1:temp,m) = -1*model.ub;
    prev_id = temp;
end
Aineq2 = [Aineq2, sparse(n_rxns,n_models), speye(n_rxns,n_rxns)]; % LHS of the inequality constraint
bineq2 = zeros(n_rxns,1); % RHS of the inequality constraint
csenseineq2 = repmat('L',n_rxns,1); % 'L' refers to the lesser than

%% Setting the constraints on the uptake reactions (Exchange reactions) (equality constraints)
n_exchange = numel(EX_.rxns); % number of exchange metabolites
Aeq3 = sparse(n_exchange,n_rxns);
for e = 1:n_exchange
    ex_ids = ismember(all_rxns,EX_.rxns(e));
    Aeq3(e,ex_ids) = 1;
end
Aeq3 = [sparse(n_exchange,2*n_models),Aeq3]; % LHS of the equality constraint
beq3 = -1*EX_.rate; % RHS of the equality constraint
csenseeq3 = repmat('E',n_exchange,1); % 'E' refers to the equality constraints

%% Setting the constraint on the total biomass weight (equality constraint)
Aeq4 = [ones(1,n_models),sparse(1,n_models+n_rxns)]; % LHS of the equality constraint
beq4 = X0; % RHS of the equality constraint
csenseeq4 = repmat('E',1,1); % 'E' refers to the equality constraints

%% Setting the constraints on decision variables "Z" and the absolute difference between the X_k and X_k_exp
%  (Inequality constraint)
Aineq3 = [-1*speye(n_models,n_models),-1*speye(n_models,n_models),sparse(n_models,n_rxns)]; % LHS of the inequality constraint
bineq3 = -1*X_k_exp(:); % RHS of the inequality constraint
csenseineq3 = repmat('L',n_models,1); % 'L' refers to the lesser than

%% Setting the constraints on decision variables "Z" and the absolute difference between the X_k and X_k_exp
%  (Inequality constraint)
Aineq4 = [speye(n_models,n_models),-1*speye(n_models,n_models),sparse(n_models,n_rxns)]; % LHS of the inequality constraint
bineq4 = X_k_exp(:) ; % RHS of the inequality constraint
csenseineq4 = repmat('L',n_models,1); % 'L' refers to the lesser than

%% bounds on the decision variables
lb = [sparse(2*n_models,1);-Inf(n_rxns,1)];
ub = Inf((2*n_models)+n_rxns,1);

%% Setting up LP problem

LPproblem.A=[Aeq1;Aeq2;Aeq3;Aeq4;Aineq1;Aineq2;Aineq3;Aineq4];
LPproblem.b=[beq1;beq2;beq3;beq4;bineq1;bineq2;bineq3;bineq4];
LPproblem.lb=lb;
LPproblem.ub=ub;
LPproblem.c=f;
LPproblem.osense=1;%minimise
LPproblem.csense = [csenseeq1;csenseeq2;csenseeq3;csenseeq4;csenseineq1;csenseineq2;csenseineq3;csenseineq4];
solution = solveCobraLP(LPproblem);
if solution.stat~=1
    fprintf('%s%s\n',num2str(solution.stat),' = solution.stat')
    fprintf('%s%s\n',num2str(solution.origStat),' = solution.origStat')
    warning('LP solution may not be optimal')
end

x=solution.full;

X_k = x(1:n_models);
V = x((2*n_models)+1:end);
tbl = table('Size',[n_rxns,3],'VariableNames',{'Model','Reactions','Aggregate_Flux'},'VariableTypes',...
    {'string','string','double'});
Model = {}; Reactions = {}; 
for m=1:n_models
    model = models{m};
    Model = [Model;repmat(modelNames(m),numel(model.rxns),1)];
    Reactions = [Reactions;model.rxns];
end
tbl.Model = Model; tbl.Reactions = Reactions; tbl.Aggregate_Flux = V;
V = tbl;
