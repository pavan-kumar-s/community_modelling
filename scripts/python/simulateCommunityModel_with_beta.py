from cobra.io import load_model, load_json_model, load_matlab_model, read_sbml_model
from cobra.util.solver import linear_reaction_coefficients
import cobra
import pandas as pd

def load_and_process_models(path2models, modelNames):
    models = []
    biomass_reactions = []
    for path in path2models:
        if path.endswith('.mat'):
            model = load_matlab_model(path)
        elif path.endswith('.json'):
            model = load_json_model(path)
        elif path.endswith('.xml'):
            model = read_sbml_model(path)
        model.solver.problem.Params.FeasibilityTol = 1e-9
        model.solver.problem.Params.OptimalityTol = 1e-9
        # Rename reactions
        for reaction in model.reactions:
            reaction.id = f"{modelNames[path2models.index(path)]}_{reaction.id}"

        # Rename metabolites
        for metabolite in model.metabolites:
            metabolite.id = f"{modelNames[path2models.index(path)]}_{metabolite.id}"

        # Repair the model to update internal indexes
        model.repair()
        biomass_reactions.append([reaction for reaction in model.reactions if reaction.objective_coefficient != 0][0].id)
        models.append(model)
    return models, biomass_reactions


def merge_models(models):
    merged_model = models[0].copy()
    for model in models[1:]:
        merged_model.merge(model)
    return merged_model


def add_biomass_weight_variables_and_updating_lb_ub_constraints(merged_model, modelNames):
    # defining the variables for the biomass weights
    name2expression = {m:merged_model.problem.Variable(m+'_biomass_weight', lb=0, ub=1000, type="continuous") for m in modelNames}
    merged_model.add_cons_vars(name2expression.values())
    merged_model.solver.update()
    # setting the bounds for the reactions based on the biomass weights and default bounds
    constraints = []
    for reaction in merged_model.reactions:
        bio_id = [i for i,j in enumerate(modelNames) if j in reaction.id][0] # getting the index of the model name 
        constraint_ub = merged_model.problem.Constraint(reaction.flux_expression-(reaction.upper_bound*name2expression[modelNames[bio_id]]), ub=0, name=f"{reaction.id}_bio_weight_ub")
        constraint_lb = merged_model.problem.Constraint(reaction.flux_expression-(reaction.lower_bound*name2expression[modelNames[bio_id]]), lb=0, name=f"{reaction.id}_bio_weight_lb")
        constraints.append(constraint_ub)
        constraints.append(constraint_lb)
    merged_model.add_cons_vars(constraints)
    merged_model.solver.update()
    return merged_model


def add_biomass_flux_constraints(growthRate, merged_model, biomass_reactions, modelNames):
    # setting the biomass flux constraints
    biomass_constraints = []
    for reaction in biomass_reactions:
        modelName = modelNames[[i for i,j in enumerate(modelNames) if j in reaction][0]]
        biomass_wt = modelName+'_biomass_weight'
        constraint = merged_model.problem.Constraint(merged_model.reactions.get_by_id(reaction).flux_expression-growthRate*merged_model.variables[biomass_wt], lb=0, ub=0,name=f"{modelName}_flux_bio_weight")
        biomass_constraints.append(constraint)
    merged_model.add_cons_vars(biomass_constraints)
    merged_model.solver.update()
    return merged_model

def add_total_biomass_constraint(X0, merged_model, modelNames, alpha=0):
    constraint = merged_model.problem.Constraint(sum([merged_model.variables[modelName+'_biomass_weight'] for modelName in modelNames])-X0, lb=-alpha, ub=alpha, name='total_biomass_weight')
    merged_model.add_cons_vars(constraint)
    merged_model.solver.update()
    return merged_model

def get_exchange_rxn_met_mapping(model):
    exchange_reactions = [i.id for i in model.exchanges] # all the exchange reactions in the model
    exchange_mets = [list(i.metabolites.keys())[0].id for i in model.exchanges] # all the exchange metabolites in the model
    return pd.DataFrame({'Exchange_rxns':exchange_reactions, 'Exchange_mets':exchange_mets})


def add_metabolomics_constraints(merged_model, exchangeRate, met2rxn, growthRate):
    metabolomics_constraints = []
    consumption_rates = []
    for met in list(exchangeRate['metabolite ID']):
        com_conc = exchangeRate[exchangeRate['metabolite ID']==met]['Microbiome_conc'].values[0]
        gf_conc = exchangeRate[exchangeRate['metabolite ID']==met]['GF_conc'].values[0]

        # creating a beta variable speicific to the metabolite
        beta = merged_model.problem.Variable(f"beta_{met}", lb=0, ub=10, type="continuous")
        # Adding the variables representing the concentration of the metabolite in the community and GF mice
        s_com = merged_model.problem.Variable(f"com_conc_{met}", lb=0, ub=1000, type="continuous")
        s_gf = merged_model.problem.Variable(f"gf_conc_{met}", lb=0, ub=1000, type="continuous")
        s_gf_lb_constraint = merged_model.problem.Constraint(s_gf+beta, lb=gf_conc, name=f"{met}_gf_lb")
        s_gf_ub_constraint = merged_model.problem.Constraint(s_gf-beta, ub=gf_conc, name=f"{met}_gf_ub")
        s_com_lb_constraint = merged_model.problem.Constraint(s_com+beta, lb=com_conc, name=f"{met}_com_lb")
        s_com_ub_constraint = merged_model.problem.Constraint(s_com-beta, ub=com_conc, name=f"{met}_com_ub")

        exc_rxns = list(met2rxn[met2rxn.Exchange_mets.apply(lambda x: met in x)]['Exchange_rxns']) # getting the list of exchange reactions for the metabolite
        consumption_rate =(s_gf - s_com)*growthRate # constraint referring to the consumption rate of the metabolite
        met_constraint = merged_model.problem.Constraint(sum([merged_model.reactions.get_by_id(r).flux_expression for r in exc_rxns])+consumption_rate, lb=0, ub=0,name=f"{met}_metabolomics_constraint")
        
        metabolomics_constraints.extend([s_gf_lb_constraint, s_gf_ub_constraint, s_com_lb_constraint, s_com_ub_constraint, met_constraint])


        consumption_rates.append((gf_conc-com_conc)*growthRate) # storing the mean consumption rate of the metabolite for comparison
    exchangeRate['Consumed'] = consumption_rates # updating the exchangeRate dataframe with the mean consumption rate
    merged_model.add_cons_vars(metabolomics_constraints)
    merged_model.solver.update()
    return merged_model,exchangeRate

def add_absolute_values_constraints(merged_model, modelNames, exp_weights):
    constraints = []
    for i,weight in enumerate(exp_weights):
        z_variable = merged_model.problem.Variable(f"z_{modelNames[i]}", lb=0, ub=1000, type="continuous")
        constraint_lb = merged_model.problem.Constraint(merged_model.variables[modelNames[i]+'_biomass_weight']-weight+z_variable, lb=0, name=f"{modelNames[i]}_exp_weight_lb")
        constraint_ub = merged_model.problem.Constraint(merged_model.variables[modelNames[i]+'_biomass_weight']-weight-z_variable, ub=0, name=f"{modelNames[i]}_exp_weight_ub")
        constraints.append(constraint_lb)
        constraints.append(constraint_ub)
    merged_model.add_cons_vars(constraints)
    merged_model.solver.update()
    return merged_model

def set_model_objective(merged_model, modelNames, exchangeRate,wt=0.5):
    objective = merged_model.problem.Objective(wt*sum([merged_model.variables['z_'+modelName] for modelName in modelNames])+(1-wt)*sum([merged_model.variables['beta_'+met] for met in list(exchangeRate['metabolite ID'])]), direction='min')
    merged_model.objective = objective
    return merged_model

def compare_actual_vs_predicted_biomass(merged_model,exp_weights, modelNames):
    predicted_biomass = [merged_model.variables[modelName+'_biomass_weight'].primal for modelName in modelNames]
    df = pd.DataFrame({'Model':modelNames, 'Predicted':predicted_biomass, 'Actual':exp_weights})
    return df

def compare_actual_vs_predicted_metabolomics_consumption(merged_model,exchangeRate, met2rxn):
    predicted_consumption, actual_consumption = [],[]
    for met in list(exchangeRate['metabolite ID']):
        exc_rxns = list(met2rxn[met2rxn.Exchange_mets.apply(lambda x: met in x)]['Exchange_rxns'])
        actual_consumption.append(exchangeRate[exchangeRate['metabolite ID']==met]['Consumed'].values[0])
        predicted_consumption.append(-sum([merged_model.reactions.get_by_id(r).flux for r in exc_rxns]))
    df = pd.DataFrame({'Metabolite':list(exchangeRate['metabolite ID']), 'Predicted':predicted_consumption, 'Actual':actual_consumption})
    return df

def compare_actual_vs_predicted_metabolomics_production(merged_model, productionRate, met2rxn, growthRate):
    predicted_production, actual_production = [],[]
    for met in list(productionRate['metabolite ID']):
        exc_rxns = list(met2rxn[met2rxn.Exchange_mets.apply(lambda x: met in x)]['Exchange_rxns'])
        gf_conc = productionRate[productionRate['metabolite ID']==met]['GF_conc'].values[0]
        com_conc = productionRate[productionRate['metabolite ID']==met]['Microbiome_conc'].values[0]
        produced_rate = (com_conc-gf_conc)*growthRate
        actual_production.append(produced_rate)
        predicted_production.append(sum([merged_model.reactions.get_by_id(r).flux for r in exc_rxns]))
    df = pd.DataFrame({'Metabolite':list(productionRate['metabolite ID']), 'Predicted':predicted_production, 'Actual':actual_production})
    return df

def simulateCommunityModel_with_beta(modelDetails, exchangeRate, mu, solver='gurobi', alpha=0,productionRate=None, tradeoff=0.5):
    cobra_config = cobra.Configuration()
    cobra_config.solver = solver
    
    modelPaths = list(modelDetails.Path)
    modelNames = list(modelDetails.Name)
    exp_weights = list(modelDetails.weight)
    X0 = sum(exp_weights)
    growthRate = mu

    ## running the simulation
    models, biomass_reactions =load_and_process_models(modelPaths, modelNames)
    merged_model = merge_models(models)
    print("Models are merged")
    merged_model = add_biomass_weight_variables_and_updating_lb_ub_constraints(merged_model, modelNames)
    print("Lower and upper bounds are updated based on the biomass weights")
    merged_model = add_biomass_flux_constraints(growthRate, merged_model, biomass_reactions, modelNames)
    print("Biomass flux constraints are added based on the growth rate")
    merged_model = add_total_biomass_constraint(X0, merged_model, modelNames, alpha)
    print("Total biomass constraint is added")
    met2rxn = get_exchange_rxn_met_mapping(merged_model)
    merged_model, exchangeRate = add_metabolomics_constraints(merged_model, exchangeRate, met2rxn, growthRate)
    print("Metabolomics constraints are added")
    merged_model = add_absolute_values_constraints(merged_model, modelNames, exp_weights)
    print("Absolute value constraints are added that are needed for the objective function")
    merged_model = set_model_objective(merged_model, modelNames, exchangeRate, tradeoff)
    print("Optimizing the model")
    solution = merged_model.optimize(objective_sense=None)
    print("Optimization is done")
    df_biomass = compare_actual_vs_predicted_biomass(merged_model,exp_weights, modelNames)
    print('Comparing the actual and predicted biomass weights')
    print(df_biomass)
    print('Comparing the actual and predicted metabolomics consumption rates by the community')
    df_metabolomics_consumed = compare_actual_vs_predicted_metabolomics_consumption(merged_model,exchangeRate, met2rxn)
    print(df_metabolomics_consumed)
    if productionRate is not None:
        print('Comparing the actual and predicted metabolomics production rates by the community')
        df_metabolomics_produced = compare_actual_vs_predicted_metabolomics_production(merged_model,productionRate, met2rxn, growthRate)
        print(df_metabolomics_produced)
        return solution, merged_model, df_biomass, df_metabolomics_consumed, df_metabolomics_produced
    else:
        return solution, merged_model, df_biomass, df_metabolomics_consumed

    