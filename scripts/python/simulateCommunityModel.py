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


def add_biomass_weight_and_update_constraints(merged_model, modelNames):
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
    # constraint = merged_model.problem.Constraint(sum([merged_model.variables[modelName+'_biomass_weight'] for modelName in modelNames]), lb=0, ub=X0, name='total_biomass_weight')
    merged_model.add_cons_vars(constraint)
    merged_model.solver.update()
    return merged_model

def get_exchange_rxn_met_mapping(model):
    exchange_reactions = [i.id for i in model.exchanges] # all the exchange reactions in the model
    exchange_mets = [list(i.metabolites.keys())[0].id for i in model.exchanges] # all the exchange metabolites in the model
    return pd.DataFrame({'Exchange_rxns':exchange_reactions, 'Exchange_mets':exchange_mets})

def add_metabolomics_constraints(merged_model, exchangeRate, met2rxn, beta=0):
    metabolomics_constraints = []
    for met in list(exchangeRate['metabolite ID']):
        exc_rxns = list(met2rxn[met2rxn.Exchange_mets.apply(lambda x: met in x)]['Exchange_rxns'])
        consumption_rate = exchangeRate[exchangeRate['metabolite ID']==met]['Consumed'].values[0]
        sd = exchangeRate[exchangeRate['metabolite ID']==met]['SD'].values[0]
        # constraint = merged_model.problem.Constraint(sum([merged_model.reactions.get_by_id(r).flux_expression for r in exc_rxns])+consumption_rate, lb=-beta, ub=beta,name=f"{met}_metabolomics_constraint")
        # constraint = merged_model.problem.Constraint(sum([merged_model.reactions.get_by_id(r).flux_expression for r in exc_rxns]), lb=-consumption_rate, ub=0,name=f"{met}_metabolomics_constraint")
        constraint = merged_model.problem.Constraint(sum([merged_model.reactions.get_by_id(r).flux_expression for r in exc_rxns]), lb=-consumption_rate-sd, ub=-consumption_rate+sd,name=f"{met}_metabolomics_constraint")
        metabolomics_constraints.append(constraint)
    merged_model.add_cons_vars(metabolomics_constraints)
    merged_model.solver.update()
    return merged_model

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

def set_model_objective(merged_model, modelNames):
    objective = merged_model.problem.Objective(sum([merged_model.variables['z_'+modelName] for modelName in modelNames]), direction='min')
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

def simulateCommunityModel(modelDetails, exchangeRate, mu, solver='gurobi',alpha=0, beta=0):
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
    merged_model = add_biomass_weight_and_update_constraints(merged_model, modelNames)
    merged_model = add_biomass_flux_constraints(growthRate, merged_model, biomass_reactions, modelNames)
    
    merged_model = add_total_biomass_constraint(X0, merged_model, modelNames, alpha)
    met2rxn = get_exchange_rxn_met_mapping(merged_model)
    merged_model = add_metabolomics_constraints(merged_model, exchangeRate, met2rxn, beta)
    merged_model = add_absolute_values_constraints(merged_model, modelNames, exp_weights)
    merged_model = set_model_objective(merged_model, modelNames)
    solution = merged_model.optimize(objective_sense=None)
    df_biomass = compare_actual_vs_predicted_biomass(merged_model,exp_weights, modelNames)
    print('Comparing the actual and predicted biomass weights')
    print(df_biomass)
    df_metabolomics = compare_actual_vs_predicted_metabolomics_consumption(merged_model,exchangeRate, met2rxn)
    print('Comparing the actual and predicted metabolomics consumption rates by the community')
    print(df_metabolomics)
    return solution, merged_model, df_biomass, df_metabolomics

    