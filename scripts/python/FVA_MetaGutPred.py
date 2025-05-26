import pandas as pd
from tqdm import tqdm
def addBiomassWeightConstraints(merged_model, modelNames):
    # getting the predicted biomass weight for each model
    predicted_biomass = [merged_model.variables[modelName+'_biomass_weight'].primal for modelName in modelNames]

    # modifying the lb and ub of the biomass weight variables
    for modelName, biomass_weight in zip(modelNames, predicted_biomass):
        merged_model.variables[modelName+'_biomass_weight'].bounds = (biomass_weight, biomass_weight)
    merged_model.solver.update()

    return merged_model

def FVA_MetaGutPred(merged_model, reactions, modelNames):

    # Add biomass weight constraints
    merged_model = addBiomassWeightConstraints(merged_model, modelNames)

    df = pd.DataFrame()
    minFlux, maxFlux = [], []
    for reaction in tqdm(reactions, desc='FVA MetaGutPred'):
        # optimizing the model for min flux
        objective = merged_model.problem.Objective(merged_model.reactions.get_by_id(reaction).flux_expression, direction='min')
        merged_model.objective = objective
        solution = merged_model.optimize()
        minFlux.append(solution.objective_value)
        
        # optimizing the model for max flux
        objective = merged_model.problem.Objective(merged_model.reactions.get_by_id(reaction).flux_expression, direction='max')
        merged_model.objective = objective
        solution = merged_model.optimize()
        maxFlux.append(solution.objective_value)
    # creating a dataframe with the results
    df['reactions'] = reactions
    df['minFlux'] = minFlux
    df['maxFlux'] = maxFlux
    return df