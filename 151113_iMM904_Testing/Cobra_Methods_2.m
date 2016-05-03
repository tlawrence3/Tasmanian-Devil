function Cobra_Methods(model, description, expressionData)
  warning('off');
  initCobraToolbox;
  %import the models and data
  model_string = sprintf('data/models/%s.mat',model);
  expressionData_string = sprintf('data/geneRules/%s.mat',description);
  load(model_string);
  load(expressionData_string);
  %perform pFBA
  fluxes_ub_pFBA_max = optimizeCbModel(model,'max','one');
  %perform GIMME
  [model_GIMME,Rxns_GIMME] = createTissueSpecificModel(model,expressionData,1,1,[],'GIMME',[find(model.c) 0.9],1);
  model_GIMME_taxicab = optimizeCbModel(model_GIMME,'max','one');
  %perform iMAT
  [model_iMAT,Rxns_iMAT] = createTissueSpecificModel(model,expressionData,1,1,[],'Shlomi',[],1);
  model_iMAT_taxicab = optimizeCbModel(model_iMAT,'max','one');
  save_string = sprintf('data/COBRAResults/%s_%s.mat',model,description);
  save(save_string);
