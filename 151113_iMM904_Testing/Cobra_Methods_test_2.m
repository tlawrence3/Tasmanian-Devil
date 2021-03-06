function Cobra_Methods(model, description, lb_biomass, biomass_reaction)
  if nargin < 4
    biomass_reaction = '';
  end
  if nargin < 3
    lb_biomass = 0;
  end 
  warning('off');
  initCobraToolbox;
  %import the models and data
  model_string = sprintf('data/models/%s',model);
  expressionData_string = sprintf('data/geneRules/%s.mat',description);
  m = load(model_string);
  data = load(expressionData_string);
  %perform pFBA
  fluxes_pFBA_max = optimizeCbModel(m,'max','one');
  %perform GIMME
  [model_GIMME,Rxns_GIMME] = createTissueSpecificModel(m,data,1,1,[],'GIMME',[find(m.c) 0.9],1);
  model_GIMME_taxicab = optimizeCbModel(model_GIMME,'max','one');
  %perform iMAT
  [model_iMAT,Rxns_iMAT] = createTissueSpecificModel(m,data,1,1,[],'Shlomi',[],1);
  model_iMAT_taxicab = optimizeCbModel(model_iMAT,'max','one');
  save_string = sprintf('data/COBRAResults/Results_%s_%s',description,model);
  if lb_biomass ~= 0
    m_lb = m;
    m_lb.lb(find(strcmp(m_lb.rxns,biomass_reaction),1)) = lb_biomass;
    %perform pFBA
    fluxes_lb_pFBA_max = optimizeCbModel(m_lb,'max','one');
    %perform GIMME
    [model_lb_GIMME,Rxns_lb_GIMME] = createTissueSpecificModel(m_lb,data,1,1,[],'GIMME',[find(m_lb.c) 0.9],1);
    model_lb_GIMME_taxicab = optimizeCbModel(model_lb_GIMME,'max','one');
    %perform iMAT
    [model_lb_iMAT,Rxns_lb_iMAT] = createTissueSpecificModel(m_lb,data,1,1,[],'Shlomi',[],1);
    model_lb_iMAT_taxicab = optimizeCbModel(model_lb_iMAT,'max','one');

    m_ub = m;
    m_ub.ub(find(strcmp(m_ub.rxns,biomass_reaction),1)) = lb_biomass;
    %perform pFBA
    fluxes_ub_pFBA_max = optimizeCbModel(m_ub,'max','one');
    %perform GIMME
    [model_ub_GIMME,Rxns_ub_GIMME] = createTissueSpecificModel(m_ub,data,1,1,[],'GIMME',[find(m_ub.c) 0.9],1);
    model_ub_GIMME_taxicab = optimizeCbModel(model_ub_GIMME,'max','one');
    %perform iMAT
    [model_ub_iMAT,Rxns_ub_iMAT] = createTissueSpecificModel(m_ub,data,1,1,[],'Shlomi',[],1);
    model_ub_iMAT_taxicab = optimizeCbModel(model_ub_iMAT,'max','one');

    save_string = sprintf('data/COBRAResults/Results_%s_%s_lb.mat',description,model(1:length(model)-4));
  end
  clear model description model_string expressionData_string lb_biomass biomass_reaction;
  save(save_string);
