%%MATLAB
function test = test_func(args)
warning('off');
initCobraToolbox;
#import the models and data
arg1 = args.arg1;
arg2 = args.arg2;
m = load('%s.mat' % arg1);
data = load('%s.mat' % arg2);
#perform pFBA
fluxes_ub_pFBA_max = optimizeCbModel(m,'max','one');
#perform GIMME
[model_GIMME,Rxns_GIMME] = createTissueSpecificModel(m,data,1,1,[],'GIMME',[find(m.c) 0.9],1);
Aerobic_fluxes_model_GIMME_taxicab = optimizeCbModel(model_GIMME,'max','one');
#perform iMAT
[model_iMAT,Rxns_iMAT] = createTissueSpecificModel(m,data,1,1,[],'Shlomi',[],1);
Aerobic_fluxes_model_iMAT_taxicab = optimizeCbModel(model_iMAT,'max','one')
save('160422_test.mat');
end
