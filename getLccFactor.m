% Copyright Jonathan L. Bradshaw 2018

function [lccFactor_constr_treat, lccFactor_om_treat, lccFactor_constr_convey_pump, lccFactor_constr_convey_pipe, lccFactor_om_convey, lccFactor_energy_convey] = getLccFactor(assessPeriod)

%% assumed parameters needed for computeLCCFactor
discEff = 0.03; % real discount rate (low)
lifetime = [50,20];
fracComponent_treat = [0.5,0.5];
fracComponent_pump = [0.5,0.5];
fracComponent_pipe = [1.0,0];
con2capFactor = 1.69;
isLoan = 1;
loanRate = 0.055;
payback = 25;
costEsc = 0.03; % cost escalation or inflation rate.
discNom = (1+discEff)*(1+costEsc)-1;

[lccFactor_constr_treat, lccFactor_om_treat] = computeLccFactor(discEff,lifetime,assessPeriod,fracComponent_treat,con2capFactor,isLoan,loanRate,payback,discNom);
[lccFactor_constr_convey_pump, lccFactor_om_convey] = computeLccFactor(discEff,lifetime,assessPeriod,fracComponent_pump,con2capFactor,isLoan,loanRate,payback,discNom);
[lccFactor_constr_convey_pipe, lccFactor_om_convey] = computeLccFactor(discEff,lifetime,assessPeriod,fracComponent_pipe,con2capFactor,isLoan,loanRate,payback,discNom);

lccFactor_constr_convey_pump = lccFactor_constr_convey_pump*con2capFactor;
lccFactor_constr_convey_pipe = lccFactor_constr_convey_pipe*con2capFactor;

lccFactor_energy_convey = lccFactor_om_convey; % lcc factor for energy is the same as for O&M

end