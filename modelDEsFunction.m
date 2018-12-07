function modelDEsFunction = modelDEsFunction(t,y)

% Here are the equations for the HV strain epidemic

global beta
global alpha
global mu

SI = y(1);
SN = y(2);
I = y(3);
R = y(4);

modelDEsFunction = zeros(4, 1);

modelDEsFunction(1) = -beta*(1-alpha)*SI*I;
modelDEsFunction(2) = -beta*SN*I;
modelDEsFunction(3) = beta*(1-alpha)*SI*I + beta*SN*I - mu*I;
modelDEsFunction(4) = mu*I;

end

