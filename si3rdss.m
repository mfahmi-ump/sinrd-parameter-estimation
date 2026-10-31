function ss = si3rdss(k,data)
Ptotal = 3.153*10^7;%malaysia
time = data(:,1);
I1obs = data(:,2);
I2obs = data(:,3);
I3obs = data(:,4);
Robs = data(:,5);
Dobs = data(:,6);
Sobs = Ptotal - I1obs - I2obs - I3obs - Robs - Dobs;

s0 = Sobs(1);
i10 = I1obs(1);
i20 = I2obs(1);
i30 = I3obs(1);
r0 = Robs(1);
d0 = Dobs(1);

beta1_start = k(1);
beta2_start = k(2);
beta3_start = k(3);
r1_start = k(4);
r2_start = k(5);
r3_start = k(6);
d1_start = k(7);
d2_start = k(8);
d3_start = k(9);
p1_start = k(10);
p2_start = k(11);
p3_start = k(12);
gamma_start = k(13);

[Smodel, I1model, I2model, I3model, Rmodel, Dmodel] = si3rdsrk4(time,s0,i10,i20,i30,r0,d0,Ptotal,beta1_start,beta2_start,beta3_start,r1_start,r2_start,r3_start,d1_start,d2_start,d3_start,p1_start,p2_start,p3_start,gamma_start);
 
% Serr = (sum(Sobs - Smodel.').^2);
% I1err = (sum(I1obs - I1model.').^2);
% I2err = (sum(I2obs - I2model.').^2);
% I3err = (sum(I3obs - I3model.').^2);
% Rerr = (sum(Robs - Rmodel.').^2);
% Derr = (sum(Dobs - Dmodel.').^2);

SRerr = (sum(Sobs - Smodel.').^2)/max(Smodel);
I1Rerr = (sum(I1obs - I1model.').^2)/max(I1model);
I2Rerr = (sum(I2obs - I2model.').^2)/max(I2model);
I3Rerr = (sum(I3obs - I3model.').^2)/max(I3model);
RRerr = (sum(Robs - Rmodel.').^2)/max(Rmodel);
DRerr = (sum(Dobs - Dmodel.').^2)/max(Dmodel);

SGerr = ((Sobs(end) - Smodel(end)).^2)/max(Smodel);
I1Gerr = ((I1obs(end) - I1model(end)).^2)/max(I1model);
I2Gerr = ((I2obs(end) - I2model(end)).^2)/max(I2model);
I3Gerr = ((I3obs(end) - I3model(end)).^2)/max(I3model);
RGerr = ((Robs(end) - Rmodel(end)).^2)/max(Rmodel);
DGerr = ((Dobs(end) - Dmodel(end)).^2)/max(Dmodel);

% Serr = ((Sobs(end) - Smodel(end)).^2);
% Ierr = ((Iobs(end) - Imodel(end)).^2);
% Rerr = ((Robs(end) - Rmodel(end)).^2);
% Derr = ((Dobs(end) - Dmodel(end)).^2);

ss = (SRerr + I1Rerr + I2Rerr + I3Rerr + DRerr + RRerr) + (SGerr + I1Gerr*1.5 + I2Gerr + I3Gerr + DGerr*1.5 + RGerr*1.5);
end
