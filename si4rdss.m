function ss = si4rdss(k,data)
Ptotal = 3.153*10^7;%malaysia
time = data(:,1);
I1obs = data(:,2);
I2obs = data(:,3);
I3obs = data(:,4);
I4obs = data(:,5);
Robs = data(:,6);
Dobs = data(:,7);
Sobs = data(:,8);
ESobs = data(:,9);

s0 = Sobs(1);
i10 = I1obs(1);
i20 = I2obs(1);
i30 = I3obs(1);
i40 = I4obs(1);
r0 = Robs(1);
d0 = Dobs(1);

p1_start = k(1);
p2_start = k(2);
p3_start = k(3);
p4_start = k(4);
beta1_start = k(5);
beta2_start = k(6);
beta3_start = k(7);
beta4_start = k(8);
r1_start = k(9);
r2_start = k(10);
r3_start = k(11);
r4_start = k(12);
d1_start = k(13);
d2_start = k(14);
d3_start = k(15);
d4_start = k(16);
gamma_start = k(17);
p_vax_start = k(18);

[Smodel, I1model, I2model, I3model, I4model, Rmodel, Dmodel] = si4rdsrk4(time,s0,i10,i20,i30,i40,r0,d0,Ptotal,p1_start,p2_start,p3_start,p4_start,beta1_start,beta2_start,beta3_start,beta4_start,r1_start,r2_start,r3_start,r4_start,d1_start,d2_start,d3_start,d4_start,gamma_start,p_vax_start);

% Serr = (sum(Sobs - Smodel.').^2);
% I1err = (sum(I1obs - I1model.').^2);
% I2err = (sum(I2obs - I2model.').^2);
% I3err = (sum(I3obs - I3model.').^2);
% Rerr = (sum(Robs - Rmodel.').^2);
% Derr = (sum(Dobs - Dmodel.').^2);

% SRerr2 = (sum(Sobs(lengthofcurrparam) - Smodel(lengthofcurrparam).').^2)/max(Smodel(lengthofcurrparam));
% I1Rerr2 = (sum(I1obs(lengthofcurrparam) - I1model(lengthofcurrparam).').^2)/max(I1model(lengthofcurrparam));
% I2Rerr2 = (sum(I2obs(lengthofcurrparam) - I2model(lengthofcurrparam).').^2)/max(I2model(lengthofcurrparam));
% I3Rerr2 = (sum(I3obs(lengthofcurrparam) - I3model(lengthofcurrparam).').^2)/max(I3model(lengthofcurrparam));
% I4Rerr2 = (sum(I4obs(lengthofcurrparam) - I4model(lengthofcurrparam).').^2)/max(I4model(lengthofcurrparam));
% RRerr2 = (sum(Robs(lengthofcurrparam) - Rmodel(lengthofcurrparam).').^2)/max(Rmodel(lengthofcurrparam));
% DRerr2 = (sum(Dobs(lengthofcurrparam) - Dmodel(lengthofcurrparam).').^2)/max(Dmodel(lengthofcurrparam));

SRerr = (sum(Sobs - Smodel.').^2)/max(Smodel);
I1Rerr = (sum(I1obs - I1model.').^2)/max(I1model);
I2Rerr = (sum(I2obs - I2model.').^2)/max(I2model);
I3Rerr = (sum(I3obs - I3model.').^2)/max(I3model);
I4Rerr = (sum(I4obs - I4model.').^2)/max(I4model);
RRerr = (sum(Robs - Rmodel.').^2)/max(Rmodel);
DRerr = (sum(Dobs - Dmodel.').^2)/max(Dmodel);
VRerr = (mean(p_vax_start) - mean(ESobs./Sobs)).^2;


SGerr = ((Sobs(end) - Smodel(end)).^2)/mean(Smodel);
I1Gerr = ((I1obs(end) - I1model(end)).^2)/mean(I1model);
I2Gerr = ((I2obs(end) - I2model(end)).^2)/mean(I2model);
I3Gerr = ((I3obs(end) - I3model(end)).^2)/mean(I3model);
I4Gerr = ((I4obs(end) - I4model(end)).^2)/mean(I4model);
RGerr = ((Robs(end) - Rmodel(end)).^2)/mean(Rmodel);
DGerr = ((Dobs(end) - Dmodel(end)).^2)/mean(Dmodel);

% Serr = ((Sobs(end) - Smodel(end)).^2);
% Ierr = ((Iobs(end) - Imodel(end)).^2);
% Rerr = ((Robs(end) - Rmodel(end)).^2);
% Derr = ((Dobs(end) - Dmodel(end)).^2);
% 
ss = (SRerr + I1Rerr + I2Rerr + I3Rerr + I4Rerr + DRerr + RRerr) + (SGerr + I1Gerr*1.5 + I2Gerr + I3Gerr + I4Gerr + DGerr + RGerr + VRerr);
end
