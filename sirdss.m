function ss = sirdss(k,data)
Ptotal = 3.153*10^7;%malaysia
time = data(:,1);
Iobs = data(:,2);
Robs = data(:,3);
Dobs = data(:,4);

Sobs = Ptotal - Iobs - Robs - Dobs;

s0 = Sobs(1);
i0 = Iobs(1);
r0 = Robs(1);
d0 = Dobs(1);

beta_start = abs(k(1));
r_start = abs(k(2));
d_start = abs(k(3));
gamma_start = abs(k(4));

[Smodel, Imodel, Rmodel, Dmodel] = sirdsrk4(time,s0,i0,r0,d0,Ptotal,beta_start,r_start,d_start,gamma_start,0.05);
 
% Serr = (sum(Sobs - Smodel.').^2);
% I1err = (sum(I1obs - I1model.').^2);
% I2err = (sum(I2obs - I2model.').^2);
% I3err = (sum(I3obs - I3model.').^2);
% Rerr = (sum(Robs - Rmodel.').^2);
% Derr = (sum(Dobs - Dmodel.').^2);

SRerr = (sum(Sobs - Smodel.').^2)/max(Smodel);
IRerr = (sum(Iobs - Imodel.').^2)/max(Imodel);
RRerr = (sum(Robs - Rmodel.').^2)/max(Rmodel);
DRerr = (sum(Dobs - Dmodel.').^2)/max(Dmodel);

SGerr = ((Sobs(end) - Smodel(end)).^2)/max(Sobs);
IGerr = ((Iobs(end) - Imodel(end)).^2)/max(Iobs);
RGerr = ((Robs(end) - Rmodel(end)).^2)/max(Robs);
DGerr = ((Dobs(end) - Dmodel(end)).^2)/max(Dobs);

% Serr = ((Sobs(end) - Smodel(end)).^2);
% Ierr = ((Iobs(end) - Imodel(end)).^2);
% Rerr = ((Robs(end) - Rmodel(end)).^2);
% Derr = ((Dobs(end) - Dmodel(end)).^2);

ss = (SRerr + IRerr + DRerr + RRerr) + (SGerr + IGerr + DGerr + RGerr);
end
