%%
clear model data params options
datafile = readmatrix('sird_data.csv');
Ptotal = 3.153*10^7;%malaysia

Idata = datafile(:,3);
Rdata = datafile(:,2);
Ddata = datafile(:,4);
Sdata = datafile(:,5);
tdata = datafile(:,1);
t = 1:1:length(Idata);
% mu = 0.000021;
mu = 0;
%%
% slopedetection
dIdt = diff([eps; Idata(:)])./diff([eps; t(:)]);
sdIdt = smoothdata(dIdt,'gaussian',12);
d2Idt2 = diff([eps; sdIdt(:)])./diff([eps; t(:)]);
sd2Idt2 = smoothdata(d2Idt2,'gaussian',15);
[pks,locs] = findpeaks(abs(sd2Idt2));
tdays = datetime(2020,3,17) + caldays(1:length(t));
plot(tdays,Idata,'.');
xtickformat('dd-MM-yy')
xline([locs]);
startdatelist = [1,locs.'];
enddatelist = [locs.',t(end)];
paramestimate = zeros(4,length(Idata));
xlabel('time (days)');
ylabel('Number of Infected Population');

writematrix(paramestimate,'mcmc_result.csv')
%%
clear model data params options
mco_number = length(startdatelist);
% create panel
p = panel();

% layout a variety of sub-panels
% p.pack('h', {1/4 []})
p.pack(4, 1);

% set margins
p.de.margin = 2;
p.marginleft = 5;
p.margin = [13 10 8 8];

% and some properties
p.fontsize = 8;

p(1,1).select();
ylabel('S')
p(2,1).select();
ylabel('I')
p(3,1).select();
ylabel('R')
p(4,1).select();
ylabel('D')

mcmcresult = zeros(length(startdatelist),6);
% for mcophase =  1:length(startdatelist)
for mcophase =  1:length(startdatelist)

    startdate = startdatelist(length(startdatelist) - mcophase);
    enddate = enddatelist(length(startdatelist));
    mcoenddate = enddatelist(length(startdatelist) - mcophase);
%     startdate = startdatelist(mcophase);
%     enddate = enddatelist(mcophase);
%     mcoenddate = enddatelist(mcophase);
    Idata2 = Idata(startdate:enddate);
    Rdata2 = Rdata(startdate:enddate);
    Ddata2 = Ddata(startdate:enddate);
    Idata3 = Idata(startdate:mcoenddate);
    Rdata3 = Rdata(startdate:mcoenddate);
    Ddata3 = Ddata(startdate:mcoenddate);
    Sdata3 = Sdata(startdate:mcoenddate);
    Sdata2 = Sdata(startdate:enddate);
    tdata2 = t(startdate:enddate);
    
    data = [[tdata2.'],[Idata2],[Rdata2],[Ddata2]];

    %%
    % Initial concentrations are saved in |data| to be used in sum of
    % squares function.

    I0 = data(1,2);
    R0 = data(1,3); D0 = data(1,4);
    S0 = Ptotal - R0 - D0 - I0;

    %%
    % Refine the first guess for the parameters with |fminseacrh| and
    % calculate residual variance as an estimate of the model error variance.
    Cdata3 = Idata3 + Rdata3 + Ddata3;
    beta_start = mean(diff(Cdata3)./Idata3(1:end-1).*Sdata3(1:end-1)./Ptotal);
    r_start = mean((diff(Rdata3)+ mu*mean(Rdata3))./Idata3(1:end-1));
    d_start = mean(diff(Ddata3)./Idata3(1:end-1));
    gamma_start = 0.0000000001;
    k00 = [beta_start,r_start,d_start,gamma_start]';
    [k0,ss0] = fminsearch(@sirdss,k00,[],data);
    mse = ss0/(length(data));

    %%
    params = {
        {'beta_start', k0(1), 0,1}
        {'r_start', k0(2), 0,1}
        {'d_start', k0(3), 0,1}
        {'gamma_start', k0(4), 0,0.001}
        };

    model.ssfun = @sirdss;
    model.sigma2 = mse;
    options.nsimu = 5000;
    options.updatesigma = 1;
    options.method = 'dram';

    %%
    [results,chain,s2chain] = mcmcrun(model,data,params,options);
%     
%     results = [];
%     [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
%     chain_stat = chainstats(chain,results);
%     [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
%     chain_stat = chainstats(chain,results);
%     [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
%     chain_stat = chainstats(chain,results);

    figure(2); clf
    mcmcplot(chain,[],results,'chainpanel')
    subplot(2,2,4)
    mcmcplot(sqrt(s2chain),[],[],'dens',2)
    title('error std')

    %%
    % Function |chainstats| lists some statistics, including the
    % estimated Monte Carlo error of the estimates.
    chain_stat = chainstats(chain,results);
    beta_est = chain_stat(1,1);
    r_est = chain_stat(2,1);
    d_est = chain_stat(3,1);
    gamma_est = chain_stat(4,1);
    
    [Smodel,Imodel,Rmodel, Dmodel] = sirdsrk4(tdata2.',S0,I0,R0,D0,Ptotal,beta_est,r_est,d_est,gamma_est,0.05);
%     figure(3); clf
%     set(gca,'xticklabel',[])
    tdayssplit = tdays(startdate:enddate);
    p(1,1).select();
    plot(tdayssplit,Sdata2,'ks',tdayssplit,Smodel,'k-');
    p(2,1).select();
    plot(tdayssplit,Idata2,'bs',tdayssplit,Imodel,'b-');
    p(3,1).select();
    plot(tdayssplit,Rdata2,'cs',tdayssplit,Rmodel,'c-');
    p(4,1).select();
    plot(tdayssplit,Ddata2,'rs',tdayssplit,Dmodel,'r-');
%     hold on
    prestartdate = sum(paramestimate(1,:) == 0);
    paramestimate(1,startdate:prestartdate) = beta_est;
    paramestimate(2,startdate:prestartdate) = r_est;
    paramestimate(3,startdate:prestartdate) = d_est;
    paramestimate(4,startdate:prestartdate) = gamma_est;
    paramestimate(5,startdate:prestartdate) = 0.05;
    writematrix(paramestimate,'mcmc_result.csv')
end

hold off
%%
% ARMSE = sqrt(mean((Adata2 - Amodel).^2))
% RRMSE = sqrt(mean((Rdata2 - Rmodel).^2))
% DRMSE = sqrt(mean((Ddata2 - Dmodel).^2))