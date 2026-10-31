%%
clear model data params options
datafile = readmatrix('si3rd_data.csv');
Ptotal = 3.153*10^7;%malaysia
global paramestimate
I1data = datafile(:,3);
I2data = datafile(:,5);
I3data = datafile(:,6);
Rdata = datafile(:,2);
Ddata = datafile(:,4);
Sdata = datafile(:,7);
tdata = datafile(:,1);
t = 1:1:length(I1data);
% mu = 0.000021;
mu = 0;
%%
% slopedetection
dIdt = diff([eps; I1data(:)])./diff([eps; t(:)]);
sdIdt = smoothdata(dIdt,'gaussian',32);
d2Idt2 = diff([eps; sdIdt(:)])./diff([eps; t(:)]);
sd2Idt2 = smoothdata(d2Idt2,'gaussian',25);
[pks,locs] = findpeaks(abs(sd2Idt2));
dIdt = diff([eps; I2data(:)])./diff([eps; t(:)]);
sdIdt = smoothdata(dIdt,'gaussian',32);
d2Idt2 = diff([eps; sdIdt(:)])./diff([eps; t(:)]);
sd2Idt2 = smoothdata(d2Idt2,'gaussian',25);
[pks2,locs2] = findpeaks(abs(sd2Idt2));
locs = union(locs,locs2);
tdays = datetime(2020,3,17) + caldays(1:length(t));
% plot(tdays,Idata,'.');
% xtickformat('dd-MM-yy')
% xline([locs]);
startdatelist = [1,locs.'];
enddatelist = [locs.',t(end)];
paramestimate = zeros(15,length(I1data));
writematrix(paramestimate,'mcmc_result.csv')
%%
clear model data params options
mco_number = length(startdatelist);
% create panel
p = panel();

% layout a variety of sub-panels
% p.pack('h', {1/4 []})
p.pack(6, 1);

% set margins
p.de.margin = 2;
p.marginleft = 5;
p.margin = [13 10 8 8];

% and some properties
p.fontsize = 8;

p(1,1).select();
ylabel('S')
p(2,1).select();
ylabel('I1')
p(3,1).select();
ylabel('I2')
p(4,1).select();
ylabel('I3')
p(5,1).select();
ylabel('R')
p(6,1).select();
ylabel('D')

mcmcresult = zeros(length(startdatelist),15);
% for mcophase =  1:length(startdatelist)
for mcophase =  1:length(startdatelist)

    startdate = startdatelist(length(startdatelist) - mcophase);
    enddate = enddatelist(length(startdatelist));
    mcoenddate = enddatelist(length(startdatelist) - mcophase);
%     startdate = startdatelist(mcophase);
%     enddate = enddatelist(mcophase);
%     mcoenddate = enddatelist(mcophase);
    I1data2 = I1data(startdate:enddate);
    I2data2 = I2data(startdate:enddate);
    I3data2 = I3data(startdate:enddate);
    Rdata2 = Rdata(startdate:enddate);
    Ddata2 = Ddata(startdate:enddate);
    I1data3 = I1data(startdate:mcoenddate);
    I2data3 = I2data(startdate:mcoenddate);
    I3data3 = I3data(startdate:mcoenddate);
    Rdata3 = Rdata(startdate:mcoenddate);
    Ddata3 = Ddata(startdate:mcoenddate);
    Sdata3 = Sdata(startdate:mcoenddate);
    Sdata2 = Sdata(startdate:enddate);
    tdata2 = t(startdate:enddate);
    
    data = [[tdata2.'],[I1data2],[I2data2],[I3data2],[Rdata2],[Ddata2]];

    %%
    % Initial concentrations are saved in |data| to be used in sum of
    % squares function.

    I10 = data(1,2); I20 = data(1,3); I30 = data(1,4);
    R0 = data(1,5); D0 = data(1,6);
    S0 = Ptotal - R0 - D0 - I10 - I20 - I30;

    %%
    % Refine the first guess for the parameters with |fminseacrh| and
    % calculate residual variance as an estimate of the model error variance.
    Cdata3 = I1data3 + Rdata3 + Ddata3;
    p1_start = mean(I1data3)/mean(I1data3 + I2data3 + I3data3);
    p2_start = mean(I2data3)/mean(I1data3 + I2data3 + I3data3);
    p3_start = mean(I3data3)/mean(I1data3 + I2data3 + I3data3);
    beta1_start = mean(diff(Cdata3)./I1data3(1:end-1).*Sdata3(1:end-1)./Ptotal);
    r1_start = mean((diff(Rdata3)+mu*mean(Rdata3))./I1data3(1:end-1));
    d1_start = mean(diff(Ddata3)./I1data3(1:end-1));
    r2_start = r1_start*0.7;
    r3_start = r1_start*0.8;
    d2_start = d1_start*1.2;
    d3_start = d1_start*2;
    beta2_start = mean((diff(I2data3) + r2_start.*p2_start.*diff(Rdata3) + d2_start.*p2_start.*diff(Ddata3))./I1data3(1:end-1).*Sdata3(1:end-1)./Ptotal);
    beta3_start = mean((diff(I3data3) + r3_start.*p3_start.*diff(Rdata3) + d3_start.*p3_start.*diff(Ddata3))./I1data3(1:end-1).*Sdata3(1:end-1)./Ptotal);
    gamma_start = 0.0001;
    k00 = [beta1_start,beta2_start,beta3_start,r1_start,r2_start,r3_start,d1_start,d2_start,d3_start,p1_start,p2_start,p3_start,gamma_start]';
    [k0,ss0] = fminsearch(@si3rdss,k00,[],data);
    mse = ss0/(length(data));
    
    

    %%
    params = {
        {'beta1_start', k0(1), 0,1}
        {'beta2_start', k0(2), 0,1}
        {'beta3_start', k0(3), 0,1}
        {'r1_start', k0(4), 0,1}
        {'r2_start', k0(5), 0,1}
        {'r3_start', k0(6), 0,1}
        {'d1_start', k0(7), 0,1}
        {'d2_start', k0(8), 0,1}
        {'d3_start', k0(9), 0,1}
        {'p1_start', k0(10), 0,1}
        {'p2_start', k0(11), 0,1}
        {'p3_start', k0(12), 0,1}
        {'gamma_start', k0(13), 0,0.001}
        };

    model.ssfun = @si3rdss;
    model.sigma2 = mse;
    options.nsimu = 5000;
    options.updatesigma = 1;
    options.method = 'dram';

    %%
%     [results,chain,s2chain] = mcmcrun(model,data,params,options);
%     
    results = [];
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
    chain_stat = chainstats(chain,results);
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
    chain_stat = chainstats(chain,results);
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
    chain_stat = chainstats(chain,results);

%     figure(4); clf
%     mcmcplot(chain,[],results,'chainpanel')
%     subplot(4,3,12)
%     mcmcplot(sqrt(s2chain),[],[],'dens',2)
%     title('error std')

    %%
    % Function |chainstats| lists some statistics, including the
    % estimated Monte Carlo error of the estimates.
    chain_stat = chainstats(chain,results);
    beta1_est = chain_stat(1,1);
    beta2_est = chain_stat(2,1);
    beta3_est = chain_stat(3,1);
    r1_est = chain_stat(4,1);
    r2_est = chain_stat(5,1);
    r3_est = chain_stat(6,1);
    d1_est = chain_stat(7,1);
    d2_est = chain_stat(8,1);
    d3_est = chain_stat(9,1);
    p1_est = chain_stat(10,1);
    p2_est = chain_stat(11,1);
    p3_est = chain_stat(12,1);
    gamma_est = chain_stat(13,1);
    
    [Smodel,I1model,I2model,I3model,Rmodel, Dmodel] = si3rdsrk4(tdata2.',S0,I10,I20,I30,R0,D0,Ptotal,beta1_est,beta2_est,beta3_est,r1_est,r2_est,r3_est,d1_est,d2_est,d3_est,p1_est,p2_est,p3_est,gamma_est);
%     figure(3); clf
%     set(gca,'xticklabel',[])
    tdayssplit = tdays(startdate:enddate);
    p(1,1).select();
    plot(tdayssplit,Sdata2,'ks',tdayssplit,Smodel,'k-');
    p(2,1).select();
    plot(tdayssplit,I1data2,'bs',tdayssplit,I1model,'b-');
    p(3,1).select();
    plot(tdayssplit, I2data2, 'ms', tdayssplit, I2model, 'm-'); 
    p(4,1).select();
    plot(tdayssplit, I3data2, 'ms', tdayssplit, I3model, 'm-');
    p(5,1).select();
    plot(tdayssplit,Rdata2,'gs',tdayssplit,Rmodel,'g-');
    p(6,1).select();
    plot(tdayssplit,Ddata2,'rs',tdayssplit,Dmodel,'r-');
%     hold on
    prestartdate = sum(paramestimate(1,:) == 0);
    paramestimate(1,startdate:prestartdate) = beta1_est;
    paramestimate(2,startdate:prestartdate) = beta2_est;
    paramestimate(3,startdate:prestartdate) = beta3_est;
    paramestimate(4,startdate:prestartdate) = r1_est;
    paramestimate(5,startdate:prestartdate) = r2_est;
    paramestimate(6,startdate:prestartdate) = r3_est;
    paramestimate(7,startdate:prestartdate) = d1_est;
    paramestimate(8,startdate:prestartdate) = d2_est;
    paramestimate(9,startdate:prestartdate) = d3_est;
    paramestimate(10,startdate:prestartdate) = p1_est;
    paramestimate(11,startdate:prestartdate) = p2_est;
    paramestimate(12,startdate:prestartdate) = p3_est;
    paramestimate(13,startdate:prestartdate) = gamma_est;
    paramestimate(14,startdate:prestartdate) = 0.05;
    
end

hold off
%%
% ARMSE = sqrt(mean((Adata2 - Amodel).^2))
% RRMSE = sqrt(mean((Rdata2 - Rmodel).^2))
% DRMSE = sqrt(mean((Ddata2 - Dmodel).^2))