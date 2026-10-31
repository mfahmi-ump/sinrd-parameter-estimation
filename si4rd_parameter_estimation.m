%%
clear model data params options
datafile = readmatrix('si4rd_data_adjusted.csv');
Ptotal = 3.153*10^7;%malaysia
global paramestimate
% 1:500
% 480:1000
% 980:1377

% I1data = datafile(980:1377,2);
% I2data = datafile(980:1377,5);
% I3data = datafile(980:1377,6);
% I4data = datafile(980:1377,7);
% Rdata = datafile(980:1377,3);
% Ddata = datafile(980:1377,4);
% Sdata = datafile(980:1377,8);
% tdata = datafile(980:1377,1);
% Hdata = datafile(980:1377,9);

% I1data = datafile(480:1000,2);
% I2data = datafile(480:1000,5);
% I3data = datafile(480:1000,6);
% I4data = datafile(480:1000,7);
% Rdata = datafile(480:1000,3);
% Ddata = datafile(480:1000,4);
% Sdata = datafile(480:1000,8);
% tdata = datafile(480:1000,1);
% Hdata = datafile(480:1000,9);

% I1data = datafile(1:500,2);
% I2data = datafile(1:500,5);
% I3data = datafile(1:500,6);
% I4data = datafile(1:500,7);
% Rdata = datafile(1:500,3);
% Ddata = datafile(1:500,4);
% Sdata = datafile(1:500,8);
% tdata = datafile(1:500,1);
% Hdata = datafile(1:500,9);

% I1data = datafile(712:1088,2);
% I2data = datafile(712:1088,5);
% I3data = datafile(712:1088,6);
% I4data = datafile(712:1088,7);
% Rdata = datafile(712:1088,3);
% Ddata = datafile(712:1088,4);
% Sdata = datafile(712:1088,8);
% tdata = datafile(712:1088,1);
% Hdata = datafile(712:1088,9);

I1data = datafile(:,2);
I2data = datafile(:,5);
I3data = datafile(:,6);
I4data = datafile(:,7);
Rdata = datafile(:,3);
Ddata = datafile(:,4);
Sdata = datafile(:,8);
tdata = datafile(:,1);
Hdata = datafile(:,9);

I1data = max(I1data,0.01);
I2data = max(I2data,0.01);
I3data = max(I3data,0.01);
I4data = max(I4data,0.01);

t = 1:1:length(I1data);
mu = 0;
%%
% slopedetection
dIdt = diff([eps; I1data(:)])./diff([eps; t(:)]);
sdIdt = smoothdata(dIdt,'gaussian',42);
d2Idt2 = diff([eps; sdIdt(:)])./diff([eps; t(:)]);
sd2Idt2 = smoothdata(d2Idt2,'gaussian',30);
[pks,locs] = findpeaks(abs(sd2Idt2));
dIdt = diff([eps; I2data(:)])./diff([eps; t(:)]);
sdIdt = smoothdata(dIdt,'gaussian',42);
d2Idt2 = diff([eps; sdIdt(:)])./diff([eps; t(:)]);
sd2Idt2 = smoothdata(d2Idt2,'gaussian',30);
[pks2,locs2] = findpeaks(abs(sd2Idt2));
locs = union(locs,locs2);
tdays = datetime(2020,3,24) + caldays(1:length(t));
% plot(tdays,I1data,'.');
% xtickformat('dd-MM-yy')
% xline([locs]);
startdatelist = [1,locs.'];
enddatelist = [locs.',t(end)];
paramestimate = zeros(19,length(I1data));
%%
clear model data params options
% create panel
p = panel();

% layout a variety of sub-panels
% p.pack('h', {1/4 []})
p.pack(7, 1);

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
ylabel('I4')
p(6,1).select();
ylabel('R')
p(7,1).select();
ylabel('D')
for mcophase =  1:length(startdatelist)
% for mcophase = 2:3

    startdate = startdatelist(length(startdatelist) - mcophase);
    enddate = enddatelist(length(startdatelist));
    mcoenddate = enddatelist(length(startdatelist) - mcophase);
    
%     startdate = startdatelist(mcophase);
%     enddate = enddatelist(mcophase);
%     mcoenddate = enddatelist(mcophase);

    I1data2 = I1data(startdate:enddate);
    I2data2 = I2data(startdate:enddate);
    I3data2 = I3data(startdate:enddate);
    I4data2 = I4data(startdate:enddate);
    Rdata2 = Rdata(startdate:enddate);
    Ddata2 = Ddata(startdate:enddate);
    Hdata2 = Hdata(startdate:enddate);
    
    I1data3 = I1data(startdate:mcoenddate);
    I2data3 = I2data(startdate:mcoenddate);
    I3data3 = I3data(startdate:mcoenddate);
    I4data3 = I4data(startdate:mcoenddate);
    Rdata3 = Rdata(startdate:mcoenddate);
    Ddata3 = Ddata(startdate:mcoenddate);
    Sdata3 = Sdata(startdate:mcoenddate);
    Hdata3 = Hdata(startdate:mcoenddate);
    
    Sdata2 = Sdata(startdate:enddate);
    tdata2 = t(startdate:enddate);
        
    data = [[tdata2.'],[I1data2],[I2data2],[I3data2],[I4data2],[Rdata2],[Ddata2],[Sdata2]];

    %%
    % Initial concentrations are saved in |data| to be used in sum of
    % squares function.

    I10 = data(1,2); 
    I20 = data(1,3); I30 = data(1,4); I40 = data(1,5);
    R0 = data(1,6); D0 = data(1,7);
    S0 = data(1,8);

    Cdata3 = I1data3 + I2data3 + I3data3 + I4data3 + Rdata3 + Ddata3;
    Itotal3 = I1data3 + I2data3 + I3data3 + I4data3;
    p1_start = mean(diff(Cdata3) - diff(Hdata3))/mean(diff(Cdata3));
    p2_start = (1-p1_start)*mean(I2data3)/mean(I2data3 + I3data3 + I4data3);
    p3_start = (1-p1_start)*mean(I3data3)/mean(I2data3 + I3data3 + I4data3);
    p4_start = (1-p1_start)*mean(I4data3)/mean(I2data3 + I3data3 + I4data3);
    beta1_start = p1_start.*mean(diff(Cdata3)./Itotal3(1:end-1)./Sdata3(1:end-1).*Ptotal);
    beta2_start = p2_start.*mean(diff(Cdata3)./Itotal3(1:end-1)./Sdata3(1:end-1).*Ptotal);
    beta3_start = p3_start.*mean(diff(Cdata3)./Itotal3(1:end-1)./Sdata3(1:end-1).*Ptotal);
    beta4_start = p4_start.*mean(diff(Cdata3)./Itotal3(1:end-1)./Sdata3(1:end-1).*Ptotal);
    
    gamma_start = 0.0000001;
    r1_start = mean(diff(-I1data3) + p1_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Rdata3./(gamma_start+1).* I1data3);
    r2_start = mean(diff(-I2data3) + p2_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Rdata3./(gamma_start+1).* I2data3);
    r3_start = mean(diff(-I3data3) + p3_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Rdata3./(gamma_start+1).* I3data3);
    r4_start = mean(diff(-I4data3) + p4_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Rdata3./(gamma_start+1).* I4data3);
    d1_start = mean(diff(-I1data3) + p1_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Ddata3.* I1data3);
    d2_start = mean(diff(-I2data3) + p2_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Ddata3.* I2data3);
    d3_start = mean(diff(-I3data3) + p3_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Ddata3.* I3data3);
    d4_start = mean(diff(-I4data3) + p4_start.*diff(Cdata3))./ mean((Rdata3.*(gamma_start+1) + Ddata3)./ Ddata3.* I4data3);
    
    sigma_start = 0.05;
    k00 = [p1_start,p2_start,p3_start,p4_start,beta1_start,beta2_start,beta3_start,beta4_start,r1_start,r2_start,r3_start,r4_start,d1_start,d2_start,d3_start,d4_start,gamma_start]';
    [k0,ss0] = fminsearch(@si4rdss,k00,[],data);
    mse = ss0/7;

    %%
    params = {
        {'p1_start', max(k0(1),0), 0,1}
        {'p2_start', max(k0(2),0), 0,1}
        {'p3_start', max(k0(3),0), 0,1}
        {'p4_start', max(k0(4),0), 0,1}
        {'beta1_start', max(k0(5),0), 0,1}
        {'beta2_start', max(k0(6),0), 0,1}
        {'beta3_start', max(k0(7),0), 0,1}
        {'beta4_start', max(k0(8),0), 0,1}
        {'r1_start', max(k0(9),0), 0,1}
        {'r2_start', max(k0(10),0), 0,1}
        {'r3_start', max(k0(11),0), 0,1}
        {'r4_start', max(k0(12),0), 0,1}
        {'d1_start', max(k0(13),0), 0,1}
        {'d2_start', max(k0(14),0), 0,1}
        {'d3_start', max(k0(15),0), 0,1}
        {'d4_start', max(k0(16),0), 0,1}
        {'gamma_start', max(k0(17),0), 0,0.000001}
        };

    model.ssfun = @si4rdss;
    model.sigma2 = mse;
    options.nsimu = 1000;
    options.updatesigma = 1;
    options.method = 'dram';

    %%
%     [results,chain,s2chain] = mcmcrun(model,data,params,options);
%     results = [];
    try
        [results,chain,s2chain] = mcmcrun(model,data,params,options);
        chain_stat = chainstats(chain,results);
        [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
        chain_stat = chainstats(chain,results);
        [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
        chain_stat = chainstats(chain,results);
    catch
        [results,chain,s2chain] = mcmcrun(model,data,params,options);
    end
    p1_est = chain_stat(1,1);
    p2_est = chain_stat(2,1);
    p3_est = chain_stat(3,1);
    p4_est = chain_stat(4,1);
    beta1_est = chain_stat(5,1);
    beta2_est = chain_stat(6,1);
    beta3_est = chain_stat(7,1);
    beta4_est = chain_stat(8,1);
    r1_est = chain_stat(9,1);
    r2_est = chain_stat(10,1);
    r3_est = chain_stat(11,1);
    r4_est = chain_stat(12,1);
    d1_est = chain_stat(13,1);
    d2_est = chain_stat(14,1);
    d3_est = chain_stat(15,1);
    d4_est = chain_stat(16,1);
    gamma_est = chain_stat(17,1);
    %%
    % Function |chainstats| lists some statistics, including the
    % estimated Monte Carlo error of the estimates.
    chain_stat = chainstats(chain,results);
    %%
%     ymax = max(max(Idata2),max(Rdata2));
    xlength = length(I1data2);
    [Smodel,I1model,I2model,I3model,I4model,Rmodel, Dmodel] = si4rdsrk4(tdata2.',S0,I10,I20,I30,I40,R0,D0,Ptotal,p1_est,p2_est,p3_est,p4_est,beta1_est,beta2_est,beta3_est,beta4_est,r1_est,r2_est,r3_est,r4_est,d1_est,d2_est,d3_est,d4_est,gamma_est,0.05);
%     figure(3); clf
%     set(gca,'xticklabel',[])
    tdayssplit = tdays(startdate:enddate);
    p(1,1).select();
    plot(tdayssplit,Sdata2,'ks',tdayssplit,Smodel,'k-');
    p(2,1).select();
    plot(tdayssplit,I1data2,'bs',tdayssplit,I1model,'b-');
    p(3,1).select();
    plot(tdayssplit,I2data2,'gs',tdayssplit,I2model,'g-');
    p(4,1).select();
    plot(tdayssplit,I3data2,'cs',tdayssplit,I3model,'c-');
    p(5,1).select();
    plot(tdayssplit,I4data2,'ms',tdayssplit,I4model,'m-');
    p(6,1).select();
    plot(tdayssplit,Rdata2,'cs',tdayssplit,Rmodel,'c-');
    p(7,1).select();
    plot(tdayssplit,Ddata2,'rs',tdayssplit,Dmodel,'r-');
%     hold on
    prestartdate = sum(paramestimate(1,:) == 0);
    paramestimate(1,startdate:prestartdate) = chain_stat(1,1);
    paramestimate(2,startdate:prestartdate) = chain_stat(2,1);
    paramestimate(3,startdate:prestartdate) = chain_stat(3,1);
    paramestimate(4,startdate:prestartdate) = chain_stat(4,1);
    paramestimate(5,startdate:prestartdate) = chain_stat(5,1);
    paramestimate(6,startdate:prestartdate) = chain_stat(6,1);
    paramestimate(7,startdate:prestartdate) = chain_stat(7,1);
    paramestimate(8,startdate:prestartdate) = chain_stat(8,1);
    paramestimate(9,startdate:prestartdate) = chain_stat(9,1);
    paramestimate(10,startdate:prestartdate) = chain_stat(10,1);
    paramestimate(11,startdate:prestartdate) = chain_stat(11,1);
    paramestimate(12,startdate:prestartdate) = chain_stat(12,1);
    paramestimate(13,startdate:prestartdate) = chain_stat(13,1);
    paramestimate(14,startdate:prestartdate) = chain_stat(14,1);
    paramestimate(15,startdate:prestartdate) = chain_stat(15,1);
    paramestimate(16,startdate:prestartdate) = chain_stat(16,1);
    paramestimate(17,startdate:prestartdate) = chain_stat(17,1);
    paramestimate(18,startdate:prestartdate) = 0.05;
end


