clc
clear all
close all

% Load Data
load('./Figure4_Sup2_Data.mat')

%% Fit Model (linear)
tbl = table(m_m_PM,log_gSize,'VariableNames',{'m_m_PM','log_gSize'});
lm = fitlm(tbl,'m_m_PM~log_gSize');
xl = [min(log_gSize),max(log_gSize)];
yl = lm.Coefficients{1,1} + lm.Coefficients{2,1}.*xl;
[aic0, bic0] = aicbic(lm.LogLikelihood,lm.NumCoefficients,lm.NumObservations);

% Fit Model with additional categorical information
% phylyum
tbl1 = tbl;
tbl1.phy = nominal(phylum);
lm1 = fitlm(tbl1,'m_m_PM~phy+log_gSize');
[aic(1), bic(1)] = aicbic(lm1.LogLikelihood,lm1.NumCoefficients,lm1.NumObservations);
% class
tbl2 = tbl;
tbl2.cla = nominal(class);
lm2 = fitlm(tbl2,'m_m_PM~cla+log_gSize');
[aic(2), bic(2)] = aicbic(lm2.LogLikelihood,lm2.NumCoefficients,lm2.NumObservations);
% order
tbl3 = tbl;
tbl3.ord = nominal(order);
lm3 = fitlm(tbl3,'m_m_PM~ord+log_gSize');
[aic(3), bic(3)] = aicbic(lm3.LogLikelihood,lm3.NumCoefficients,lm3.NumObservations);
% family
tbl4 = tbl;
tbl4.fam = nominal(family);
lm4 = fitlm(tbl4,'m_m_PM~fam+log_gSize');
[aic(4), bic(4)] = aicbic(lm4.LogLikelihood,lm4.NumCoefficients,lm4.NumObservations);
% genus
tbl5 = tbl;
tbl5.gen = nominal(genus);
lm5 = fitlm(tbl5,'m_m_PM~gen+log_gSize');
[aic(5), bic(5)] = aicbic(lm5.LogLikelihood,lm5.NumCoefficients,lm5.NumObservations);
% species
tbl6 = tbl;
tbl6.spe = nominal(species);
lm6 = fitlm(tbl6,'m_m_PM~spe+log_gSize');
[aic(6), bic(6)] = aicbic(lm6.LogLikelihood,lm6.NumCoefficients,lm6.NumObservations);
% strain (model)
tbl7 = tbl;
tbl7.str = nominal(strain);
lm7 = fitlm(tbl7,'m_m_PM~str+log_gSize');
[aic(7), bic(7)] = aicbic(lm7.LogLikelihood,lm7.NumCoefficients,lm7.NumObservations);


%% Fit Model (quadratic)
tblb = table(m_m_PM,log_gSize,'VariableNames',{'m_m_PM','log_gSize'});
lmb = fitlm(tblb,'m_m_PM~log_gSize+log_gSize^2');
xc = [min(log_gSize):0.01:max(log_gSize)];
yc = lmb.Coefficients{1,1} + lmb.Coefficients{2,1}.*xc + lmb.Coefficients{3,1}.*xc.^2;
[aic0b, bic0b] = aicbic(lmb.LogLikelihood,lmb.NumCoefficients,lmb.NumObservations);

% Fit Model with additional categorical information
% phylyum
tbl1 = tbl;
tbl1.phy = nominal(phylum);
lm1b = fitlm(tbl1,'m_m_PM~phy+log_gSize+log_gSize^2');
[aicb(1), bicb(1)] = aicbic(lm1b.LogLikelihood,lm1b.NumCoefficients,lm1b.NumObservations);
% class
tbl2 = tbl;
tbl2.cla = nominal(class);
lm2b = fitlm(tbl2,'m_m_PM~cla+log_gSize+log_gSize^2');
[aicb(2), bicb(2)] = aicbic(lm2b.LogLikelihood,lm2b.NumCoefficients,lm2b.NumObservations);
% order
tbl3 = tbl;
tbl3.ord = nominal(order);
lm3b = fitlm(tbl3,'m_m_PM~ord+log_gSize+log_gSize^2');
[aicb(3), bicb(3)] = aicbic(lm3b.LogLikelihood,lm3b.NumCoefficients,lm3b.NumObservations);
% family
tbl4 = tbl;
tbl4.fam = nominal(family);
lm4b = fitlm(tbl4,'m_m_PM~fam+log_gSize+log_gSize^2');
[aicb(4), bicb(4)] = aicbic(lm4b.LogLikelihood,lm4b.NumCoefficients,lm4b.NumObservations);
% genus
tbl5 = tbl;
tbl5.gen = nominal(genus);
lm5b = fitlm(tbl5,'m_m_PM~gen+log_gSize+log_gSize^2');
[aicb(5), bicb(5)] = aicbic(lm5b.LogLikelihood,lm5b.NumCoefficients,lm5b.NumObservations);
% species
tbl6 = tbl;
tbl6.spe = nominal(species);
lm6b = fitlm(tbl6,'m_m_PM~spe+log_gSize+log_gSize^2');
[aicb(6), bicb(6)] = aicbic(lm6b.LogLikelihood,lm6b.NumCoefficients,lm6b.NumObservations);
% strain (model)
tbl7 = tbl;
tbl7.str = nominal(strain);
lm7b = fitlm(tbl7,'m_m_PM~str+log_gSize+log_gSize^2');
[aicb(7), bicb(7)] = aicbic(lm7b.LogLikelihood,lm7b.NumCoefficients,lm7b.NumObservations);

%% Combined Figure
% data
figure()
hold on
% plot data
plot(log_gSize,m_m_PM,'k.')
ylabel('average PM')
xlabel('log genome size')
% plot linear fit
plot(xl,yl,'b','linewidth',2);
% plot quadratic fit
plot(xc,yc,'r','linewidth',2);

figure()
aic1 = [aic0,aic];
plot(aic1,'bo')
hold on
plot(aic1,'b-')
aic1b = [aic0b,aicb];
plot(aic1b,'ro')
hold on
plot(aic1b,'r-')
ylabel('AIC')
xticklabels({'none','phylum','class','order','family','genus','species','strain'})

figure()
bic1 = [bic0,bic];
plot(bic1,'bo')
hold on
plot(bic1,'b-')
bic1b = [bic0b,bicb];
plot(bic1b,'ro')
hold on
plot(bic1b,'r-')
ylabel('BIC')
xticklabels({'none','phylum','class','order','family','genus','species','strain'})