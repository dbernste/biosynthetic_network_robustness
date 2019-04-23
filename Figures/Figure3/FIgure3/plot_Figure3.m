clear all
close all
clc

[num,txt,raw] = xlsread('./Figure3_Data.xlsx');

f_max(1,1:50) = num(1:50,70);
f_max(2,1:50) = num(51:100,70);
f_max(3,1:50) = num(101:150,70);
f_max(4,1:50) = num(151:200,70);
f_max(5,1:50) = num(201:250,70);

f_min(1,1:50) = num(1:50,69);
f_min(2,1:50) = num(51:100,69);
f_min(3,1:50) = num(101:150,69);
f_min(4,1:50) = num(151:200,69);
f_min(5,1:50) = num(201:250,69);

m_PM = num(1:end,1:68);

% normalize and calculate statistics for FBA
f_max_n = f_max./max(max(f_max));
f_min_n = f_min./max(max(f_min));
m_f_max_n = mean(f_max_n,2);
s_f_max_n = std(f_max_n,'',2);
m_f_min_n = mean(f_min_n,2);
s_f_min_n = std(f_min_n,'',2);

% calculate error rates for PM
% clear zeros
tmp = zeros(size(m_PM,2),1);
for I = 1:size(m_PM,2)
    if m_PM(:,I)==zeros(size(m_PM,1),1)
        tmp(I) = 1;
    end
end
nz_ind = find(tmp==0);
clear tmp

m_PM_max_nz = max(m_PM(:,nz_ind));
% calculate error rates for PM
m_PM_nz = m_PM(:,nz_ind);
%m_PM_nz_n = m_PM_nz./m_PM_max_nz;
n_nz = zeros(size(m_PM_nz,1),1);
for I = 1:size(m_PM_nz,1)
    m_PM_n_nz = m_PM_nz(I,:);
    n_nz(I) = norm(ones(size(m_PM_max_nz))-m_PM_n_nz,1);
end
n_n_nz = n_nz./length(m_PM_max_nz);
m_n_n_nz = zeros(5,1);
s_n_n_nz = zeros(5,1);
for I = 1:5
    m_n_n_nz(I) = mean(n_n_nz(((I-1)*50+1):(I)*50));
    s_n_n_nz(I) = std(n_n_nz(((I-1)*50+1):(I)*50));
end


%% Plot Results
x = [0;4;16;64;256;1024]./2583;

figure(1)
h(2) = semilogx(x,[1;m_f_max_n],'r.','markersize',20);
hold on
yp = [[m_f_max_n]-[s_f_max_n]/sqrt(50);flipud([m_f_max_n]+[s_f_max_n]/sqrt(50));[m_f_max_n(1)]-[s_f_max_n(1)]/sqrt(50)];
xp = [x(2:end);flipud(x(2:end));x(2)];
patch(xp,yp,[1 0.9 0.9],'EdgeColor','none')
h(2) = semilogx(x,[1;m_f_max_n],'r.','markersize',20);
semilogx(x,[1;m_f_max_n],'r')
semilogx(x,[1;m_f_max_n]-[0;s_f_max_n]/sqrt(50),'r:')
semilogx(x,[1;m_f_max_n]+[0;s_f_max_n]/sqrt(50),'r:')

yp = [[m_f_min_n]-[s_f_min_n]/sqrt(50);flipud([m_f_min_n]+[s_f_min_n]/sqrt(50));[m_f_min_n(1)]-[s_f_min_n(1)]/sqrt(50)];
xp = [x(2:end);flipud(x(2:end));x(2)];
patch(xp,yp,[0.9 0.9 1],'EdgeColor','none')
h(2) = semilogx(x,[1;m_f_max_n],'r.','markersize',20);
h(1) = semilogx(x,[1;m_f_min_n],'b.','markersize',20);
semilogx(x,[1;m_f_min_n],'b')
semilogx(x,[1;m_f_min_n]-[0;s_f_min_n]/sqrt(50),'b:')
semilogx(x,[1;m_f_min_n]+[0;s_f_min_n]/sqrt(50),'b:')

yp = [1-[m_n_n_nz]+[s_n_n_nz]/sqrt(50);flipud(1-[m_n_n_nz]-[s_n_n_nz]/sqrt(50));1-[m_n_n_nz(1)]+[s_n_n_nz(1)]/sqrt(50)];
xp = [x(2:end);flipud(x(2:end));x(2)];
patch(xp,yp,[1 0.9 1],'EdgeColor','none')
h(3) = semilogx(x,1-[0;m_n_n_nz],'m.','markersize',20);
semilogx(x,1-[0;m_n_n_nz],'m')
semilogx(x,1-[0;m_n_n_nz]+[0;s_n_n_nz]/sqrt(50),'m:')
semilogx(x,1-[0;m_n_n_nz]-[0;s_n_n_nz]/sqrt(50),'m:')

ylabel('quantitative difference accuracy','fontsize',16)
xlabel('fraction of reactions randomly removed','fontsize',16)

%legend(h,'FBA (minimal media)','FBA (complete media)','Producibility Metric')

axis([0 1 0 1])

%% Plot Number of Correct
thresh = 0.1;
p_max = sum(f_max_n>thresh,2);
p_min = sum(f_min_n>thresh,2);

p_n_n_nz = zeros(5,1);
p_n_n_nz6 = zeros(5,1);
for I = 1:5 % for each set
    for J = 1:50 % for each model
        vect = m_PM_nz((I-1)*50+J,:);
        if min(vect) > 0.1
            p_n_n_nz(I) = p_n_n_nz(I) + 1;
        end
        if min(vect) > 0.6
            p_n_n_nz6(I) = p_n_n_nz6(I) + 1;
        end
    end
end

figure(2)
h(2) = semilogx(x,[50;p_max/50],'r.','markersize',20);
hold on
semilogx(x,[50;p_max]./50,'r')
h(1) = semilogx(x,[50;p_min]./50,'b.','markersize',20);
semilogx(x,[50;p_min]./50,'b')
h(3) = semilogx(x,[50;p_n_n_nz]./50,'m.','markersize',20);
semilogx(x,[50;p_n_n_nz]./50,'m')
h(4) = semilogx(x,[50;p_n_n_nz6]./50,'m.','markersize',20);
semilogx(x,[50;p_n_n_nz6]./50,'m:')

ylabel('biomass production accuracy','fontsize',16)
xlabel('fraction of reactions randomly removed','fontsize',16)

legend(h(1:3),'FBA (minimal medium)','FBA (complete medium)','Producibility Metric')

axis([0 1 0 1])