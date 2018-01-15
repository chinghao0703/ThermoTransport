% This code generates the phase diagram both in the main text and SM. The
% data is pregenerated by solving the coupled chemical kinetics equations
% using computing clusters over the parameter regimes of interest.

clear all;



%load('combinedData_vol5theta10_c2m2NewICnew');  cindex=2;  %so C_tot = 20 nM
load('combinedData_vol5theta10_c5m2NewICnew');  cindex=5;  %so C_tot = 50 nM
%load('combinedData_vol5theta10_c10m2NewICnew'); cindex=10; %so C_tot = 100 nM



if cindex==2
    str2='[C]_{all}=20';
elseif cindex==5
    str2='[C]_{all}=50';
elseif cindex==10
    str2='[C]_{all}=100';
end
    
  


[Trans, N2]=meshgrid(trans_linear,N2_linear); % (transporter, NTF2)

%FigHandle = figure('Position', [100, 100, 1049, 895]);
hFig = figure(1);
set(hFig, 'Position', [100, 96, 1086, 454])
f1=subplot(121)
pcolor(Trans,N2,CargoNCmaxc10');
line([0,200],[100,100],'LineStyle','--','Color','w','LineWidth',4)
shading flat;
colorbar;
colormap(f1,parula);
xlabel('[Importin] (nM)', 'fontsize', 20);
ylabel('[NTF2] (nM)', 'fontsize', 20);
title(['[C]_{nu}/[C]_{cyto} maximum']);
set(gca,'Fontsize',18);

f2=subplot(122)
pcolor(Trans(5:151,5:101),N2(5:151,5:101),EPallc10max(5:101,5:151)');
%pcolor(Trans(5:151,5:101),N2(5:151,5:101),EPallc10(5:101,5:151)');
line([20,200],[100,100],'LineStyle','--','Color','w','LineWidth',4)
shading interp;
colorbar;
colormap(f2,jet);
xlabel('[Importin] (nM)', 'fontsize', 20);
ylabel('[NTF2] (nM)', 'fontsize', 20);
title(['EP/k_BT']);
axis([20 200 10 300])
set(gca,'Fontsize',18);


str2


