% This code solves the chemical kinetics equations and generates FIG. A2
% in SI. Other figures are generated using a parallelized
% verion of this code (to facilite running on computing clusters).

clear all;


global Xdim; % number of species: 1~14
global tdim; 
global k;    % all rate constants: k(i,1)=k_i^+, k(i,2)=k_i^-

Xdim=14; % 14 species
k=ones(11,2);
tdim=20;
tinterval=[0 tdim];



% X0: initial conditions

% for NTF2 > C  (labeled as 50N)
%X0=[    50   50    0    0    5    0   50   50   25    0   0    5    0   50];
      % X1   X2   X3   X4   X5   X6   X7   X8   X9  X10  X11  X12  X13  X14
      
%for NTF2 = C (labeled as 50)   
%X0=[    5    50    0    0    5    0   50   5   25    0   0    5    0   50];
      % X1   X2   X3   X4   X5   X6   X7   X8   X9  X10  X11  X12  X13  X14

theta = linspace(1,100,50);
speciesequvalue50N=zeros(length(theta),14);
speciesequvalue50=zeros(length(theta),14);



% theta = 10.1. The exact value of theta is not critical to our analysis,
% as shown in FIG. S4.
theta2=1/10.1;
theta5=10.1;

% exp(\Delta F) = exp(\Delta\tilde{F}) in Eq.(S35) (S43)
expFdiff= 50; 


% below specify the kinetics rate constants (SM section II). Note that
% cytosolic and nuclear volume effects will be taken care of in reaction2.m

k=[0.1 ,2.5;...
    expFdiff*theta2,1;...  
    1,1;...
    0.1, 1;...
    expFdiff*theta5, 1;...
    1, 1;...
    20,1;...
    1,1;...
    1,1;...
    1, 20;...
    1,0.2];


for ind = 1:length(theta)
    
    theta2 = 1/ theta(ind);
    theta5 = theta(ind);
    
    k(2,1) = expFdiff*theta2;
    k(5,1) = expFdiff*theta5;
    
    %for NTF2 > C  (labeled as 50N)
    X0=[    50   50    0    0    5    0   50   50   25    0   0    5    0   50];
          % X1   X2   X3   X4   X5   X6   X7   X8   X9  X10  X11  X12  X13  X14
      
    %for NTF2 = C (labeled as 50)   
    X00=[    5    50    0    0    5    0   50   5   25    0   0    5    0   50];
           % X1   X2   X3   X4   X5   X6   X7   X8   X9  X10  X11  X12  X13  X14

    
    
    [T,X] = ode45(@reaction2,tinterval,X0);
    [TT,XX] = ode45(@reaction2,tinterval,X00);
    
    speciesequvalue50N(ind,:)=X(length(T),:);
    speciesequvalue50(ind,:)=XX(length(TT),:);
    
    
    
    
end


totRan50N = speciesequvalue50N(:,2)+speciesequvalue50N(:,3)+speciesequvalue50N(:,10) + ...
    speciesequvalue50N(:,6)+speciesequvalue50N(:,9)+speciesequvalue50N(:,13);

RanGDPfraction50N = (speciesequvalue50N(:,2)+speciesequvalue50N(:,3)+speciesequvalue50N(:,10))./totRan50N;
RanGTPfraction50N = (speciesequvalue50N(:,6)+speciesequvalue50N(:,9)+speciesequvalue50N(:,13))./totRan50N;




totRan50 = speciesequvalue50(:,2)+speciesequvalue50(:,3)+speciesequvalue50(:,10) + ...
    speciesequvalue50(:,6)+speciesequvalue50(:,9)+speciesequvalue50(:,13);

RanGDPfraction50 = (speciesequvalue50(:,2)+speciesequvalue50(:,3)+speciesequvalue50(:,10))./totRan50;
RanGTPfraction50 = (speciesequvalue50(:,6)+speciesequvalue50(:,9)+speciesequvalue50(:,13))./totRan50;

hFig1 = figure(1);
%set(hFig2, 'Position', [100, 96, 1086, 454])
title('[NTF2]_{tot}=100 nM, [C]_{tot}=10 nM')
subplot(2,2,[1,2]);
plot(theta,speciesequvalue50N(:,12)./speciesequvalue50N(:,5),...
     theta,speciesequvalue50N(:,14)./speciesequvalue50N(:,7),...
     theta,speciesequvalue50N(:,11)./speciesequvalue50N(:,4),'LineWidth',2.5);
grid on;
xlabel('\theta', 'fontsize', 20);
legend('[C]_{nu}/[C]_{cyto}','[Im]_{nu}/[Im]_{cyto}','[Im-C]_{nu}/[Im-C]_{cyto}');
set(gca,'Fontsize',18);

    
subplot(2,2,[3,4]);
plot(theta,RanGDPfraction50N,theta,RanGTPfraction50N,'LineWidth',2.5);
grid on;
axis([1 100 0 1]);
xlabel('\theta : free GTP to free GDP ratio', 'fontsize', 20);
legend('n_{RanGDP}','n_{RanGTP}')
ylabel('RanGD(T)P fraction');
set(gca,'Fontsize',18);



hFig2 = figure(2);
%set(hFig2, 'Position', [100, 96, 1086, 454])
title('[NTF2]_{tot}=10 nM, [C]_{tot}=10 nM')
subplot(2,2,[1,2]);
plot(theta,speciesequvalue50(:,12)./speciesequvalue50(:,5),...
     theta,speciesequvalue50(:,14)./speciesequvalue50(:,7),...
     theta,speciesequvalue50(:,11)./speciesequvalue50(:,4),'LineWidth',2.5);
grid on;
xlabel('\theta', 'fontsize', 20);
legend('[C]_{nu}/[C]_{cyto}','[Im]_{nu}/[Im]_{cyto}','[Im-C]_{nu}/[Im-C]_{cyto}');
set(gca,'Fontsize',18);



    
subplot(2,2,[3,4]);
plot(theta,RanGDPfraction50,theta,RanGTPfraction50,'LineWidth',2.5);
grid on;
axis([1 100 0 1]);
xlabel('\theta : free GTP to free GDP ratio', 'fontsize', 20);
legend('n_{RanGDP}','n_{RanGTP}')
ylabel('RanGD(T)P fraction');
set(gca,'Fontsize',18);


hFig3 = figure(3);
set(hFig3, 'Position', [100, 96, 1086, 454])
title('[NTF2]_{tot}=100 nM, [C]_{tot}=10 nM')
subplot(2,2,1);
plot(theta,speciesequvalue50N(:,12)./speciesequvalue50N(:,5),...
     theta,speciesequvalue50N(:,14)./speciesequvalue50N(:,7),...
     theta,speciesequvalue50N(:,11)./speciesequvalue50N(:,4),'LineWidth',2.5);
grid on;
xlabel('\theta', 'fontsize', 20);
legend('[C]_{nu}/[C]_{cyto}','[Im]_{nu}/[Im]_{cyto}','[Im-C]_{nu}/[Im-C]_{cyto}');
set(gca,'Fontsize',18);

    
subplot(2,2,3);
plot(theta,RanGDPfraction50N,theta,RanGTPfraction50N,'LineWidth',2.5);
grid on;
axis([1 100 0 1]);
xlabel('\theta : free GTP to free GDP ratio', 'fontsize', 20);
legend('n_{RanGDP}','n_{RanGTP}')
ylabel('RanGD(T)P fraction');
set(gca,'Fontsize',18);

subplot(2,2,2);
plot(theta,speciesequvalue50(:,12)./speciesequvalue50(:,5),...
     theta,speciesequvalue50(:,14)./speciesequvalue50(:,7),...
     theta,speciesequvalue50(:,11)./speciesequvalue50(:,4),'LineWidth',2.5);
grid on;
xlabel('\theta', 'fontsize', 20);
legend('[C]_{nu}/[C]_{cyto}','[Im]_{nu}/[Im]_{cyto}','[Im-C]_{nu}/[Im-C]_{cyto}');
set(gca,'Fontsize',18);

subplot(2,2,4);
plot(theta,RanGDPfraction50,theta,RanGTPfraction50,'LineWidth',2.5);
grid on;
axis([1 100 0 1]);
xlabel('\theta : free GTP to free GDP ratio', 'fontsize', 20);
legend('n_{RanGDP}','n_{RanGTP}')
ylabel('RanGD(T)P fraction');
set(gca,'Fontsize',18);



