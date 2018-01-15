
%clear all;
clc;

global Xdim; % number of species: 1~14
global tdim;
global k;    % all rate constants: k(i,1)=k_i^+, k(i,2)=k_i^-

vol=1; % 0: neglect volume factor
expFdiff=50;
Xdim=14; % 14 species
k=ones(11,2);
tdim=20;
tinterval=[0 tdim];
% X0: initial conditions
X0=[    50   50    0    0    5    0   50   50   25    0   0    5    0   50];
      % X1   X2   X3   X4   X5   X6   X7   X8   X9  X10  X11  X12  X13  X14

%index=0.04:.005:.25;
theta=0:1:120; % free [T]/[D]
RanDfraction=zeros(1,length(theta));
RanTfraction=zeros(1,length(theta));
Cargonufraction=zeros(1,length(theta));
Cargocytofraction=zeros(1,length(theta));
Transporternu2cyto=zeros(1,length(theta));
%TransportEffic=zeros(1,length(theta));
NTF2nufraction=zeros(1,length(theta));
NTF2cytofraction=zeros(1,length(theta));
TCnu2cyto=zeros(1,length(theta));
EP=zeros(1,length(theta));


for in2=1:length(theta)
    
    %theta2=1/20; % free parameter Eq.(41): ratio of free [D]/ free [T] at cyto
    %theta5=1/20; % free parameter Eq.(33): ratio of free [T]/ free [D] at nucl
    
    theta2=1/theta(in2);
    theta5=theta(in2);
    
   % below specify the kinetics rate constants (SM section II). Note that
   % cytosolic and nuclear volume effects will be taken care of in reaction2.m
   % We simulated cases both with and without volume factor. In fact, apart from the
   % change in the exact value of the steady-state nuclear localizatio measures (i.e. NL, EP etc.),
   % the qualitative behavior remains robust (i.e. same functional behavior of NL, EP
   % as one varies the importin, cargo, NTF2 concentrations)

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
    
    if vol==0
        [T,X] = ode45(@reaction,tinterval,X0);
    elseif vol==1
        [T,X] = ode45(@reaction2,tinterval,X0);
    end
        
    
    
    
    RanDfraction(in2)=X(length(T),2)/(X(length(T),2)+X(length(T),9)); % @SS
    RanTfraction(in2)=X(length(T),9)/(X(length(T),2)+X(length(T),9)); % @SS
    Cargonufraction(in2)=X(length(T),12)./(X(length(T),12)+X(length(T),5)); % @SS [C]_nu/[C]_cyto
    Cargocytofraction(in2)=X(length(T),5)./(X(length(T),12)+X(length(T),5));
    NTF2nufraction(in2)=X(length(T),8)/(X(length(T),8)+X(length(T),1)); 
    NTF2cytofraction(in2)=X(length(T),1)/(X(length(T),8)+X(length(T),1));
    TCnu2cyto(in2)=X(length(T),11)/X(length(T),4);
   
    
    Transporternu2cyto(in2)=X(length(T),14)/X(length(T),7); % time it takes to reach half max of [C]_nu
    
    normalization=sum(X(length(T),:));
    Xtemp=X(length(T),:);
    Pss=X(length(T),:)/normalization;
    
    
    Kdd=zeros(1,5);
    Kdd(2)=(k(2,1)/(k(2,1)+k(3,2)))/(k(2,2)*Xtemp(7)*Xtemp(2)/(k(2,2)*Xtemp(7)*Xtemp(2)+k(1,1)*Xtemp(1)*Xtemp(2)));
    Kdd(1)=(k(1,2)/(k(1,2)+k(6,1)))/(k(1,1)*Xtemp(1)*Xtemp(2)/(k(1,1)*Xtemp(1)*Xtemp(2)+k(2,2)*Xtemp(7)*Xtemp(2)));
    Kdd(4)=(k(4,2)/(k(4,2)+k(3,1)))/(k(4,1)*Xtemp(14)*Xtemp(9)/(k(4,1)*Xtemp(14)*Xtemp(9)+k(5,2)*Xtemp(8)*Xtemp(9)));
    Kdd(5)=(k(5,1)/(k(5,1)+k(6,2)))/(k(5,2)*Xtemp(8)*Xtemp(9)/(k(5,2)*Xtemp(8)*Xtemp(9)+k(4,1)*Xtemp(14)*Xtemp(9)));
    
    EP(in2)= Pss(6)*(k(2,1)/(k(2,1)+k(3,2)))*log(Kdd(2))...
        +Pss(6)*(k(3,2)/(k(2,1)+k(3,2)))*log((k(3,2)/(k(2,1)+k(3,2)))/(k(3,1)/(k(3,1)+k(4,2))))...
        + (Pss(2)+Pss(1))*(k(1,1)*Xtemp(1)*Xtemp(2)/(k(1,1)*Xtemp(1)*Xtemp(2)+k(2,2)*Xtemp(7)*Xtemp(2)))*log(1/Kdd(1))...
        +(Pss(2)+Pss(7))*(k(2,2)*Xtemp(7)*Xtemp(2)/(k(2,2)*Xtemp(7)*Xtemp(2)+k(1,1)*Xtemp(1)*Xtemp(2)))*log(1/Kdd(2))...
        +Pss(3)*(k(6,1)/(k(6,1)+k(1,2)))*log((k(6,1)/(k(6,1)+k(1,2)))/(k(6,2)/(k(6,2)+k(5,1))))...
        +Pss(3)*(k(1,2)/(k(1,2)+k(6,1)))*log(Kdd(1))...
        +Pss(10)*(k(5,1)/(k(5,1)+k(6,2)))*log(Kdd(5))...
        +Pss(10)*(k(6,2)/(k(5,1)+k(6,2)))*log((k(6,2)/(k(5,1)+k(6,2)))/(k(6,1)/(k(6,1)+k(1,2))))...
        +(Pss(14)+Pss(9))*(k(4,1)*Xtemp(14)*Xtemp(9)/(k(4,1)*Xtemp(14)*Xtemp(9)+k(5,2)*Xtemp(8)*Xtemp(9)))*log(1/Kdd(4))...
        +(Pss(8)+Pss(9))*(k(5,2)*Xtemp(8)*Xtemp(9)/(k(5,2)*Xtemp(8)*Xtemp(9)+k(4,1)*Xtemp(14)*Xtemp(9)))*log(1/Kdd(5))...
        +Pss(13)*(k(3,1)/(k(3,1)+k(4,2)))*log((k(3,1)/(k(3,1)+k(4,2)))/(k(3,2)/(k(3,2)+k(2,1))))...
        +Pss(13)*k(4,2)/(k(4,2)+k(3,1))*log(Kdd(4));
    
 
    
    
    
    
end


%[row,col] = find(RanD2RanT==max(max(RanD2RanT)));


figure(1);
subplot(221);
plot(theta,RanDfraction,theta,RanTfraction,'LineWidth',2.5);
axis([0 max(theta) 0 1]);
grid on;
xlabel('\theta', 'fontsize', 20);
legend('n_{RanGDP}','n_{RanGTP}')
title('fraction of RanGD(T)P');
set(gca,'Fontsize',18);

subplot(2,2,[3,4]);
plot(theta,Cargonufraction./Cargocytofraction,...
     theta,Transporternu2cyto,...
     theta,TCnu2cyto,'LineWidth',2.5);
grid on;
xlabel('\theta', 'fontsize', 20);
title('Nu to Cyto ratio');
legend('[C]_{nu}/[C]_{cyto}','[T]_{nu}/[T]_{cyto}','[TC]_{nu}/[TC]_{cyto}')
set(gca,'Fontsize',18);


subplot(222);
plot(theta,EP,'LineWidth',2.5);
axis([0 max(theta) min(EP) 1.03*max(EP)]);
grid on;
xlabel('\theta', 'fontsize', 20);
title('EP/k_B T');
set(gca,'Fontsize',18);

str=['expF' num2str(expF) 'theta' num2str(min(theta)) 'to' num2str(max(theta))];

