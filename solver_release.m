% This code solves the chemical kinetics equations and generates FIG.
% 3B,C in the main text. Other figures are generated using a parallelized
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
X0=[    50   50    0    0    5    0   50   50   25    0   0    5    0   50];
      % X1   X2   X3   X4   X5   X6   X7   X8   X9  X10  X11  X12  X13  X14


Transporter=(0:5:100);
freeCargo=5:5:50;
speciesequvalue=zeros(length(Transporter),14,length(freeCargo));


% theta = 10.1. The exact value of theta is not critical to our analysis,
% as shown in FIG. S4.
theta2=1/10.1;
theta5=10.1;

% exp(\Delta F) = exp(\Delta\tilde{F}) in Eq.(S35) (S43)
expFdiff= 50; 


% below specify the kinetics rate constants (SM section II). Note that the
% cytosolic and nuclear volume effects will be taken care of in reaction2.m
% We simulated cases both with and without volume factor. In fact, apart from the
% change in the exact value of the steady-state nuclear localizatio measures (NL, EP etc.),
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


% looping over variables and solve the coupled DE
for cargoint = 1: length(freeCargo)
    
    X0(5)  = freeCargo(cargoint);
    X0(12) = freeCargo(cargoint);
    
    for transint =1:length(Transporter)
        
        
        X0(7)  = Transporter(transint);
        X0(14) = Transporter(transint);
        
        
        [T,X] = ode45(@reaction2,tinterval,X0);
        
        speciesequvalue(transint,:,cargoint)=X(length(T),:);
        
    end
    
end



% now plot Figure 3 B,C:
hFig1 = figure(1);
set(hFig1, 'Position', [100, 96, 1086, 454])

subplot(121)
% note: k4^+ = 0.1
plot(Transporter,0.1*speciesequvalue(:,9,2),Transporter,0.1*speciesequvalue(:,9,5),...
     Transporter,0.1*speciesequvalue(:,9,8),Transporter,0.1*speciesequvalue(:,9,10),...
    'LineWidth',2.5, 'MarkerSize',8); 
grid on;
xlabel('X_{14}: [Im]_{nu} nM');
ylabel('scaled flux $\tilde{\Phi}_4^+$', 'Interpreter','latex')
h= legend('$[C]_{tot}=10 nM$','$[C]_{tot}=40 nM$','$[C]_{tot}=70 nM$','$[C]_{tot}=100 nM$');
set(h,'Interpreter','latex')
set(gca,'Fontsize', 20);

subplot(122)
plot(Transporter,speciesequvalue(:,12,2),Transporter,speciesequvalue(:,12,5),...
     Transporter,speciesequvalue(:,12,8),Transporter,speciesequvalue(:,12,10),...
    'LineWidth',2.5, 'MarkerSize',8); 
grid on;
xlabel('X_{14}: [Im]_{nu} nM');
ylabel('scaled flux $\tilde{\Phi}_7^-$', 'Interpreter','latex')
h= legend('$[C]_{tot}=10 nM$','$[C]_{tot}=40 nM$','$[C]_{tot}=70 nM$','$[C]_{tot}=100 nM$');
set(h,'Interpreter','latex')
set(gca,'Fontsize', 20);



% here find the maximum difference betwee \phi_4^+ and \phi_7^-
deltaflux=zeros(length(Transporter),4); 
deltaflux(:,1)= speciesequvalue(:,12,2)-0.1*speciesequvalue(:,9,2);
deltaflux(:,2)= speciesequvalue(:,12,5)-0.1*speciesequvalue(:,9,5);
deltaflux(:,3)= speciesequvalue(:,12,8)-0.1*speciesequvalue(:,9,8);
deltaflux(:,4)= speciesequvalue(:,12,10)-0.1*speciesequvalue(:,9,10);
[value,index]=max(deltaflux);
imstar = Transporter(index);



figure(2)

plot(Transporter,speciesequvalue(:,12,2)-0.1*speciesequvalue(:,9,2),Transporter,speciesequvalue(:,12,5)-0.1*speciesequvalue(:,9,5),...
     Transporter,speciesequvalue(:,12,8)-0.1*speciesequvalue(:,9,8),Transporter,speciesequvalue(:,12,10)-0.1*speciesequvalue(:,9,10),...
    'LineWidth',2.5, 'MarkerSize',8); 
hold on;
scatter(imstar,value,100,'filled','d');
grid on;
xlabel('X_{14}: [Im]_{nu} nM');
ylabel(' $\tilde{\Phi}_7^--\tilde{\Phi}_4^+$', 'Interpreter','latex')
h= legend('$[C]_{tot}=10 nM$','$[C]_{tot}=40 nM$','$[C]_{tot}=70 nM$','$[C]_{tot}=100 nM$');
set(h,'Interpreter','latex')
set(gca,'Fontsize', 20);
