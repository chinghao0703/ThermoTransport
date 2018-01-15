clear all;


%addpath data; 



trans_linear=2*(0:2:100);
N2_linear=2*(0:1:150);
TatwhichNCmaxN50=zeros(1,3);
TatwhichEPmaxN50=zeros(1,3);
TatwhichNCmaxN100=zeros(1,3);
TatwhichEPmaxN100=zeros(1,3);
TatwhichNCmaxN200=zeros(1,3);
TatwhichEPmaxN200=zeros(1,3);
TatwhichEPminN200=zeros(1,3);
cargo = [20,50,100];
t2= linspace(0,200);





% for c2, NTF2=50,100
load('globalphase_m2_dataset1_c2_vol5_NEWIC')
CargoNCmaxtmp2=maxcargoNCratio(:,26);
[M,I]= max(CargoNCmaxtmp2);
TatwhichNCmaxN50(1)= trans_linear(I);


%EP2=EPhat(:,26)';
EP2max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP2max = EP2max(:,26)';

[M,I]=max(EP2max);
TatwhichEPmaxN50(1)= trans_linear(I);

[M,I]=min(EP2max);
TatwhichEPminN50(1)= trans_linear(I);



CargoNCmaxtmp2=maxcargoNCratio(:,51);
[M,I]= max(CargoNCmaxtmp2);
TatwhichNCmaxN100(1)= trans_linear(I);


%EP2=EPhat(:,26)';
EP2max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP2max = EP2max(:,26)';

[M,I]=max(EP2max);
TatwhichEPmaxN100(1)= trans_linear(I);

[M,I]=min(EP2max);
TatwhichEPminN100(1)= trans_linear(I);




% for c2, NTF2=200
load('globalphase_m2_dataset2_c2_vol5_NEWIC')
CargoNCmaxtmp2=maxcargoNCratio(:,50);
[M,I]= max(CargoNCmaxtmp2);
TatwhichNCmaxN200(1)= trans_linear(I);


%EP2=EPhat(:,50)';
EP2max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP2max = EP2max(:,50)';

[M,I]=max(EP2max);
TatwhichEPmaxN200(1)= trans_linear(I);

[M,I]=min(EP2max);
TatwhichEPminN200(1)= trans_linear(I);




% for c5, NTF2=50,100
load('globalphase_m2_dataset1_c5_vol5_NEWIC')
CargoNCmaxtmp5=maxcargoNCratio(:,26);
[M,I]= max(CargoNCmaxtmp5);
TatwhichNCmaxN50(2)= trans_linear(I);


%EP5=EPhat(:,26)';
EP5max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP5max = EP5max(:,26)';

[M,I]=max(EP5max);
TatwhichEPmaxN50(2)= trans_linear(I);

[M,I]=min(EP5max);
TatwhichEPminN50(2)= trans_linear(I);



CargoNCmaxtmp5=maxcargoNCratio(:,51);
[M,I]= max(CargoNCmaxtmp5);
TatwhichNCmaxN100(2)= trans_linear(I);


%EP2=EPhat(:,26)';
EP5max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP5max = EP5max(:,26)';

[M,I]=max(EP5max);
TatwhichEPmaxN100(2)= trans_linear(I);

[M,I]=min(EP5max);
TatwhichEPminN100(2)= trans_linear(I);

% for c5, NTF2=200
load('globalphase_m2_dataset2_c5_vol5_NEWIC')
CargoNCmaxtmp5=maxcargoNCratio(:,50);
[M,I]= max(CargoNCmaxtmp5);
TatwhichNCmaxN200(2)= trans_linear(I);

%EP5=EPhat(:,50)';
EP5max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP5max = EP5max(:,50)';
[M,I]=max(EP5max);
TatwhichEPmaxN200(2)= trans_linear(I);

[M,I]=min(EP5max);
TatwhichEPminN200(2)= trans_linear(I);


% for c10, NTF2=50,100
load('globalphase_m2_dataset1_c10_vol5_NEWIC')
CargoNCmaxtmp10=maxcargoNCratio(:,26);
[M,I]= max(CargoNCmaxtmp10);
TatwhichNCmaxN50(3)= trans_linear(I);


%EP5=EPhat(:,26)';
EP10max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP10max = EP10max(:,26)';

[M,I]=max(EP10max);
TatwhichEPmaxN50(3)= trans_linear(I);

[M,I]=min(EP10max);
TatwhichEPminN50(3)= trans_linear(I);



CargoNCmaxtmp10=maxcargoNCratio(:,51);
[M,I]= max(CargoNCmaxtmp10);
TatwhichNCmaxN100(3)= trans_linear(I);


%EP2=EPhat(:,26)';
EP10max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP10max = EP10max(:,26)';

[M,I]=max(EP10max);
TatwhichEPmaxN100(3)= trans_linear(I);

[M,I]=min(EP10max);
TatwhichEPminN100(3)= trans_linear(I);

% for c10, NTF2=200
load('globalphase_m2_dataset2_c10_vol5_NEWIC')
CargoNCmaxtmp10=maxcargoNCratio(:,50);
[M,I]= max(CargoNCmaxtmp10);
TatwhichNCmaxN200(3)= trans_linear(I);


%EP10=EPhat(:,50)';
EP10max = reshape(max(EP,[],2),size(EP,1),size(EP,3));
EP10max = EP10max(:,50)';
[M,I]=max(EP10max);
TatwhichEPmaxN200(3)= trans_linear(I);

[M,I]=min(EP10max);
TatwhichEPminN200(3)= trans_linear(I);



%%%
[M,I]=max(EP2max(1:40));
TatwhichEPmaxN100(1)= trans_linear(I);
EPmaxN100(1) = M;

[M,I]=min(EP2max(1:40));
TatwhichEPminN100(1)= trans_linear(I);
EPminN100(1) = M;

[M,I]=max(EP5max(1:40));
TatwhichEPmaxN100(2)= trans_linear(I);
EPmaxN100(2) = M;

[M,I]=min(EP5max(1:40));
TatwhichEPminN100(2)= trans_linear(I);
EPminN100(2) = M;

[M,I]=max(EP10max(1:40));
TatwhichEPmaxN100(3)= trans_linear(I);
EPmaxN100(3) = M;

[M,I]=min(EP10max(1:40));
TatwhichEPminN100(3)= trans_linear(I);
EPminN100(3) = M;



%%%
NLmax(1) = max(CargoNCmaxtmp2);
NLmax(2) = max(CargoNCmaxtmp5);
NLmax(3) = max(CargoNCmaxtmp10);

hFig1= figure(1);
set(hFig1, 'Position', [100, 96, 1086, 454])
subplot(1,2,1)
plot(trans_linear,CargoNCmaxtmp2, trans_linear,CargoNCmaxtmp5, trans_linear,CargoNCmaxtmp10, 'LineWidth',2.5);
hold on;
plot(TatwhichNCmaxN100, NLmax, 'kd','MarkerSize',12)
grid on;
xlabel('[Importin] (nM) ', 'fontsize', 20);
ylabel('NL ratio');
legend('[C]_{tot}=20 nM','[C]_{tot}=50 nM','[C]_{tot}=100 nM')
set(gca,'Fontsize',18);

subplot(1,2,2)
%plot(trans_linear,EP2curve, trans_linear,EP5curve,trans_linear,EP10curve,trans_linear,EP2,'o', trans_linear,EP5,'o',trans_linear,EP10,'o', 'LineWidth',2.5);
plot(trans_linear,EP2max, trans_linear,EP5max,trans_linear,EP10max, 'LineWidth',2.5);
hold on;
plot(TatwhichEPminN100, EPminN100, 'ko', TatwhichEPmaxN100, EPmaxN100, 'ks','MarkerSize',12)
grid on;
xlabel('[Importin] (nM)', 'fontsize', 20);
ylabel('EP/k_BT');
legend('[C]_{tot}=20 nM','[C]_{tot}=50 nM','[C]_{tot}=100 nM')
set(gca,'Fontsize',18);







figure(2)
scatter(cargo, TatwhichNCmaxN100, 150,'filled', 'bd','LineWidth',9);
hold on;
scatter(cargo, TatwhichEPminN100, 150,'filled', 'ro','LineWidth',9); 
hold on;
scatter(cargo, TatwhichEPmaxN100, 150,'filled', 'gs','LineWidth',9);
grid on;
legend('NL max, [NTF2]=100', 'EP min, [NTF2]=100', 'EP max, [NTF2]=100');
xlabel('[C]', 'fontsize', 28);
ylabel('[Im^*]','fontsize', 28)
set(gca,'Fontsize',28);





