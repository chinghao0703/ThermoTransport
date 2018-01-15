function dX = reaction2(t,X)

global Xdim; % number of species: 1~14
global tdim; 
global k;    % all rate constants: k(i,1)=k_i^+, k(i,2)=k_i^-

volfactor=5; % i.e. \nu_{C}/\nu_{N} = 5 or \nu_{C} =500 um^3 and \nu_{N} =100 um^3, a= 100 um^3/sec;



dX = zeros(Xdim,1);% a column vector, dX(i,:) : time series of species i



dX(1) = -k(1,1)*X(1)*X(2)+k(1,2)*X(3)+k(11,1)*X(8)/volfactor-k(11,2)*X(1)/volfactor;

dX(2) = -k(1,1)*X(1)*X(2)+k(1,2)*X(3)+k(2,1)*X(6)-k(2,2)*X(7)*X(2);

dX(3) = k(1,1)*X(1)*X(2)-k(1,2)*X(3)-k(6,1)*X(3)/volfactor+k(6,2)*X(10)/volfactor;

dX(4) = -k(9,1)*X(4)/volfactor+k(9,2)*X(11)/volfactor+k(10,1)*X(5)*X(7)-k(10,2)*X(4);

dX(5) = -k(10,1)*X(5)*X(7)+k(10,2)*X(4);

dX(6) = -k(2,1)*X(6)+k(2,2)*X(2)*X(7)+k(3,1)*X(13)/volfactor-k(3,2)*X(6)/volfactor;

dX(7) = -k(2,2)*X(2)*X(7)+k(2,1)*X(6)+k(8,2)*X(14)/volfactor-k(8,1)*X(7)/volfactor ...
        +k(10,2)*X(4)-k(10,1)*X(7)*X(5);

dX(8) = -k(5,2)*X(8)*X(9)+k(5,1)*X(10)+k(11,2)*X(1)-k(11,1)*X(8);

dX(9) = -k(5,2)*X(8)*X(9)+k(5,1)*X(10)+k(4,2)*X(13)-k(4,1)*X(14)*X(9);

dX(10) = k(5,2)*X(8)*X(9)-k(5,1)*X(10)-k(6,2)*X(10)+k(6,1)*X(3);

dX(11) = k(7,2)*X(12)*X(14)-k(7,1)*X(11)-k(9,2)*X(11)+k(9,1)*X(4);

dX(12) = -k(7,2)*X(12)*X(14)+k(7,1)*X(11);

dX(13) = k(3,2)*X(6)-k(3,1)*X(13)-k(4,2)*X(13)+k(4,1)*X(14)*X(9);

dX(14) = k(4,2)*X(13)-k(4,1)*X(14)*X(9)-k(7,2)*X(14)*X(12)+k(7,1)*X(11)...
         -k(8,2)*X(14)+k(8,1)*X(7);
