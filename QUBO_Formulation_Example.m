
vol=[6.26,5.26,5.16,16.2,12.3,23.8,7.30,6.30,6.20]';
mu=[0.350,0.250,0.240,4.09,2.50,6.02,1.08,0.980,0.970]';
v=[0.96,0.15,0.97,0.39,0.28,0.08,0.14,0.42,0.91]';

CorMat=[1	0.70	0.70	0.17	0.06	0.23	0.52	0.52	0.52
0.70	1	0.70	0.17	0.06	0.23	0.52	0.52	0.52
0.70	0.70	1	0.17	0.06	0.23	0.52	0.52	0.52
0.17	0.17	0.17	1	0.29	0.57	0.50	0.50	0.50
0.06	0.06	0.06	0.29	1	0.20	0.25	0.25	0.25
0.23	0.23	0.23	0.57	0.20	1	0.44	0.44	0.44
0.52	0.52	0.52	0.50	0.25	0.44	1	0.70	0.70
0.52	0.52	0.52	0.50	0.25	0.44	0.70	1	0.70
0.52	0.52	0.52	0.50	0.25	0.44	0.70	0.70	1];

Cov=CorMat.*(vol*vol');%covariance matrix

V=diag(v);

%Constraint 1
%use sum(mu'*V) to come up with sensible value for a
a=2;

%Constraint 2
%use sum(v) to come up with sensible value for v0
v0=3;

%Constraint 3
k=4;

n=length(vol);
N=15;%number of resulting qubits
P=100;

% % %Test Q matrix construction using the example in the tutorial paper (Glover 2018)
% % Q0a=diag([3,2,1,1,3,2]); 
% % Q1a=P*GetPartialMatrixFromSquaredPenaltyTerm([1;0;1;0;0;1;-1]);
% % Q2a=P*GetPartialMatrixFromSquaredPenaltyTerm([0;1;1;0;1;1;-1]);
% % Q3a=P*GetPartialMatrixFromSquaredPenaltyTerm([0;0;1;1;1;0;-1]);
% % Q4a=P*GetPartialMatrixFromSquaredPenaltyTerm([1;1;0;1;0;1;-1])
% % Qa=Q0a+Q1a+Q2a+Q3a+Q4a;
% % %Get Q matrix that matches that shown in the paper

%Now generate our Q matrix
Q0=zeros(N,N);
Q0(1:n,1:n)=V*Cov*V;

Q1=P*GetPartialMatrixFromSquaredPenaltyTerm([V*mu;-1;-2;0;0;0;0;-a]);
Q2=P*GetPartialMatrixFromSquaredPenaltyTerm([v;0;0;-1;-2;0;0;-v0]);
Q3=P*GetPartialMatrixFromSquaredPenaltyTerm([ones(9,1);0;0;0;0;-1;-2;-k]);

Q=Q0+Q1+Q2+Q3;%in symmetric form
Q_upper=zeros(size(Q))+2*triu(Q,1)+diag(diag(Q));% in upper triangular form if that is required

function Q=GetPartialMatrixFromSquaredPenaltyTerm(x)
    %Must pass in a column vector of coefficients equated to zero, where
    %the last element is the constant term
    temp=x*x';
    Q=temp(1:end-1,1:end-1)+2*diag(temp(end,1:end-1));
end

