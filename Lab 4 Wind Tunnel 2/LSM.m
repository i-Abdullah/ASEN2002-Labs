function [ m_sl b_int sig_y sig_b sig_m Q ] = LSM(t,y)
% ---------------------------------------------------------------
% October 15th, 2018
% Done by :
%              1- Evan Shults
%              2- Abdulla Alameri
% ---------------------------------------------------------------
%
% This function is will try to provide a least square sum for a first
% degree polynomial using matrices theory. Part of a Project 1: Calorimetry
% for ASEN2012: Experimntal and Computatioanl Methods, CU Boulder, Fall
% 2018.
% ---------------------------------------------------------------
% INPUTS:
%           - t : span of dependent variable values
%           - y : span of independent variable values
%           - NOTE : t and y must be the same length and orientation
%           otherwise the function will throw an error.
%
% ---------------------------------------------------------------
% OUTPUTS:
%           - m_sl : best estimated slope.
%           - b_int : best y-intercept.
%           - sig_y : uncertainty in fitted y values.
%           - sig_b : uncertainty in y intercept.
%           - sig_m : uncertainty in the slope.
%           - Q : matrix with (1/sig_y^2) in the diagonal, used to find the
%           uncertitntiy in new x values that aren't included in the fit later .
%           Also used to find sig_b, sig_m.


% set up the matrices for operation
b_coef = ones(1,length(t))';
A = [ t b_coef ];
b = y;

%fidn the slope and intercept.

x = inv(transpose(A) * (A)) * transpose(A) * b;
m_sl = x(1);
b_int = x(2);

%% sigma y is equation  from the book [ ? (1/N-2) sum(yi-mti-b) ]


% the under square root
for i=1:length(b)
    
    sq(i) = (b(i) - (x(1)*A(i,1)) - x(2))^2;
    
end

% sig y 

sig_y = sqrt( (1/(length(b)-2)) * sum(sq));

% weight matrix 
weight_matrix = zeros(length(b),length(b));

%put the sig/y in the daiagonal 

for i = 1:length(b)
    
    weight_matrix(i,i) = 1 / (sig_y)^2;
    
end

%find Q.
Q = inv( transpose(A) * weight_matrix * A ) ;

%Extract sig_m, sig_b.

sig_m = sqrt(Q(1,1));
sig_b = sqrt(Q(2,2));


end