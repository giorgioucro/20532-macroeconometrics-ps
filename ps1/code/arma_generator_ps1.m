function [y,e,AR_polynomial,MA_polynomial] = arma_generator_ps1( ...
    T,sigma2_e,AR_input,MA_input,input_type,burn)
% ARMA_GENERATOR_PS1 Generate T observations from an ARMA(p,q) with loops.
%
% The function implements the lag-polynomial equation
%
%       A(L)y_t = B(L)e_t,
%
% where
%
%       A(L) = 1 + A_1 L + ... + A_p L^p,
%       B(L) = 1 + B_1 L + ... + B_q L^q,
%       e_t  is Gaussian white noise with variance sigma2_e.
%
% Solving the equation for y_t gives the recursion
%
%       y_t = -A_1 y_(t-1) - ... - A_p y_(t-p)
%             + e_t + B_1 e_(t-1) + ... + B_q e_(t-q).
%
% REQUIRED INPUTS
%   T          Number of observations returned after burn-in.
%   sigma2_e   Variance of the Gaussian white-noise innovation e_t.
%   AR_input   AR polynomial coefficients or AR polynomial roots.
%   MA_input   MA polynomial coefficients or MA polynomial roots.
%
% OPTIONAL INPUTS
%   input_type Use 'coefficients' when AR_input and MA_input already contain
%              polynomial coefficients.  Use 'roots' when they contain the
%              roots from which MATLAB should construct the polynomials.
%              The default is 'coefficients'.
%   burn       Number of initial observations discarded.  The default is
%              500, which reduces the effect of zero starting conditions.
%
% OUTPUTS
%   y             T-by-1 vector containing the ARMA realization.
%   e             T-by-1 vector containing the corresponding innovations.
%   AR_polynomial Row vector [1 A_1 ... A_p] actually used.
%   MA_polynomial Row vector [1 B_1 ... B_q] actually used.
%
% EXAMPLE USING POLYNOMIAL COEFFICIENTS
%   A(L)=1-0.8L and B(L)=1+0.4L+0.2L^2:
%
%   [y,e] = arma_generator_ps1(500,0.4,[1 -0.8],[1 0.4 0.2], ...
%                               'coefficients',500);
%
% EXAMPLE USING THE SAME POLYNOMIAL ROOTS
%   AR_roots = roots([1 -0.8]);
%   MA_roots = roots([1 0.4 0.2]);
%   [y,e] = arma_generator_ps1(500,0.4,AR_roots,MA_roots,'roots',500);

% If input_type is omitted, interpret the input vectors as coefficients.
if nargin<5 || isempty(input_type)
    input_type = 'coefficients';
end

% If burn is omitted, discard 500 initial observations by default.
if nargin<6 || isempty(burn)
    burn = 500;
end

% Check the scalar inputs before starting the simulation.
if T<=0 || T~=round(T)
    error('T must be a positive integer.');
end
if sigma2_e<0
    error('The innovation variance sigma2_e cannot be negative.');
end
if burn<0 || burn~=round(burn)
    error('burn must be a nonnegative integer.');
end

% Convert the two possible input formats into polynomial coefficients.
if strcmpi(input_type,'coefficients')
    % Empty input means that the corresponding polynomial has order zero.
    if isempty(AR_input)
        AR_polynomial = 1;
    else
        AR_polynomial = AR_input(:)';
    end

    if isempty(MA_input)
        MA_polynomial = 1;
    else
        MA_polynomial = MA_input(:)';
    end
elseif strcmpi(input_type,'roots')
    % poly constructs the coefficient vector of a polynomial from its roots.
    if isempty(AR_input)
        AR_polynomial = 1;
    else
        AR_polynomial = real(poly(AR_input));
    end

    if isempty(MA_input)
        MA_polynomial = 1;
    else
        MA_polynomial = real(poly(MA_input));
    end
else
    error('input_type must be either ''coefficients'' or ''roots''.');
end

% The leading coefficients must be nonzero.  Normalize both polynomials so
% that their first coefficient is one, as in the equation A(L)y=B(L)e.
if AR_polynomial(1)==0 || MA_polynomial(1)==0
    error('The leading coefficient of each polynomial must be nonzero.');
end
AR_polynomial = AR_polynomial/AR_polynomial(1);
MA_polynomial = MA_polynomial/MA_polynomial(1);

% p and q are the AR and MA orders.
p = length(AR_polynomial)-1;
q = length(MA_polynomial)-1;

% Simulate extra observations, which will be removed as burn-in.
total_observations = T+burn;
e_all = sqrt(sigma2_e)*randn(total_observations,1);
y_all = zeros(total_observations,1);

% Generate each observation from earlier y values and current/past shocks.
for t = 1:total_observations
    AR_part = 0;

    % Add -A_i*y_(t-i) whenever the lagged observation already exists.
    for i = 1:p
        if t-i>=1
            AR_part = AR_part-AR_polynomial(i+1)*y_all(t-i);
        end
    end

    MA_part = 0;

    % j=0 includes the current innovation e_t.  Larger j values add the
    % q lagged innovations according to the MA polynomial.
    for j = 0:q
        if t-j>=1
            MA_part = MA_part+MA_polynomial(j+1)*e_all(t-j);
        end
    end

    y_all(t) = AR_part+MA_part;
end

% Return only the final T observations and their contemporaneous shocks.
y = y_all(burn+1:end);
e = e_all(burn+1:end);

end
