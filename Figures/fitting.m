function [Ps,invZResidues,invZConstant,invZProportional]=fitting(Fs,f,Npar,Niter,Ka,weight)
Ns   = length(f);               % Number of samples.
Npol = Npar;                    % Number of poles. Many poles as number of parallel branches Npar.
Ps   = InitialPoles(f,Npol);    % Set initial poles.
s    = 1i*2*pi*f.';              % Vector of values of variable "s"

for khg=1:Niter-1
    % Fit the trace of Yc (Poles)
    [Ps]=Poles(Fs,s,Ps,Ka,weight);
end

[invZResidues,invZConstant,invZProportional] = Residue(Fs,s,Ps,Ns,Ka,weight);

end


function [Ps]=InitialPoles(f,Npol)

MIN = min(f);
MAX = max(f);
if mod(Npol,2) == 0 % Even # of poles
    Ps = [];
else % Odd # of poles
    Ps = 0.5*(max(f)-min(f));
end

% Set the complex initial poles
beta = linspace(MIN,MAX,fix(Npol/2));
for n=1:length(beta)
    alfa=beta(n)*1e-2;
    Ps=[Ps alfa+1i*beta(n) alfa-1i*beta(n)];
end

Ps = Ps.'; % Column vector of initial poles
end


function [A]=Poles(Fs,s,Pi,Ka,weight)

% Fs : transfer function to be fitted
% s  : j*w
% Pi : poles derived from previous iteration.
% Ka : definies the structure of the function to be fitted.

Ns  = length(s);  % Length of vector containing frequencies
Np  = length(Pi); % Length of vector containing starting poles

CPX = imag(Pi)~=0; % 0 for real pole and 1 for complex pole
rp  = 0; % Initialize the index for real poles
cp  = 0; % Initialize the index for complex poles
RePole = []; % Initialize the vector of real poles
CxPole = [];%Initialize the vector of complex poles

% Loop to separate real poles and complex poles
for k = 1:Np
    if CPX(k) == 0     % Real Pole
        rp = rp + 1;
        RePole(rp) = Pi(k);
    elseif CPX(k) == 1 % Complex pole
        cp = cp + 1;
        CxPole(cp) = Pi(k);
    end
end

RePole = sort(RePole);    % Sort real poles
CxPole = sort(CxPole);    % Sort complex poles
Lambda = [RePole CxPole]; % Concatenate poles
I = eye(Np);              % Unit matrix
A = [];                   % Poles
B = weight.';               % the weight factor
C = [];                   % Residues
D = 0;                    % Constant term
E = 0;                    % Proportional term

num_real = numel(find(imag(Lambda) == 0));

dix      = [repelem(0,num_real) ...
            repmat([1 2],1,(Np-num_real)/2)]; % Initialized vector of pole types;

%%  Creates matrix A divided in four parts A = [A1 A2 A3 A4]
Dk=zeros(Ns,Np); % Initialize matrix with zeros
for m = 1 : Np % Iterative cycle for all poles
    if dix(m)== 0 % For a real pole
        Dk(:,m) = B./(s-Lambda(m));
    elseif dix(m)== 1 % For the real part
        Dk(:,m)=B./(s-Lambda(m)) + B./(s-Lambda(m)');
    elseif dix(m)== 2 % For the imag part
        Dk(:,m) = 1i.*B./(s-Lambda(m-1)) -1i.*B./(s-Lambda(m-1)');
    end
end

A1 = Dk;
A2 = ones(Ns,1);
A3 = s;

for col = 1: Np
    A4(:,col) = -(Dk(:,col).*Fs.');
end

% Asigns values to A
if Ka == 1
    A = [A1 A4]; % Strictly proper rational fitting
elseif Ka == 2
    A = [A1 A2 A4]; % Proper rational fitting
elseif Ka == 3
    A = [A1 A2 A3 A4]; % Improper rational fitting
else
    disp('Ka need to be 1, 2 or 3')
end

% Creates matrix b = B*Fs
b = B.*Fs.';

%% Solve the system Ax = b

% Separating real and imaginary part
% Following procedure is necessary in case of complex poles (see [1]).
An = [real(A); imag(A)]; % Real and imaginary part of A
bn = [real(b); imag(b)]; % Real and imaginary part of b


% Routine to apply the Euclidian norm to An
[Xmax Ymax] = size(An);
for col=1:Ymax
    Euclidian(col)=norm(An(:,col),2);
    An(:,col)=An(:,col)./Euclidian(col);
end

% Solving system
Xn = An\bn;
Xn = Xn./Euclidian.';

% Put the residues into matrix C
C = Xn(Ymax-Np+1:Ymax);

% C complex when the residues are complex
for m=1:Np
    if dix(m)==1
        alfa = C(m); % real part of a complex pole
        beta = C(m+1); % imag part of a complex pole
        C(m) = alfa + 1i*beta; % the complex pole
        C(m+1) = alfa - 1i*beta; % the conjugate
    end
end

% Now calculate the zeros for sigma
BDA = zeros(Np);
KQA = ones(Np,1);
% Loop to calculate the zeros of sigma which are the new poles
for km = 1:Np
    if dix(km)== 0 % For a real pole
        BDA(km,km) = Lambda(km);
    elseif dix(km)== 1 % For a cp with - imag part
        BDA(km,km) = real(Lambda(km));
        BDA(km,km+1) = imag(Lambda(km));
        KQA(km) = 2;
        Aux = C(km);
        C(km) = real(Aux);
    elseif dix(km)== 2 % For a cp with + imag part
        BDA(km,km) = real(Lambda(km));
        BDA(km,km-1) = imag(Lambda(km));
        KQA(km) = 0;
        C(km) = imag(Aux);
    end
end
ZEROS = BDA - KQA*C.';
POLS = eig(ZEROS).';
%Forcing (flipping) unstable poles to make them stable
uns = real(POLS) > 0;
POLS(uns) = POLS(uns)-2*real(POLS(uns));
% Sort poles in ascending order. First real poles and then complex poles
CPX = imag(POLS)~=0; % Set to 0 for a real pole and to1 for a
%complex pole
rp = 0; % Initialize index for real poles
cp = 0; % Initialize index for complex poles
RePole = []; % Initialize the vector of real poles
CxPole = []; % Initialize the vector of cp
% Loop to separate real and complex poles
for k = 1:Np
    if CPX(k) == 0 % Real Pole
        rp = rp + 1;
        RePole(rp) = POLS(k);
    elseif CPX(k) == 1 % Complex pole
        cp = cp + 1;
        CxPole(cp) = POLS(k);
    end
end
RePole = sort(RePole); % Sort real poles
CxPole = sort(CxPole); % Sort complex poles
% For conjugate pairs store first the one with positive imag part
CxPole = (CxPole.')';
NewPol = [RePole CxPole];
A = NewPol.'; % Output
end


%Function Residue

function [C,D,E]=Residue(Fs,s,Pi,Ns,Ka,weight)
Np = length(Pi);
CPX = imag(Pi)~=0; % 0 for a rp and 1 for cp
rp = 0; % Initialize the index for real poles
cp = 0; % Initialize the index for complex poles
RePole = []; % Initialize the vector of real poles
CxPole=[]; %Initialize the vector of complex poles
% Loop to separate real poles and complex poles
for k = 1:Np
    if CPX(k) == 0 % Real Pole
        rp = rp + 1;
        RePole(rp) = Pi(k);
    elseif CPX(k) == 1 % Complex pole
        cp = cp + 1;
        CxPole(cp) = Pi(k);
    end
end
RePole = sort(RePole); % Sort real poles
CxPole = sort(CxPole); % Sort complex poles
CxPole = (CxPole.')';
Lambda = [RePole CxPole];
I = eye(Np); % Unit diagonal matrix
A = []; % Poles
B = weight.';  % weight factor
C = []; % Residues
D = 0; % Constant term
E = 0; % Proportional term


num_real = numel(find(imag(Lambda) == 0));
dix = [repelem(0,num_real) ...
    repmat([1 2],1,(Np-num_real)/2)]; % Initialized vector of pole types;



% Output matrices:
k=zeros(Ns,Np);
for m=1:Np
    if dix(m)==0 % Real pole
        Dk(:,m) = B./(s-Lambda(m));
    elseif dix(m)==1 % Complex pole, 1st part
        Dk(:,m) = B./(s-Lambda(m)) + B./(s-Lambda(m)');
    elseif dix(m)==2 % Complex pole, 2nd part
        Dk(:,m) = i.*B./(s-Lambda(m-1)) - i.*B./(s-Lambda(m-1)');
    end
end

% Creates work space for matrices A and b
A1=Dk;
A2=B;
A3=B.*s;

if Ka == 1
    A = A1; % Strictly proper rational fit
elseif Ka == 2
    A = [A1 A2]; % Proper rational fit
elseif Ka == 3
    A = [A1 A2 A3]; % Improper fit
else
    disp('Ka must be 1, 2 or 3')
end

b = B.*Fs.';

An = [real(A); imag(A)]; % Real and imag part of A
bn = [real(b); imag(b)]; % Real and imag part of b

% Routine to apply the Euclidian norm to An
[Xmax,Ymax] = size(An);
for col=1:Ymax
    Eeuclidian(col)=norm(An(:,col),2);
    An(:,col)=An(:,col)./Eeuclidian(col);
end
% Solving system X
Xn=An\bn;
Xn=Xn./Eeuclidian.';

%% Deriving residues
C=Xn(1:Np);
% C is complex when the residues are complex
for m=1:Np
    if dix(m)==1
        alfa   = C(m); % real part of a complex pole
        beta   = C(m+1); % imag part of a complex pole
        C(m)   = alfa + 1i*beta; % the complex pole
        C(m+1) = alfa - 1i*beta; % the conjugate
    end
end

%% Outputs
if Ka == 1
    A = Lambda.'; % Poles
    C = C; % Residues
    D = 0; % Constant term
    E = 0; % Proportional term
elseif Ka == 2
    A = Lambda.'; % Poles
    C = C; % Residues
    D = Xn(Np+1); % Constant term
    E = 0; % Proportional term
elseif Ka == 3
    A = Lambda.'; % Poles
    C = C; % Residues
    D = Xn(Np+1); % Constant term
    E = Xn(Np+2); % Proportional term
end
end

