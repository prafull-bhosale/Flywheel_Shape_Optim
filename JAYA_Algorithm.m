% ------------------------ JAYA_Algorithm_code ----------------------------------------

function [xOpt,fOpt] = JAYA_Algorithm(objective,pop,var,x,maxGen,mini,maxi,FunctionTolerance,MaxStallGenerations,…
sp_params)
%% JAYA algorithm
%% Problem Definition

% pop = 1000;               % Population size
% var = 8;                 % Number of design variables
% maxGen = 20;            % Maximum number of iterations
% mini = 0.010*ones(1,var);  % Lower Bound of Variables
% maxi = 0.060*ones(1,var);   % Upper Bound of Variables
% objective = @myobj;      % Cost Function

%% initialize
[row,var] = size(mini);
x = zeros(pop,var);
fnew = zeros(pop,1);
f = zeros(pop,1);
fopt= zeros(1,1);
xopt=zeros(1,var);

%%  Generation and Initialize the positions
for i=1:var
    x(:,i) = mini(i)+(maxi(i)-mini(i))*rand(pop,1);
end
for i=1:pop
    f(i) = objective(x(i,:),sp_params);
end

%%  Main Loop
gen=1;
stallGenCount=0;
while(gen <= maxGen)
    
    [row,col]=size(x);
    [b,bindex]=min(f);
    Best=x(bindex,:);
    [w,windex]=max(f);
    worst=x(windex,:);
    xnew=zeros(row,col);
    
    for i=1:row
        for j=1:col
            xnew(i,j)=(x(i,j))+rand*(Best(j)-abs(x(i,j))) - (worst(j)-abs(x(i,j)));  %
        end
    end
    
    for i=1:row
        xnew(i,:) = max(min(xnew(i,:),maxi),mini);
        fnew(i,:) = objective(xnew(i,:),sp_params);
    end
    
    for i=1:pop
        if(fnew(i)<f(i))
            x(i,:) = xnew(i,:);
            f(i) = fnew(i);
        end
    end
    
    fnew = []; xnew = [];
    [fopt(gen),ind] = min(f);
    xopt(gen,:)= x(ind,:);
    if stallGenCount>=MaxStallGenerations 
        break;
    end
    
    if (gen>1) && ((abs(fopt(gen)-fopt(gen-1))/abs(fopt(gen))) < FunctionTolerance)
        stallGenCount=stallGenCount+1;                    
    else
        stallGenCount=0;
    end
    
    disp(['Iteration No. = ',num2str(gen), ',   Best Cost = ',num2str(min(f))])
    gen = gen+1;
end

%%

[fOpt,ind] = min(fopt);
xOpt=xopt(ind,:);
Fes = pop*ind;
disp(['Optimum value = ',num2str(fOpt,10)])
disp('Design variable = ');
xOpt

figure(1)
plot(fopt,'linewid',2)
xlabel('Iteration')
ylabel('Best Cost');
legend('JAYA')
disp(' ' )
end
