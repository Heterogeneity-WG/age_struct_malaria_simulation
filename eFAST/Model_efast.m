%% EFAST main function
% SI[POI,timepoints,QOI] : first order sensitivity indices
% STI[POI,timepoints,QOI] : total effect sensitivity indices
% CV[POI,timepoints,QOI]
clear;
close all;
%% eFAST configurations
NR = 5; %: no. of search curves - RESAMPLING
k = 8 + 1; % # of input factors (parameters varied) + dummy parameter
NS = 65; % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points
MI = 4; %: maximum number of fourier coefficients that may be retained in calculating the partial variances without interferences between the assigned frequencies

%% Model parameters and Numerical ODE setting
Parameter_settings_EFAST;

% Computation of the frequency for the group of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= 65 per factor.\n']);
    return;
end

%% Pre-allocation of the output matrix Y
% Y will save only the points of interest specified in the vector time_points
Y(NS,length(time_points),length(y0),length(pmin),NR)=0;  % pre-allocation

% Loop over k parameters (input factors)
for i=1:k % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);   
    % Loop over the NR search curves.
    for L=1:NR
        % Setting the vector of frequencies OM for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar variable S and the coordinates {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1        
        % Transform distributions from standard uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); % this is what assigns 'our' values rather than 0:1 dist
        
        % Do the NS model evaluations.
        for run_num=1:NS
            disp([num2str([i, run_num, L]), '//parameter run NR']) % keeps track of [parameter run NR]
            %% put in model time evolution HERE
            % ODE model
            f = @ODE_efast;
            [t,y]=ode15s(@(t,y)f(t,y,X(:,:,i,L),run_num),tspan,y0,[]); 
            %%
            Y(run_num,:,:,i,L)=y(time_points+1,:); % saves output at the time points of interest
        end 
    end 
end 
% save Model_efast.mat; % save the workspace

%% CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,1:4);
[CVsi, CVsti]=CVmethod(Si, rangeSi,Sti,rangeSti,4); % Calculate Coeff. of Var. for Si and STi for Viral load (variable 4). See online Supplement A.5 for details.
s_HIV = efast_ttest(Si,rangeSi,Sti,rangeSti,1:length(time_points),efast_var,4,y_var_label,0.05); % T-test on Si and STi for Viral load (variable 4)

