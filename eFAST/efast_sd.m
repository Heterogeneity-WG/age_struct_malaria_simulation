function [Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,time_points,var)
% Input
% Y - data matrix; OMi MI - eFast frequencies; 
% time_points - time to check the quantity (default = 1)
% var = index of QOIs to analyze (among Size_QOI) (default = 1)
% Output
% Si[POI,timepoints,QOI]
% Sti[POI,timepoints,QOI]

[NS, Size_timepts, Size_QOI, NP, NR] = size(Y);

for u=1:length(var)
    for t=1:length(time_points) %loop through time_points
        % disp([num2str([t u]), '//timepoint QOI']) 
        for i=1:NP % loop through parameters
            for L=1:NR 
                Y(:,t,u,i,L) = (Y(:,t,u,i,L)-mean(Y(:,t,u,i,L)))';
                % Fourier coeff. at [1:OMi/2].
                N=length(Y(:,t,u,i,L));
                NQ = (N-1)/2;
                N0 = NQ+1;
                COMPL = 0;
                Y_VECP = Y(N0+(1:NQ),t,u,i,L)+Y(N0-(1:NQ),t,u,i,L);
                Y_VECM = Y(N0+(1:NQ),t,u,i,L)-Y(N0-(1:NQ),t,u,i,L);
                for j=1:OMi/2
                    ANGLE = j*2*(1:NQ)*pi/N;
                    C_VEC = cos(ANGLE);
                    S_VEC = sin(ANGLE);
                    AC(j) = (Y(N0,t,u,i,L)+Y_VECP'*C_VEC')/N;
                    BC(j) = Y_VECM'*S_VEC'/N;
                    COMPL = COMPL+AC(j)^2+BC(j)^2;
                end
                % Computation of V_{(ci)}.
                Vci(L) = 2*COMPL;
                % Fourier coeff. at [P*OMi, for P=1:MI].
                COMPL = 0;
                Y_VECP = Y(N0+(1:NQ),t,u,i,L)+Y(N0-(1:NQ),t,u,i,L);
                Y_VECM = Y(N0+(1:NQ),t,u,i,L)-Y(N0-(1:NQ),t,u,i,L);
                for j=OMi:OMi:OMi*MI
                    ANGLE = j*2*(1:NQ)*pi/N;
                    C_VEC = cos(ANGLE');
                    S_VEC = sin(ANGLE');
                    AC(j) = (Y(N0,t,u,i,L)+Y_VECP'*C_VEC)/N;
                    BC(j) = Y_VECM'*S_VEC/N;
                    COMPL = COMPL+AC(j)^2+BC(j)^2;
                end
                % Computation of V_i.
                Vi(L) = 2*COMPL;
                % Computation of the total variance in the time domain.
                V(L) = Y(:,t,u,i,L)'*Y(:,t,u,i,L)/N;
            end %L
            % Computation of sensitivity indexes.
            Si(i,t,u) = mean(Vi)/mean(V);
            Sti(i,t,u) = 1-mean(Vci)/mean(V);
            rangeSi(i,t,:,u) = Vi./V;
            rangeSti(i,t,:,u) = 1-(Vci./V);
        end %i
    end %t
end