function X = efast_gensamples(X,OMi,MI,pmin,pmax,type)

[NS,~,k,NR] = size(X);

for i=1:k % Loop over POIs, including dummy
    OMci = SETFREQ(k,OMi/2/MI,i);  % selecting the set of frequencies. OMci(i), i=1:k-1, contains the set of frequencies to be used by the complementary group.
    for L=1:NR  % Loop over the NR search curves
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
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,type); % this is what assigns 'our' values rather than 0:1 dist
    end 
end 
end