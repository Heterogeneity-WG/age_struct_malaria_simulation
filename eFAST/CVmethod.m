function [CVsi, CVsti]=CVmethod(Si, rangeSi,Sti,rangeSti,out)
% Coeff. of Var.
meanSi=[];
meanSti=[];
[k, s, NR, u]=size(rangeSi);
if u==1
    out=1;
    for j=1:k
        for t=1:s
            meanSi(j,t,out)=(mean(rangeSi(j,t,:,out)));
            meanSti(j,t,out)=(mean(rangeSti(j,t,:,out)));
            stdSi(j,t,out)=(std(rangeSi(j,t,:,out)));
            stdSti(j,t,out)=(std(rangeSti(j,t,:,out)));
        end
    end
    a=Si(:,:,out)./meanSi(:,:,out);
    b=Sti(:,:,out)./meanSti(:,:,out);
else
    for j=1:k
        for t=1:s
            meanSi(j,t,out)=squeeze(mean(rangeSi(j,t,:,out)));
            meanSti(j,t,out)=squeeze(mean(rangeSti(j,t,:,out)));
            stdSi(j,t,out)=squeeze(std(rangeSi(j,t,:,out)));
            stdSti(j,t,out)=squeeze(std(rangeSti(j,t,:,out)));
        end
    end
a=stdSi(:,:,out)./meanSi(:,:,out);
b=stdSti(:,:,out)./meanSti(:,:,out);

end

CVsi=100*a;
CVsti=100*b;
end