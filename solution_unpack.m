function  [SH, EH, DH, AH, VH, UH, SM, EM, IM, Cm, Cac, Cv, Ctot, MH] = solution_unpack(soln)
SH=soln.SH; EH=soln.EH; DH=soln.DH; AH=soln.AH(:,end); VH=soln.VH; UH=soln.UH;
SM=soln.SM; EM=soln.EM; IM=soln.IM;
Cm=soln.Cm; Cac=soln.Cac; Cv=soln.Cv; Ctot=soln.Ctot;
MH=soln.MH;
end
