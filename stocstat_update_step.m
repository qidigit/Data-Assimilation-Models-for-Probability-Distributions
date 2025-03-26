function [dz,dU,dR] = stocstat_update_step(zt,Ut,Rt, FF,SS, params)
MC=params.MC; Dt=params.Dt; 
L=params.L; G=params.G; B=params.B; 
coupling_strch=params.coupling_strch; Bsis=eye(3);

zp=zt-mean(zt);
hm=Buu2(B,Ut,Ut)-G*Ut+L*Ut + FF';
Hm = Buu(B,zp,zp);
Hm_bar = mean(Hm)'; 
fm = coupling_strch*Rt(:);
for i=1:3
    for j=1:3
        Lvs(i,j)=(Buu2(B,Ut,Bsis(:,j))+Buu2(B,Bsis(:,j),Ut)+(L-G)*Bsis(:,j))'*Bsis(:,i);
    end
end
hv=Lvs*Rt+Rt'*Lvs'+SS*SS';
Hv=repmat(Hm,[1 1 3]).*reshape(zp,[MC 1 3]) + repmat(zp,[1 1 3]).*reshape(Hm,[MC 1 3]);
Hv_bar = mean(Hv); 

dz = (L*zt')'-zt*G+Buu(B,Ut',zt)+Buu(B,zt,Ut') + Hm-fm';
dU = hm + Hm_bar;
dR = hv+squeeze(Hv_bar) + params.eps*(cov(zt)-Rt);