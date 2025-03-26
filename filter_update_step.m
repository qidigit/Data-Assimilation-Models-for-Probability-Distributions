function [dz] = filter_update_step(zt,Ut,Rt, U_dot,R_dot,FF,SS, params)
MC=params.MC; %Dt=params.Dt; 
L=params.L; G=params.G; B=params.B; 
Bsis=eye(3);
sig_m=params.sig_m; sig_v=params.sig_v;
scale_m=params.scale_m; scale_v=params.scale_v;


zp=zt-mean(zt);

Hm = Buu(B,zp,zp);
Hm_bar = mean(Hm)'; %sum(Hm)'/(MC-1);
for i=1:3
    for j=1:3
        Lvs(i,j)=(Buu2(B,Ut,Bsis(:,j))+Buu2(B,Bsis(:,j),Ut)+(L-G)*Bsis(:,j))'*Bsis(:,i);
    end
end

% innovation & Kalman gains
hm=Buu2(B,Ut,Ut)-G*Ut+L*Ut + FF';
Inno_m = U_dot' - (Hm_bar'+hm');
Hm_p=Hm-mean(Hm);
zHmG=repmat(zp,[1 1 3]).*reshape(Hm_p*((sig_m)^(-2)),[MC 1 3]);
HmGHm=sum((Hm_p*((sig_m)^(-2))).*Hm_p,2);
Km_bar = mean(zHmG); %sum(zHmG)/(MC-1);

Hv=repmat(Hm,[1 1 3]).*reshape(zp,[MC 1 3]) + repmat(zp,[1 1 3]).*reshape(Hm,[MC 1 3]);
Hv_bar = mean(Hv); %sum(Hv)/(MC-1);
hv=Lvs*Rt+Rt'*Lvs'+SS*SS';
Inno_v = reshape(R_dot, [1 3 3]) - (Hv_bar+reshape(hv,[1 3 3]));
Hv_p=Hv-mean(Hv);
zHvG=repmat(zp,[1 1 3 3]).*reshape(Hv_p.*reshape(sig_v.^(-2),[1 3 3]),[MC 1 3 3]);
HvGHv=sum(sum((Hv_p.*reshape(sig_v.^(-2),[1 3 3])).*Hv_p,3),2);
Kv_bar=mean(zHvG); %sum(zHvG)/(MC-1);


dzHmG=(zHmG-Km_bar);
dzHvG=(zHvG-Kv_bar); 
% additional scaling to prevent divergence; may introduce additional error though
thd=1.*std(dzHmG); dzHmG=min(dzHmG,thd); dzHmG=max(dzHmG,-thd);  
thd=1.*std(dzHvG); dzHvG=min(dzHvG,thd); dzHvG=max(dzHvG,-thd);
dz = 1/2*sum(Km_bar.*reshape(Inno_m,[1 1 3]),3)+1/3*sum(sum(Kv_bar.*reshape(Inno_v,[1 1 3 3]),4),3) ...
     +scale_m*1/2*sum(dzHmG.*reshape(Inno_m,[1 1 3]),3)+scale_v*1/3*sum(sum(dzHvG.*reshape(Inno_v,[1 1 3 3]),4),3)  ...
     +1/2*sum(Km_bar.*reshape(Hm_bar,[1 1 3]),3)+1/3*sum(sum(Kv_bar.*reshape(Hv_bar,[1 1 3 3]),4),3) ...
     +1/4*mean(zp.*HmGHm)+1/9*mean(zp.*HvGHv);


