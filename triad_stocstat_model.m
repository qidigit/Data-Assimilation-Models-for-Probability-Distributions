function [zt,TT,z_samps,Ts, mean_u,cov_u, mom_u] = triad_stocstat_model(params, m0, var0)
MC=params.MC; Dt=params.Dt; tstep=params.tstep; Nt=params.Nt;
L=params.L; G=params.G; B=params.B; 

% Initial conditions
zt(:,1)=sqrt(var0(1))*randn(MC,1);
zt(:,2)=sqrt(var0(2))*randn(MC,1);
zt(:,3)=sqrt(var0(3))*randn(MC,1);

% model coeffs
Bsis=eye(3);
for i=1:3
    for j=1:3
        coupling_strch(:,i+(j-1)*size(Bsis,2)) =  Buu2(B,Bsis(:,i),Bsis(:,j));
    end
end



ts=1;
Ut=m0';
Rt=diag(var0);
zsamp=zt;
TT(ts)=0;
mean_u(:,ts)=Ut;
cov_u(:,:,ts)=Rt;
% save samples
ic=1; ns=floor(.5/Dt);
for ii=1:Nt
    
    t=Dt*(ii);     % time

    % noise term
    SS=params.SS+params.tmS(ii)^2*(params.SSt-params.SS);
    Winc(:,1) = randn(MC,1);
    Winc(:,2) = randn(MC,1);
    Winc(:,3) = randn(MC,1);

    % RK4
    % step 1
    fz = Buu(B,zsamp,zsamp);
    fm = coupling_strch*Rt(:);
    Hm = mean(fz)';
    for i=1:3
        for j=1:3
            Lvs(i,j)=(Buu2(B,Ut,Bsis(:,j))+Buu2(B,Bsis(:,j),Ut)+(L-G)*Bsis(:,j))'*Bsis(:,i);
            Q_f(i,j)=mean(fz(:,i).*zsamp(:,j));
        end
    end
    k1=(L*zsamp')'-zsamp*G+Buu(B,Ut',zsamp)+Buu(B,zsamp,Ut') + fz-fm';
    l1=Buu2(B,Ut,Ut)-G*Ut+L*Ut + Hm+params.tmM(ii,:)';
    s1=Lvs*Rt+Rt'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(zsamp)-Rt);

    z1 = zsamp+.5*Dt*k1;
    U1 = Ut   +.5*Dt*l1;
    R1 = Rt   +.5*Dt*s1;
    % step 2
    fz = Buu(B,z1,z1);
    fm = coupling_strch*R1(:);
    Hm = mean(fz)';
    for i=1:3
        for j=1:3
            Lvs(i,j)=(Buu2(B,U1,Bsis(:,j))+Buu2(B,Bsis(:,j),U1)+(L-G)*Bsis(:,j))'*Bsis(:,i);
            Q_f(i,j)=mean(fz(:,i).*z1(:,j));
        end
    end
    k2=(L*z1')'-z1*G+Buu(B,U1',z1)+Buu(B,z1,U1') + fz-fm';
    l2=Buu2(B,U1,U1)-G*U1+L*U1 + Hm+params.tmM(ii,:)';
    s2=Lvs*R1+R1'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(z1)-R1);

    z2 = zsamp+.5*Dt*k2;
    U2 = Ut   +.5*Dt*l2;
    R2 = Rt   +.5*Dt*s2;
    % step 3
    fz = Buu(B,z2,z2);
    fm = coupling_strch*R2(:);
    Hm = mean(fz)';
    for i=1:3
        for j=1:3
            Lvs(i,j)=(Buu2(B,U2,Bsis(:,j))+Buu2(B,Bsis(:,j),U2)+(L-G)*Bsis(:,j))'*Bsis(:,i);
            Q_f(i,j)=mean(fz(:,i).*z2(:,j));
        end
    end
    k3=(L*z2')'-z2*G+Buu(B,U2',z2)+Buu(B,z2,U2') + fz-fm';
    l3=Buu2(B,U2,U2)-G*U2+L*U2 + Hm+params.tmM(ii,:)';
    s3=Lvs*R2+R2'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(z2)-R2);

    z3 = zsamp+1*Dt*k3;
    U3 = Ut   +1*Dt*l3;
    R3 = Rt   +1*Dt*s3;
    % step 4
    fz = Buu(B,z3,z3);
    fm = coupling_strch*R3(:);
    Hm = mean(fz)';
    for i=1:3
        for j=1:3
            Lvs(i,j)=(Buu2(B,U3,Bsis(:,j))+Buu2(B,Bsis(:,j),U3)+(L-G)*Bsis(:,j))'*Bsis(:,i);
            Q_f(i,j)=mean(fz(:,i).*z3(:,j));
        end
    end
    k4=(L*z3')'-z3*G+Buu(B,U3',z3)+Buu(B,z3,U3') + fz-fm';
    l4=Buu2(B,U3,U3)-G*U3+L*U3 + Hm+params.tmM(ii,:)';
    s4=Lvs*R3+R3'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(z3)-R3);

    dUt = l1/6+l2/3+l3/3+l4/6;
    dRt = s1/6+s2/3+s3/3+s4/6;

    % stoch samples
    zsamp=zsamp+Dt*(k1/6+k2/3+k3/3+k4/6);
    zsamp=zsamp+sqrt(Dt)*Winc*SS;
    zsamp=zsamp-mean(zsamp);
    % mean field
    Ut = Ut + Dt*dUt;
    % covariance
    Rt = Rt + Dt*dRt;
    

    
    % save solution
    if ~mod(ii,tstep)
        ts=ts+1;
        TT(ts)=t;

        zt(:,1)=zsamp(:,1);
        zt(:,2)=zsamp(:,2);
        zt(:,3)=zsamp(:,3);
        
        mean_u(:,ts)=Ut;
        cov_u(:,:,ts)=Rt;
        % get higher mom
        zp=zsamp-mean(zsamp);
        Hm=Buu(B,zp,zp);
        Hv=repmat(zp,[1 1 3]).*reshape(Hm,[MC 1 3]);
        mom_u(:,:,ts) = mean(Hv);
    end
    if mod(ii,ns)==0
        z_samps(:,:,ic)=zsamp;
        Ts(ic)=t;
        ic=ic+1;

        display(['Model iter = ',num2str(ii),': mean = ',num2str(mean_u(:,ts)')]);
    end
end