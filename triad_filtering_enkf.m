function [zt,TT,z_samps,Ts, mean_u,cov_u,mom_u] = triad_filtering_enkf(params, m0,var0, obs_m,obs_v)
MC=params.MC; Dt=params.Dt; Nt=params.Nt;
L=params.L; G=params.G; B=params.B; 
tstep=params.tstep; 
sig_m=diag(params.sig_m); sig_v=params.sig_v;

scale = .1;

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
zsamp=zt-mean(zt);
TT(ts)=0;
mean_u(:,ts)=Ut;
cov_u(:,:,ts)=Rt;
% save samples
ic=1; ns=floor(.5/Dt);
for ii=1:Nt
    
    t=Dt*(ii);     % time
    SS=params.SS+params.tmS(ii)^2*(params.SSt-params.SS);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % forecast step
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
    dz1=(L*zsamp')'-zsamp*G+Buu(B,Ut',zsamp)+Buu(B,zsamp,Ut') + fz-fm';
    dU1=Buu2(B,Ut,Ut)-G*Ut+L*Ut + Hm+params.tmM(ii,:)';
    dR1=Lvs*Rt+Rt'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(zsamp)-Rt);

    z1 = zsamp+.5*Dt*dz1;
    U1 = Ut   +.5*Dt*dU1;
    R1 = Rt   +.5*Dt*dR1;
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
    dz2=(L*z1')'-z1*G+Buu(B,U1',z1)+Buu(B,z1,U1') + fz-fm';
    dU2=Buu2(B,U1,U1)-G*U1+L*U1 + Hm+params.tmM(ii,:)';
    dR2=Lvs*R1+R1'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(z1)-R1);

    z2 = zsamp+.5*Dt*dz2;
    U2 = Ut   +.5*Dt*dU2;
    R2 = Rt   +.5*Dt*dR2;
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
    dz3=(L*z2')'-z2*G+Buu(B,U2',z2)+Buu(B,z2,U2') + fz-fm';
    dU3=Buu2(B,U2,U2)-G*U2+L*U2 + Hm+params.tmM(ii,:)';
    dR3=Lvs*R2+R2'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(z2)-R2);

    z3 = zsamp+1*Dt*dz3;
    U3 = Ut   +1*Dt*dU3;
    R3 = Rt   +1*Dt*dR3;
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
    dz4=(L*z3')'-z3*G+Buu(B,U3',z3)+Buu(B,z3,U3') + fz-fm';
    dU4=Buu2(B,U3,U3)-G*U3+L*U3 + Hm+params.tmM(ii,:)';
    dR4=Lvs*R3+R3'*Lvs'+SS^2+Q_f+Q_f' +params.eps*(cov(z3)-R3);
    
    dUt = dU1/6+dU2/3+dU3/3+dU4/6;
    dRt = dR1/6+dR2/3+dR3/3+dR4/6;

    % stoch samples
    zsamp=zsamp+Dt*(dz1/6+dz2/3+dz3/3+dz4/6);
    Wz = randn(MC,3);
    zsamp=zsamp+sqrt(Dt)*Wz*SS;
    zsamp=zsamp-mean(zsamp);
    % mean field
    Ut = Ut + Dt*dUt;
    % covariance
    Rt = Rt + Dt*dRt;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % analysis step
    if ~mod(ii,tstep)
        % observation
        dU_obs=(obs_m(:,ii+1)-obs_m(:,ii)) / (Dt*tstep); 
        dR_obs=(obs_v(:,:,ii+1)-obs_v(:,:,ii)) / (Dt*tstep);

        % Kalman gains
        zp=zsamp-mean(zsamp);
        Hm=Buu(B,zp,zp); 
        Hm_p=Hm-mean(Hm);
        %for mean
        C_zHm  = Dt*(zp')*(Hm_p)/(MC-1);
        C_HmHm = (Dt^2)*(Hm_p')*(Hm_p)/(MC-1);
        Sig_m2 = Dt*diag(sig_m)^2;
        Km = C_zHm/(C_HmHm+Sig_m2);
        %for covariance
        Hv=repmat(Hm,[1 1 3]).*reshape(zp,[MC 1 3]) + repmat(zp,[1 1 3]).*reshape(Hm,[MC 1 3]);
        Hv_p=Hv-mean(Hv);
        C_zHv = Dt*(zp')*(Hv_p(:,:))/(MC-1);
        C_HvHv = (Dt^2)*(Hv_p(:,:)')*(Hv_p(:,:))/(MC-1);
        Sig_v2 = Dt*diag(sig_v(:))^2;
        Kv = C_zHv/(C_HvHv+Sig_v2);

        % Innovations
        %for mean
        hm=Buu2(B,Ut,Ut)-G*Ut+L*Ut + params.tmM(ii,:)';
        % dI_m = Dt*(dU_obs - (Hm'+hm));
        dI_m = Dt*(dU_obs - dUt - scale*Hm_p');
        dz_m = Km*dI_m;
        %for covariance
        for i=1:3
            for j=1:3
                Lvs(i,j)=(Buu2(B,Ut,Bsis(:,j))+Buu2(B,Bsis(:,j),Ut)+(L-G)*Bsis(:,j))'*Bsis(:,i);
            end
        end
        hv=Lvs*Rt+Rt'*Lvs'+SS^2;
        dI_v = Dt*(dR_obs(:) - dRt(:) - scale*Hv_p(:,:)');
        dz_v = Kv*dI_v;

        % filtering samples 
        zsamp = zsamp + dz_m' +dz_v';
        % normalization for filter mean
        zsamp=zsamp-mean(zsamp);

    
        % save solution
        ts=ts+1;
        TT(ts)=t;

        zt=zsamp;
        
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

        display(['EnKF filter, iter = ',num2str(ii),': mean = ',num2str(mean_u(:,ts)')]);
    end
end