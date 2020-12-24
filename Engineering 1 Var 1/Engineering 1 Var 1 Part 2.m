function AmrMohVar1Part2


% Newmark's method for 2D spring system with kinematic loads
% (large displacements)

clear all;
close all; clc;

k=[10 20 20 20 10 10 10];  % Xstiffness
m=[0.1 2 2 0.1 2];         % Xmasses
c=[2 10 10 10 2 2 2] ;          % Xdamping coefficients
F=-[0 m(1) 0 m(2) 0 m(3) 0 m(4) 0 m(5)]'*9.81; % Xforces
ind=[1 2 ;2 3 ; 2 4; 4 5; 5 3; 4 3;3 1];      % Xindices for elements
coordx = [1 3 3 5 5]; % XXcoordinates in X and Y (1D springs)
coordy = [2 1 4 1 3];

%IC=[0 0 0 0 1 1 0 0 0 0];   % Xpattern for kinematic loads
IC=[1 0 0 0 1 1 0 0 0 1];
% parameters of known displacements: ********
s=[0; -2;-1;0];  % Xdisplacement
%s2=[0; 0; 0 ;0];
%tz=1;          % start of kinematic load
omega=pi;
% *****************************************

nel=length(k), nmz=length(m)   % number of nodes and elements
M1(1:2:(2*nmz))= m; M1(2:2:(2*nmz))= m;
M=diag(M1);              % mass matrix

U=zeros(2*nmz,1);DU=zeros(2*nmz,1);DDU=zeros(2*nmz,1); % initial conditions

K=zeros(2*nmz);           % stiffness matrix
C=zeros(2*nmz);           % stiffness matrix
for j=1:nel
    i1=ind(j,1);i2=ind(j,2);
    K([2*i1-1, 2*i1,2*i2-1,2*i2],[2*i1-1, 2*i1,2*i2-1,2*i2]) =...
        K([2*i1-1, 2*i1,2*i2-1,2*i2],[2*i1-1, 2*i1,2*i2-1,2*i2])+...
        matrix_K(coordx([i1, i2]),coordy([i1, i2]),U([2*i1-1, 2*i1,2*i2-1,2*i2]),k(j));
end

% construct blocks of matrices:
M1=M(~IC,~IC);  f1=F(~IC); % main blocks
K12=K(~IC,find(IC)); C12=C(~IC,find(IC));
M12=M(~IC,find(IC));  % blocks respective to nodes which are loaded kinematically


TT=2; dt=0.01;
nsteps=TT/dt;
t=0.1;
Urez(1:2*nmz,1:nsteps+1)=0;                    % results

b0=dt^2/5; b1=dt/3; b2=3;       % parameters of Newmark's method
figure(1); hold on; axis([0 6 -2 5]); grid on;
h = [];
itmax = 1000;
epsmax = 1e-1;

for i=1:nsteps                  % numerical integration step
    dltDDU=zeros(2*nmz,1);
    q0=U(~IC)+dt*DU(~IC)+dt*dt/2*DDU(~IC);
    q1=DU(~IC)+dt*DDU(~IC);
    q2=DDU(~IC);
    % add forces generated due to the kinematic loads: .........
    ss=sin(omega*t);
    sv=omega*sin(omega*t);
    sa=-omega*omega*cos(omega*t);
    if t>=0.1 && t<=0.5
         f1k=f1-M12*s*sa-C12*s*sv-K12*s*ss;
        %  update known displacements, velocities and accelerations
        U(find(IC))=s*ss;
        DU(find(IC))=s*sv;
        DDU(find(IC))=s*sa;
    else, f1k=f1-K12*s*ss;
        U(find(IC))=s;
        DU(find(IC))=0;
        DDU(find(IC))=0;
    end
    % ...............................................................
    % Newton's Raphson's iterations
    for iii = 1:itmax
        
        K=zeros(2*nmz);           % stiffness matrix
        for j=1:nel
            i1=ind(j,1);i2=ind(j,2);
            K([2*i1-1, 2*i1,2*i2-1,2*i2],[2*i1-1, 2*i1,2*i2-1,2*i2]) =...
                K([2*i1-1, 2*i1,2*i2-1,2*i2],[2*i1-1, 2*i1,2*i2-1,2*i2])+...
                matrix_K(coordx([i1, i2]),coordy([i1, i2]),U([2*i1-1, 2*i1,2*i2-1,2*i2]),k(j));
        end       

         F1 = zeros(size(F));  % known forces
        F1(~IC) = f1k;

        psi = remainder(U,[coordx' coordy'],ind,k,F1, DU, c)-M*DDU;
       
        KD = b2*M+b1*C+b0*K;
        K12 = KD(~IC,find(IC));
        
        INCR=zeros(2*nmz,1);
        INCR(find(~IC))=KD(find(~IC),find(~IC))\(psi(find(~IC)));
        dltDDU = dltDDU + INCR;
        eps=norm(INCR)/(norm(dltDDU)+norm(INCR));
        U(~IC)=q0+b0*dltDDU(~IC);
        DU(~IC)=q1+b1*dltDDU(~IC);
        DDU(~IC)=q2+b2*dltDDU(~IC);
        if eps<epsmax, break,    end
        % visualization
        if ~isempty(h), delete(h); end
        h =  plotStructure(coordx+U(1:2:end)',coordy+U(2:2:end)',ind);
        fprintf(1,'time %10.6f s, iteration %d , tikslumas %10.6f\n',...
            t,iii,eps);
        title(sprintf('time %10.6f s, iteration %d ',t,iii));
        pause(0.00001);
    end
    % ...............................................................
    
    Urez(:,i+1)=U;
    
    
    t=t+dt;
    
end

figure(2);hold on;
plot([0:dt:TT],Urez(1,:),'--b');
plot([0:dt:TT],Urez(2,:),'-b');
plot([0:dt:TT],Urez(3,:),'--r');
plot([0:dt:TT],Urez(4,:),'-r');
plot([0:dt:TT],Urez(5,:),'--g');
plot([0:dt:TT],Urez(6,:),'-g');
plot([0:dt:TT],Urez(7,:),'--m');
plot([0:dt:TT],Urez(8,:),'-m');
return
end
function h =  plotStructure(coordx,coordy,Relations)
h = [];
for i = 1:size(Relations,1)
    h = [h; plot(coordx(Relations(i,:)),coordy(Relations(i,:)),'k-','LineWidth',2);];
end
h = [ h; plot(coordx, coordy, 'ko', 'MarkerFaceColor',[0 0 1], 'MarkerSize',12);];
end
function Ke=matrix_K(x,y,U,k)

L0=norm([x(2)-x(1),y(2)-y(1)]);
Lvec=[x(2); y(2)]+U(3:4)-...
    [x(1); y(1)]-U(1:2);
L=norm(Lvec); n=Lvec/L;

KTeblokas(1:2,1:2)=k*[ 1-L0/L+Lvec(1)^2*L0/L^3, Lvec(1)*Lvec(2)*L0/L^3; ...
    Lvec(1)*Lvec(2)*L0/L^3, (1-L0/L+Lvec(2)^2*L0/L^3) ];


Ke =  [KTeblokas, -KTeblokas; -KTeblokas, KTeblokas];

return
end
function psi=remainder(U,cor,ind,k,F,DU,c)

nel=size(ind,1);nmz=size(cor,1);
psi=zeros(2*nmz,1);

for i=1:nel
    
    mzr=ind(i,1); mzs=ind(i,2); r=[2*mzr-1,2*mzr]; s=[2*mzs-1,2*mzs];
    cr=cor(mzr,:);cs=cor(mzs,:);
    Lvec=cs+U(s)'-cr-U(r)'; L=norm(Lvec); n=Lvec/L;
    L0=norm(cs-cr);
    T=k(i)*(L-L0)+c(i)*dot(n, DU(s)-DU(r));
    psi([r,s])=psi([r,s])+[n'; -n']*T;
    
end

psi=psi+F;

return
end