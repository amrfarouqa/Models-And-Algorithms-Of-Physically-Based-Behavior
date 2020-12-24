function AmrMohEng2_Var1_Fin
close all; clc; 

llsk=3;
mass=[5, 1, 1, 1];   % masses
iner = [0.05, 0.05, 0.05, 0.05] * 0.05;    % moments of inertia
rad=[0.3, 0.5, 0.4, 0.2]; % radiuses
cor=[ 0 0 0; 2 0 0;  3 0 0; 4 0 0]; % initial coordinates (x,y,rotation)
nmz=length(mass); % number of nodes
NN=nmz*llsk;      % number of dofs 
            
g=9.81;           % acceleration of gravity              
F=zeros(1,nmz*llsk);F(2:llsk:end)=-mass*g;   % external forces
                         
stifp_n=[50000 50000 50000 50000];              % penalty stiffness (contacts)
dampp_n=[10 10 10 10]*5;                     % penalty normal dampers(contacts)
dampp_t=[1 1 1 1]*5;                         % penalty tangential dampers(contacts)
fric=[0.3 0.3 0.3 0.3 ];                        % coefficients of sliding friction

% boundary planes:
xmin=-2;xmax=10, ymin=-4; ymax=2;           % window boundaries
NRM=[0 1; 0 -1;1 0 ; -1 0;  1 2.5];           % line normals
PNT=[0 ymin; 0 ymax; xmin 0; xmax 0;  0 -1]; % point in line


U=zeros(NN,1);  DU=zeros(NN,1);  % initial displacements and velocities

% visualization:
figure(1); axis equal;axis ([xmin xmax ymin ymax]);grid on;hold on;
visualization(U,cor,rad, NRM,PNT);
pause
TT=20; dt=0.001;   % integration time and step
nsteps=TT/dt;   
t=0;               % initial time moment 
Urez=zeros(NN,1);  % array for resuls

for i=1:nsteps
    % calculate accerlerations at ith timestep:
    DDU=accelerations(U,DU,t,mass,stifp_n,dampp_n,dampp_t, fric,F,cor,rad, iner,NRM,PNT);
    DU=DU+dt*DDU';  % update velocities
    U=U+dt*DU;      % update displacements
    Urez(:,i+1)=U;  % save displacements 
    % visualization
    if(mod(i,10) ==0),
        cla; hold on
        visualization(U,cor,rad, NRM,PNT);
        hold off
        pause(0.01);
    end
    t=t+dt;
end
figure(2);hold on;
plot([0:dt:TT],Urez(2,:),'-b');
plot([0:dt:TT],Urez(5,:),'-r');
plot([0:dt:TT],Urez(8,:),'-g');
plot([0:dt:TT],Urez(11,:),'-m');
return
end

function DDU=accelerations(U,DU,t,mass,stifp_n,...
            dampp_n,dampp_t,fric,F,cor,rad, iner, NRM,PNT);
    %********************************************************************* 
    %*****  calculate accelerations  ******************
    %********************************************************************* 
    llsk=3;
    nmz=length(mass);NN=nmz*llsk;
    T=F;  % external forces 
    for i=1:nmz  % contact with boundaries
        r=[(i-1)*llsk+1:i*llsk];  du=DU(r); c=cor(i,:)'+U(r); 
        nconstr=size(NRM,1);
        for j=1:nconstr
            n=-NRM(j,:);n=n/norm(n); tau=[n(2), -n(1)]; % normal and tangential vectors of boundary plane
            A=PNT(j,:);                  % point in boundary
            dlt=dot(c(1:2)'-A,n)+rad(i);      % penetration of ith body to jth boundary;           
            if dlt > 0,
                rN= dlt*stifp_n(i)+dot(du(1:2),n)*dampp_n(i); if rN<0, rN=0; end 
                rT=(dot(du(1:2),tau)-du(3)*rad(i))*dampp_t(i); if abs(rT)>fric(i)*abs(rN), rT=sign(rT)*fric(i)*rN; end 
                T(r)=T(r)+[-rN*n-rT*tau,rT*rad(i)];
            end
        end
       
        for j=i+1:nmz % contact between rigid bodies
            s=[(j-1)*llsk+1:j*llsk];  duj=DU(s); cj=cor(j,:)'+U(s); 
            n=(cj(1:2)-c(1:2))/norm(cj(1:2)-c(1:2)); tau=[n(2);-n(1)]; %  direction vectors
            dlt=dot(c(1:2)-cj(1:2),n)+rad(i)+rad(j); % penetration during contact 

            if dlt > 0   % if bodies are in contact
                stifpn=min(stifp_n(i),stifp_n(j));
                damppn=min(dampp_n(i),dampp_n(j)); damppt=min(dampp_t(i),dampp_t(j));
                frc=min(fric(i),fric(j));
                rN= dlt*stifpn+dot(du(1:2),n)*damppn; if rN<0, rN=0; end
                rT=(dot(du(1:2)-duj(1:2),tau)-(du(3)*rad(i)-duj(3)*rad(j)))*damppt; if abs(rT)>frc*abs(rN), rT=sign(rT)*frc*rN; end 
                T(r)=T(r)+[-rN*n-rT*tau; rT*rad(i)]';
                T(s)=T(s)+[ rN*n+rT*tau;-rT*rad(j)]';
             end
        end
    end
   
    DDU(1:llsk:NN)=T(1:llsk:end)./mass;   % forces divided by masses
    DDU(2:llsk:NN)=T(2:llsk:end)./mass;   % forces divided by masses
    DDU(3:llsk:NN)=T(3:llsk:end)./iner;   % forces divided by inertia 
return
end

function visualization(U,cor,rad, NRM,PNT)
    %********************************************************
    % *****  visualization **********
    %********************************************************
    nmz=length(rad);llsk=3;
    for i=1:nmz
        % visualization :
        r=[(i-1)*llsk+1:i*llsk]; % dof numbers of ith rigid body
        u=U(r)+cor(i,:)'; % position of ith body
        spind=rad(i);  % radius of ith body
        % visualize a round body:
        rectangle('Position',[u(1)-spind,u(2)-spind,2*spind,2*spind],'Curvature',[1,1],'FaceColor',[0.4 0.6 1]);
        % line for orientation:
        plot([u(1),u(1)+spind*cos(u(3))], [u(2),u(2)+spind*sin(u(3))],'k-'); 
    end
    % visualize boundary lines:
for i=1:size(NRM,1)
    plot([PNT(i,1)-10*NRM(i,2),PNT(i,1)+10*NRM(i,2)],[PNT(i,2)+10*NRM(i,1),PNT(i,2)-10*NRM(i,1)],'k-') 
end
 
return
end