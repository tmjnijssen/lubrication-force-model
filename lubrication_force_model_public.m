%{
Script companion to:

A note on the modelling of lubrication forces in unresolved simulations
Tim M.J. Nijssen, Marcel Ottens, Johan T. Padding

Powder Technology, 2022

contact: t.m.j.nijssen@tudelft.nl (Tim M.J. Nijssen)

Tested with: MATLAB R2021b and GNU Octave 6.4.0

Delft University of Technology (NL), 24-03-2022
%}

%% user-defined settings
% material properties (select preset below, or define custom)

% steel ball bearing
% dp     = [12.7e-3 12.7e-3]; % diameter, m
% rhop   = [7780 7780];       % density, kg/m3
% Y      = [190e9 190e9];     % Young's modules, Pa
% nu     = [0.27 0.27];       % Poisson's ratio, -
% edry   = 0.97;              % dry restitution coefficient, -
% sigma  = [24e-9 24e-9];     % surface roughness, m
% mufDry = 0.11;              % dry friction coefficient, -
% mufLub = 0.02;              % lubricated friction coefficient, -

% steel bead
% dp     = [12.7e-3 12.7e-3]; % diameter, m
% rhop   = [7780 7780];       % density, kg/m3
% Y      = [190e9 190e9];     % % Young's modules, Pa
% nu     = [0.27 0.27];       % Poisson's ratio, -
% edry   = 0.97;              % dry restitution coefficient, -
% sigma  = [272e-9 272e-9];   % surface roughness, m
% mufDry = 0.11;              % dry friction coefficient, -
% mufLub = 0.02;              % lubricated friction coefficient, -

% glass sphere
% dp     = [12.7e-3 12.7e-3]; % diameter, m
% rhop   = [2540 2540];       % density, kg/m3
% Y      = [60e9 60e9];       % Young's modules, Pa
% nu     = [0.23 0.23];       % Poisson's ratio, -
% edry   = 0.97;              % dry restitution coefficient, -
% sigma  = [134e-9 134e-9];   % surface roughness, m
% mufDry = 0.4;               % dry friction coefficient, -
% mufLub = 0.1;               % lubricated friction coefficient, -

% delrin sphere
dp     = [12.7e-3 12.7e-3];   % diameter, m
rhop   = [1400 1400];         % density, kg/m3
Y      = [2.8e9 2.8e9];       % Young's modules, Pa
nu     = [0.35 0.35];         % Poisson's ratio, -
edry   = 0.97;                % dry restitution coefficient, -
sigma  = [796e-9 796e-9];     % surface roughness, m
mufDry = 0.2;                 % dry friction coefficient, -
mufLub = 0.1;                 % lubricated friction coefficient, -

% fluid properties
rhof = 1e3;                   % density, kg/m3
etaf = 1e-3;                  % dynamic viscosity, Pa s
lubrication = true;           % enable lubrication model

% initial conditions
St     = 50;                  % total Stokes number, -
theta0 = deg2rad(20);         % angle of impact, rad

%% derived properties
% particle properties
np = 2;                       % number of particles, -
dp_ = prod(dp)/sum(dp);       % reduced diameter, m
rp  = dp/2;                   % radius, m
rp_ = dp_/2;                  % reduced radius, m
vp  = pi/6*dp.^3;             % volume, m3
Mp  = vp.*rhop;               % mass, kg
Mp_ = prod(Mp)/sum(Mp);       % reduced mass, kg

% see: https://www.cfdem.com/media/DEM/docu/gran_model_hertz.html
G_ = 1/(2*(2-nu(1))*(1+nu(1))/Y(1) + 2*(2-nu(2))*(1+nu(2))/Y(2));   % reduced shear modulus, Pa
Y_ = 1/((1-nu(1)^2)/Y(1) + (1-nu(2)^2)/Y(2));                       % reduced Young's modulus, Pa
beta = log(edry)/sqrt(log(edry)^2+pi^2);                            % -

% initial conditions
vr0    = St*6*pi*etaf*rp_^2/Mp_; % initial speed, Eq. 2, m/s
h0     = dp_;                    % initial gap height, m

phi0   = [0 0]; % initial orientation, rad
omega0 = [0 0]; % innitial rotational velocity, rad/s

x0  = [[-cos(theta0)*(h0)-rp(1) rp(2)];
       [-sin(theta0)*(h0)       0    ]]; % initial position, m
v0  = [[vr0*cos(theta0)         0    ];
       [vr0*sin(theta0)         0    ]]; % initial velocity, m/s

vrn0 = vr0*cos(theta0); % initial normal speed, m/s
vrt0 = vr0*sin(theta0); % initial tangential speed, m/s

% lubrication properties
hco = rp_;              % cut-off distance, m
Stn = St*cos(theta0);   % normal Stokes number, -
Stt = St*sin(theta0);   % tangential Stokes number, -

% solver settings
tc   = 2.87*(Mp_^2/(rp_*Y_^2*vr0))^(1/5);   % contact time, s
dt   = tc/50;                               % time step, s
OoM  = 10^floor(log10(dt));                 % order of magnitude of dt
dt   = round(dt/OoM)*OoM;                   % round dt to 1 significant digit
tend = 2*(h0)/vr0;                          % end time, s
nt   = ceil((tend+1e-12)/dt);               % numer of time steps, -

t    = (0:nt-1)*dt;                         % time values, s
plotEvery = floor(nt/200);                  % plot interval, -

%% minimum approach
hmins = mean(sigma); % minimum approach distance due to roughness, Eq. 3, m

% see: http://dx.doi.org/10.1098/rsta.2008.0014
hmine = 0.37*((etaf*vrn0/Y_).^2 .* rp_.^3).^(1/5);   % minimum approach distance due to deformation, Eq. 5, m

hmin = max(hmine,hmins); % minimum approach distance, Eq. 6, m

% effective friction coefficient, Eq. 10, -
if hmins > hmine || ~lubrication
    mufEff = mufDry;
else
    mufEff = mufLub;
end

%% preallocate 
F  = zeros(2,np,nt);    % total force, N

Fcn = zeros(nt,1);      % normal contact force magnitude, N
Fct = zeros(nt,1);      % tangential contact force magnitude, N
Fln = zeros(nt,1);      % normal lubrication force magnitude, N
Fn  = zeros(nt,1);      % total normal force magnitude, N

v = zeros(2,np,nt);     % velocity, m/s
v(:,:,1) = v0;

x = zeros(2,np,nt);     % position, m
x(:,:,1) = x0;

phi = zeros(1,np,nt);   % orientation, rad
phi(:,:,1) = phi0;

omega = zeros(1,np,nt); % rotational velocity, rad/s
omega(:,:,1) = omega0;

tau = zeros(1,np,nt);   % torque, Nm

Epot = zeros(nt,1);     % potential energy, J
deltat = zeros(nt,1);   % tangential overlap, m

%% open figure
figure(1)
clf
set(gcf,'Color','w')
hold on
axis equal
xlim([x0(1,1)-1.5*rp(1) x0(1,2)+1.5*rp(2)]);
ylim([-1 1]*(abs(x0(2,1))+1.5*max(rp)));
xticks(0)
yticks(0)
xticklabels({'x'})
yticklabels({'y'})
box on
set(gca,'FontSize',11)
set(gca,'LineWidth',0.75)
set(gca,'TickLength',[0 0])

%% main loop
for i = 2:nt
    % velocity verlet integration, position
    x(:,:,i)   =   x(:,:,i-1) +     v(:,:,i-1)*dt +   F(:,:,i-1)./(Mp           )*dt^2/2;
    phi(:,:,i) = phi(:,:,i-1) + omega(:,:,i-1)*dt + tau(:,:,i-1)./(2/5*Mp.*rp.^2)*dt^2/2;

    rij = x(:,2,i)-x(:,1,i);        % distance vector, m
    nij = rij/vecnorm(rij);         % normal unit vector, -
    tij = [-nij(2); nij(1)];        % tangential unit vector, -
    h   = vecnorm(rij)-sum(rp);     % gap height, m

    vr  = v(:,1,i-1) - v(:,2,i-1);  % relative velocity, m/s
    vrn = dot(nij,vr);              % normal relative velocity, m/s
    vrt = dot(tij,vr) + omega(1,1,i-1)*rp(1) + omega(1,2,i-1)*rp(2); % tangential relative velocity, m/s

    % overlap, Eq. 7, m
    if lubrication
        deltan = max(0, hmin-h);    % overlap, m
    else
        deltan = max(0, -h);        
    end

    % contact force
    if deltan>0
        % normal contact force, see: https://www.cfdem.com/media/DEM/docu/gran_model_hertz.html
        sn = 2*Y_*sqrt(rp_*deltan);                 % Pa m = N/m
        kn = 4/3*Y_*sqrt(rp_*deltan);               % Pa m = N/m
        gamman = -2*sqrt(5/6)*beta*sqrt(sn*Mp_);    % sqrt(Pa m kg) = kg/s
        Fcn(i) = -(kn*deltan + gamman*vrn);         % N/m*m = kg/s*m/s = N

        Epot(i) = 1/2*kn*deltan^2;
    else
        Fcn(i) = 0;
    end

    if h<0
        % tangential contact force, see: https://www.cfdem.com/media/DEM/docu/gran_model_hertz.html
        %                           and: https://www.cfdem.com/media/DEM/docu/gran_tangential_history.html
        deltat(i) = deltat(i-1) + vrt*dt;
        st = 8*G_*sqrt(rp_*deltan);                 % Pa m = N/m
        kt = 8*G_*sqrt(rp_*deltan);                 % Pa m = N/m
        gammat = -2*sqrt(5/6)*beta*sqrt(st*Mp_);    % sqrt(Pa m kg) = kg/s
        Fct(i) = -kt*deltat(i);
        if abs(Fct(i)) > mufEff*abs(Fcn(i))
            Fct(i) = mufEff*abs(Fcn(i))*sign(Fct(i));
            deltat(i) = abs(Fct(i))/kt*sign(deltat(i));
        else
            Fct(i) = Fct(i) - gammat*vrt;
        end

        Epot(i) = Epot(i) + 1/2*kt*deltat(i)^2;
    else
        Fct(i) = 0;
    end

    % lubrication force
    if lubrication
        hlub = max(h,hmin);
        Fln(i) = -6*pi*etaf*vrn*rp_^2/hlub; % lubrication force, Eq. 8, N
    end

    % total normal force, Eq. 9, N
    if lubrication
        if h>hco
            Fn(i) = 0;
        elseif h<=hco && h>hmin
            Fn(i) = Fln(i);
        elseif h<=hmin && h>0
            Fn(i) = (1 - h/hmin) * Fcn(i) + h/hmin * Fln(i);
        elseif h<=0
            Fn(i) = Fcn(i);
        end
    else
        Fn(i) = Fcn(i);
    end

    % apply force and torque
    F(:,1,i) =  Fn(i)*nij +  Fct(i)*tij;
    F(:,2,i) = -Fn(i)*nij + -Fct(i)*tij;
    tau(1,1,i) = Fct(i) * rp(1);
    tau(1,2,i) = Fct(i) * rp(2);

    % velocity verlet integration, velocity
    v(:,:,i)     =     v(:,:,i-1) + (  F(:,:,i-1) +   F(:,:,i))/2./(Mp           )*dt;
    omega(:,:,i) = omega(:,:,i-1) + (tau(:,:,i-1) + tau(:,:,i))/2./(2/5*Mp.*rp.^2)*dt;

    % draw spheres
    if mod(i,plotEvery)==0 || i==nt
        cla
        plot([x0(1,1)-1.5*rp(1) x0(1,2)+1.5*rp(2)],[0 0],'Color',0.7*ones(1,3),'LineWidth',0.75); % x axis
        plot([0 0], [-1 1]*(abs(x0(2,1))+1.5*max(rp)),   'Color',0.7*ones(1,3),'LineWidth',0.75); % y axis
        for j = 1:2
            rectangle('Position',[x(1,j,i)-rp(j) x(2,j,i)-rp(j) dp(j) dp(j)],'Curvature',[1 1],'FaceColor','white','LineWidth',0.75) % particle
            plot(x(1,j,i),x(2,j,i),'k+') % centre
            plot(squeeze(x(1,j,1:i)),squeeze(x(2,j,1:i)),'k-') % trajectory
            for k = (0:3)*pi/2
                plot(x(1,j,i)+cos(phi(1,j,i)+k)*rp(j)/2, x(2,j,i)+sin(phi(1,j,i)+k)*rp(j)/2,'k.') % orientation markers
            end    
        end
        % xCoM = squeeze(sum(Mp.*x, 2)/sum(Mp)); % centre of mass position, m
        % plot(xCoM(1,  i),xCoM(2,  i),'kx') % centre of mass
        % plot(xCoM(1,1:i),xCoM(2,1:i),'k-') % centre of mass trajectory
        text(0.95*x0(1,1)-1.5*rp(1),-0.9*(abs(x0(2,1))+1.5*max(rp)),['t* = ' num2str(t(i)/(rp(1)/vr0),'%03.2f')],'FontWeight','normal','FontSize',11)
        set(gca,'Layer', 'Top')
        drawnow
    end
end

%% output variables
h  = squeeze(vecnorm(x(:,2,:) - x(:,1,:),1,1))-sum(rp);  % gap height, m
vr = squeeze(vecnorm(v(:,1,:) - v(:,2,:)));              % relative velocity, m/s

xCoM = sum(Mp.*x, 2)/sum(Mp);  % contre of mass position, m
vCoM = sum(Mp.*v, 2)/sum(Mp);  % centre of mass velocity, m/s

ECoM    = 1/2*(sum(Mp)).*squeeze(vecnorm(vCoM)) .^2;     % centre of mass kinetic energy, J
Ekin    = 1/2*(Mp).*squeeze(vecnorm(v-repmat(vCoM,[1 2 1])))'.^2; % particle kinetic energy relative to centre of mass, J  
Erot    = 1/2*(2/5*Mp.*rp.^2).*squeeze(omega)'.^2;       % rotational energy, J
Etot    = sum(Ekin,2) + sum(Erot,2) + Epot;              % total energy, J

e  = vr(end)/vr(1);                                      % total restitution coefficient, -
en = abs((v(1,1,end) - v(1,2,end))/(v(1,1,1)-v(1,2,1))); % normal restitution coefficient, -
et = abs((v(2,1,end) - v(2,2,end))/(v(2,1,1)-v(2,2,1))); % tangential restitution coefficient, -

thetai0 = atan(v(2,1,1  )/v(1,1,1  )); % impact angle i, rad
thetaj0 = 0;                           % impact angle j, rad
thetai1 = atan(v(2,1,end)/v(1,1,end)); % rebound angle i, rad
thetaj1 = atan(v(2,2,end)/v(1,2,end)); % rebound angle j, rad

Sigma = abs(thetaj1-thetai1)/(pi/2);   % effective rebound angle, -, see: http://dx.doi.org/10.1063/1.2396925

%% display results
disp(['St = ' num2str(St,4) '; theta0 = ' num2str(rad2deg(theta0),2) ' deg; e = ' num2str(e,2) '; en = ' num2str(en,2) '; et = ' num2str(et,2) '; Sigma = ' num2str(Sigma,2)])