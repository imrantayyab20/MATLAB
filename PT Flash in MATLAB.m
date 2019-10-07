Nc=input('How many components?');
Pres=input('Pressure for calcultaion [bar]?');
disp(' Use Temperature in Kelvin');
Temp=input('Temperature for calculation [K]?');
A = importdata('dataf.txt', ' ', 1);
BC=importdata('BIP.txt',' ');

Z=A.data(:, 1);
Tc=A.data(:, 2);
Pc=A.data(:, 3);
omega=A.data(:, 4);
BP=BC.data;

RVal = 0.00008314 ;

disp('Performing Stability Analysis');

 % Calculating Reduced Properties
for i = 1: 1: Nc
    Tr(i) = Temp/Tc(i);
    Pr (i) = Pres/Pc(i);
end
% Using If Condition (Kappa Values), Alpha values, ai, bi, BM
for i=1:1:Nc
if omega (i) < 0.49
     kappa (i) = 0.37464+1.54226*omega (i)-0.26992*omega (i)^2;
    else 
     kappa (i) = 0.379642+1.48503*omega(i)-0.1644*omega(i)^2+0.016667*omega(i)^3;
    end
end
for i=1:1:Nc
    Alpha (i) = (1+(kappa (i)*(1-sqrt(Tr (i)))))^2;
    ai (i) = 0.457236*RVal^2* Tc (i)^2* Alpha (i)/Pc (i);
end
for i= 1:1:Nc
    bi (i) = 0.0778*RVal*Tc (i)/Pc(i);
    Xb(i)= Z(i) * bi(i);
end
for i= 1:1:Nc
    Bi (i) = bi(i)*Pres/RVal/Temp ;
end
 
bmix = sum (Xb);
for i = 1:1: Nc
    for j = 1:1:Nc
    ami(i,j) = (1- BP (i,j) )*sqrt (ai (i)* ai(j));
    cm (i,j)= Z (i)* Z(j);
    end
end
for i=1:1:Nc
    for j=1:1:Nc
    amcm (i, j) = ami(i,j) * cm (i,j);
    end
end   
SUMamcm = sum (amcm, 1) ;
amix= sum (SUMamcm) ;
Ami = amix*Pres/RVal^2/Temp^2;
Bmi = bmix*Pres/RVal/Temp ;

%cubic PR EOS Z^3+(-1+B)Z^2+(A-3B^2-2B)Z+(-AB+B^2+B^3)=0
    
%coeffcient of z^2
A1=-1+Bmi;
%coefficient of z
B1=Ami-3*(Bmi)^2-2*Bmi;
%constant
C1=-Ami*Bmi+(Bmi)^2+(Bmi)^3;

D=(A1/3)^3-((A1*B1)/6)+(C1/2);

E=(B1/3)-(A1/3)^2;

Delta=(D)^2+(E)^3;

if Delta==0
x1=2*nthroot((-1*D),3)-(A1/3);
x2=-1*nthroot((-1*D),3)-(A1/3);
x3=-1*nthroot((-1*D),3)-(A1/3);

elseif Delta>0
    sqrtdel=sqrt(Delta);
    
    F=nthroot(((-1*D)+sqrtdel),3);
    
    G=nthroot(((-1*D)-sqrtdel),3);
    
    x1=F+G-(A1/3);
    x2=-1*(0.5*(F+G)+(A1/3))+((sqrt(3)/2)*(F-G))*1i;
    x3=-1*(0.5*(F+G)+(A1/3))-((sqrt(3)/2)*(F-G))*1i;
    else
        theta=acos((-1*D)/sqrt((-E)^3));
        x1=2*sqrt(-1*E)*cos(theta/3)-(A1/3);
        x2=2*sqrt(-1*E)*cos((theta/3)+(2/3)*pi)-(A1/3);
        x3=2*sqrt(-1*E)*cos((theta/3)+(4/3)*pi)-(A1/3);
        
end 

roots=[x1,x2,x3];
roots=roots(imag(roots)==0);

if Delta>0
    r=x1;
else
    r1=max(roots);
    r3=min(roots);

    for i = 1:1:Nc
    VAR(i) = 0;
    for j = 1:Nc
        VAR(i) = VAR(i) + Z(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end
end
% Calculating Fugacities 
for i=1:1:Nc
    lnphi1(i) = (r1-1)*Bi(i)/Bmi - log(r1-Bmi) + ((-Ami)/Bmi/(2*sqrt(2))) * (2*VAR(i)/Ami - Bi(i)/Bmi) * log ((r1+(1+sqrt(2))*Bmi)/(r1+(1-sqrt(2))*Bmi));
    lnphi3(i) = (r3-1)*Bi(i)/Bmi - log(r3-Bmi) + ((-Ami)/Bmi/(2*sqrt(2))) * (2*VAR(i)/Ami - Bi(i)/Bmi) * log ((r3+(1+sqrt(2))*Bmi)/(r3+(1-sqrt(2))*Bmi));
end

DelGdimen=0;
for i=1:1:Nc
    DelGdimen = DelGdimen + (Z(i) * (lnphi1(i) - lnphi3(i)));
    
    if DelGdimen>0
        r=r3;
    else
        r=r1;
    end
end
end
 %now we have our roots
 
 for i=1:1:Nc
   VAR(i) = 0;
    for j = 1:Nc
        VAR(i) = VAR(i) + Z(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end  
 end
 
 for i=1:1:Nc
    lnphi(i) = (r-1)*Bi(i)/Bmi - log(r-Bmi) + ((-Ami)/Bmi/(2*sqrt(2))) * (2*VAR(i)/Ami - Bi(i)/Bmi) * log ((r+(1+sqrt(2))*Bmi)/(r+(1-sqrt(2))*Bmi));
    fugacity(i) = Pres*Z(i)*exp(lnphi(i));
    phi(i)=exp(lnphi(i));
 end

 disp(' The Fugacities for each component in the mixture are =');
 disp(fugacity);
 disp( 'The fugacity coefficients of each component in mixture are =');
 disp(lnphi);
 
% For Vapor Phase

for i=1:Nc
    K(i)= exp(5.37*(1+omega(i))*(1-Tr(i)^-1))/Pr(i);
end

for i=1:Nc
    Y(i)=K(i)*Z(i);
end

criteriav=1;
iterationlimit=25;
u=0;
while criteriav>0.000001

u=u+1;
if u>=iterationlimit
    break
end
sumY=sum(Y);

for i=1:Nc
    y(i)=Y(i)/sumY;
end


% ****Calculating Vapor Phase Fugacity****

for i=1:Nc
    Xbv(i)=y(i)*bi(i);
end

bmixv=sum(Xbv);

for i=1:Nc
    for j=1:Nc
        cmv(i,j)=y(i)*y(j);
    end
end

for i=1:Nc
    for j=1:Nc
        amcmv(i,j)=ami(i,j)*cmv(i,j);
    end
end

sumamcmv=sum(amcmv,1);
amixv=sum(sumamcmv);

Amiv=amixv*Pres/RVal^2/Temp^2;
Bmiv=bmixv*Pres/RVal/Temp;

% Calculating Roots

%cubic PR EOS Z^3+(-1+B)Z^2+(A-3B^2-2B)Z+(-AB+B^2+B^3)=0
    
%coeffcient of z^2
A1v=-1+Bmiv;
%coefficient of z
B1v=Amiv-3*(Bmiv)^2-2*Bmiv;
%constant
C1v=-Amiv*Bmiv+(Bmiv)^2+(Bmiv)^3;

Dv=(A1v/3)^3-((A1v*B1v)/6)+(C1v/2);

Ev=(B1v/3)-(A1v/3)^2;

Deltav=(Dv)^2+(Ev)^3;

if Deltav==0
x1v=2*nthroot((-1*Dv),3)-(A1v/3);
x2v=-1*nthroot((-1*Dv),3)-(A1v/3);
x3v=-1*nthroot((-1*Dv),3)-(A1v/3);

elseif Deltav>0
    sqrtdelv=sqrt(Deltav);
    
    Fv=nthroot(((-1*Dv)+sqrtdelv),3);
    
    Gv=nthroot(((-1*Dv)-sqrtdelv),3);
    
    x1v=Fv+Gv-(A1v/3);
    x2v=-1*(0.5*(Fv+Gv)+(A1v/3))+((sqrt(3)/2)*(Fv-Gv))*1i;
    x3v=-1*(0.5*(Fv+Gv)+(A1v/3))-((sqrt(3)/2)*(Fv-Gv))*1i;
    else
        thetav=acos((-1*Dv)/sqrt((-Ev)^3));
        x1v=2*sqrt(-1*Ev)*cos(thetav/3)-(A1v/3);
        x2v=2*sqrt(-1*Ev)*cos((thetav/3)+(2/3)*pi)-(A1v/3);
        x3v=2*sqrt(-1*Ev)*cos((thetav/3)+(4/3)*pi)-(A1v/3);
        
end 

rootsv=[x1v,x2v,x3v];
rootsv=rootsv(imag(rootsv)==0);

if Deltav>0
    rv=x1v;
else
    r1v=max(rootsv);
    r3v=min(rootsv);

    for i = 1:1:Nc
    VARv(i) = 0;
    for j = 1:Nc
        VARv(i) = VARv(i) + y(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end
end
% Calculating Fugacities 
for i=1:1:Nc
    lnphi1v(i) = (r1v-1)*Bi(i)/Bmiv - log(r1v-Bmiv) + ((-Amiv)/Bmiv/(2*sqrt(2))) * (2*VARv(i)/Amiv - Bi(i)/Bmiv) * log ((r1v+(1+sqrt(2))*Bmiv)/(r1v+(1-sqrt(2))*Bmiv));
    lnphi3v(i) = (r3v-1)*Bi(i)/Bmiv - log(r3v-Bmiv) + ((-Amiv)/Bmiv/(2*sqrt(2))) * (2*VARv(i)/Amiv - Bi(i)/Bmiv) * log ((r3v+(1+sqrt(2))*Bmiv)/(r3v+(1-sqrt(2))*Bmiv));
end

DelGdimenv=0;
for i=1:1:Nc
    DelGdimenv = DelGdimenv + (y(i) * (lnphi1v(i) - lnphi3v(i)));
    
    if DelGdimenv>0
        rv=r3v;
    else
        rv=r1v;
    end
end
end
 %now we have our roots
 
 for i=1:1:Nc
   VARv(i) = 0;
    for j = 1:Nc
        VARv(i) = VARv(i) + y(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end  
 end
 
 for i=1:1:Nc
    lnphiv(i) = (rv-1)*Bi(i)/Bmiv - log(rv-Bmiv) + ((-Amiv)/Bmiv/(2*sqrt(2))) * (2*VARv(i)/Amiv - Bi(i)/Bmiv) * log ((rv+(1+sqrt(2))*Bmiv)/(rv+(1-sqrt(2))*Bmiv));
    fugacityv(i) = Pres*y(i)*exp(lnphiv(i));
    phiv(i)=exp(lnphiv(i));
 end
 
 for i=1:Nc
     convvapor(i)=abs((log(Y(i))+lnphiv(i))-(log(Z(i))+lnphi(i)));
 end
 
 criteriav=max(convvapor);
 
 for i=1:Nc
     Y(i)=Z(i)*phi(i)/phiv(i);
 end
 
end

disp('the number of iterations to acheive convergence for vapor phase are =');
disp(u);
disp('the value of convergence criteria for vapor phase at end of loop was =');
disp(criteriav);
disp('The sum of vapor mole numbers after acheiving convergence is=');
disp(sumY);

if sumY>1
    
    disp('Phase is unstable');
    
else

% Liquid phase

for i=1:Nc
    Kl(i)= exp(5.37*(1+omega(i))*(1-Tr(i)^-1))/Pr(i);
end

for i=1:Nc
    X(i)=Z(i)/Kl(i);
end

criterial=1;
iterationlimit=50;
t=0;
while criterial>0.000001

t=t+1;
if t>=iterationlimit
    break
end
sumX=sum(X);

for i=1:Nc
    x(i)=X(i)/sumX;
end


% ****Calculating liquid Phase Fugacity****

for i=1:Nc
    Xbl(i)=x(i)*bi(i);
end

bmixl=sum(Xbl);

for i=1:Nc
    for j=1:Nc
        cml(i,j)=x(i)*x(j);
    end
end

for i=1:Nc
    for j=1:Nc
        amcml(i,j)=ami(i,j)*cml(i,j);
    end
end

sumamcml=sum(amcml,1);
amixl=sum(sumamcml);

Amil=amixl*Pres/RVal^2/Temp^2;
Bmil=bmixl*Pres/RVal/Temp;

% Calculating Roots

%cubic PR EOS Z^3+(-1+B)Z^2+(A-3B^2-2B)Z+(-AB+B^2+B^3)=0
    
%coeffcient of z^2
A1l=-1+Bmil;
%coefficient of z
B1l=Amil-3*(Bmil)^2-2*Bmil;
%constant
C1l=-Amil*Bmil+(Bmil)^2+(Bmil)^3;

Dl=(A1l/3)^3-((A1l*B1l)/6)+(C1l/2);

El=(B1l/3)-(A1l/3)^2;

Deltal=(Dl)^2+(El)^3;

if Deltal==0
x1l=2*nthroot((-1*Dl),3)-(A1l/3);
x2l=-1*nthroot((-1*Dl),3)-(A1l/3);
x3l=-1*nthroot((-1*Dl),3)-(A1l/3);

elseif Deltal>0
    sqrtdell=sqrt(Deltal);
    
    Fl=nthroot(((-1*Dl)+sqrtdell),3);
    
    Gl=nthroot(((-1*Dl)-sqrtdell),3);
    
    x1l=Fl+Gl-(A1l/3);
    x2l=-1*(0.5*(Fl+Gl)+(A1l/3))+((sqrt(3)/2)*(Fl-Gl))*1i;
    x3l=-1*(0.5*(Fl+Gl)+(A1l/3))-((sqrt(3)/2)*(Fl-Gl))*1i;
    else
        thetal=acos((-1*Dl)/sqrt((-El)^3));
        x1l=2*sqrt(-1*El)*cos(thetal/3)-(A1l/3);
        x2l=2*sqrt(-1*El)*cos((thetal/3)+(2/3)*pi)-(A1l/3);
        x3l=2*sqrt(-1*El)*cos((thetal/3)+(4/3)*pi)-(A1l/3);
        
end 

rootsl=[x1l,x2l,x3l];
rootsl=rootsl(imag(rootsl)==0);

if Deltal>0
    rl=x1l;
else
    r1l=max(rootsl);
    r3l=min(rootsl);

    for i = 1:1:Nc
    VARl(i) = 0;
    for j = 1:Nc
        VARl(i) = VARl(i) + x(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end
end
% Calculating Fugacities 
for i=1:1:Nc
    lnphi1l(i) = (r1l-1)*Bi(i)/Bmil - log(r1l-Bmil) + ((-Amil)/Bmil/(2*sqrt(2))) * (2*VARl(i)/Amil - Bi(i)/Bmil) * log ((r1l+(1+sqrt(2))*Bmil)/(r1l+(1-sqrt(2))*Bmil));
    lnphi3l(i) = (r3l-1)*Bi(i)/Bmil - log(r3l-Bmil) + ((-Amil)/Bmil/(2*sqrt(2))) * (2*VARl(i)/Amil - Bi(i)/Bmil) * log ((r3l+(1+sqrt(2))*Bmil)/(r3l+(1-sqrt(2))*Bmil));
end

DelGdimenl=0
for i=1:1:Nc
    DelGdimenl = DelGdimenl + (x(i) * (lnphi1l(i) - lnphi3l(i)));
    
    if DelGdimenl>0
        rl=r3l;
    else
        rl=r1l;
    end
end
end
 %now we have our roots
 
 for i=1:1:Nc
   VARl(i) = 0;
    for j = 1:Nc
        VARl(i) = VARl(i) + x(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end  
 end
 
 for i=1:1:Nc
    lnphil(i) = (rl-1)*Bi(i)/Bmil - log(rl-Bmil) + ((-Amil)/Bmil/(2*sqrt(2))) * (2*VARl(i)/Amil - Bi(i)/Bmil) * log ((rl+(1+sqrt(2))*Bmil)/(rl+(1-sqrt(2))*Bmil));
    fugacityl(i) = Pres*x(i)*exp(lnphil(i));
    phil(i)=exp(lnphil(i));
 end
 
 for i=1:Nc
     convliquid(i)=abs((log(X(i))+lnphil(i))-(log(Z(i))+lnphi(i)));
 end
 
 criterial=max(convliquid);
 
 for i=1:Nc
     X(i)=Z(i)*phi(i)/phil(i);
 end
 
end

disp('the number of iterations to acheive convergence for liquid phase are =');
disp(t);
disp('the value of convergence criteria for liquid phase at end of loop was =');
disp(criterial);
disp('The sum of liquid mole numbers after acheiving convergence is=');
disp(sumX);

end

if sumY>1 || sumX>1
    disp('The single phase is unstable');
elseif sumY==1 || sumX==1
    disp('Trivial solution is reached for one of the phases');
else
    disp('the phase is stable');
end

% *************************Two Phase P-T flash****************************

if sumY>1 || sumX>1
    
disp('Proceeding with Two-phase PT flash');

criteria2p=1;
iterationlimit2p=100;
m=0;
while criteria2p>0.000001

m=m+1;
if m>=iterationlimit2p
    break
end
    
Kmin=min(K)
Kmax=max(K)

betamin=1/(1-Kmax);
betamax=1/(1-Kmin);

beta = (betamin +betamax) / 2 

convergance=1;

while convergance>1e-6

% .....Using newton Raphson Method.....%
for i=1:1:Nc
    Fbeta(i)=(Z(i)*(K(i)-1)/(1+(beta*(K(i)-1))));
end
Fbeta1=sum(Fbeta);

for i=1:1:Nc
    dFdbeta(i)=(Z(i)*((K(i)-1)^2)/(((beta*(K(i)-1))+1)^2));
end
dFdbeta1=-(sum(dFdbeta));
betanew=beta-(Fbeta1/dFdbeta1)
convergance = abs (Fbeta1)
% ......Bisection Method.......... 
if betanew<betamin
    betanew=0.5*(beta+betamin);
elseif betanew>betamax
    betanew=0.5*(beta+betamax);
end

beta = betanew;
end

% Calculating Xi and Yi
for i=1:1:Nc
    x2p(i) =Z(i)/(betanew*(K(i)-1)+1)
    y2p(i)= K(i)* x2p(i)
end

% for Vapor Phase 

for i=1:Nc
    Xb2pv(i)=y2p(i)*bi(i);
end

bmix2pv=sum(Xb2pv);

for i=1:Nc
    for j=1:Nc
        cm2pv(i,j)=y2p(i)*y2p(j);
    end
end

for i=1:Nc
    for j=1:Nc
        amcm2pv(i,j)=ami(i,j)*cm2pv(i,j);
    end
end

sumamcm2pv=sum(amcm2pv,1);
amix2pv=sum(sumamcm2pv);

Ami2pv=amix2pv*Pres/RVal^2/Temp^2;
Bmi2pv=bmix2pv*Pres/RVal/Temp;

% Calculating Roots

%cubic PR EOS Z^3+(-1+B)Z^2+(A-3B^2-2B)Z+(-AB+B^2+B^3)=0
    
%coeffcient of z^2
A12pv=-1+Bmi2pv;
%coefficient of z
B12pv=Ami2pv-3*(Bmi2pv)^2-2*Bmi2pv;
%constant
C12pv=-Ami2pv*Bmi2pv+(Bmi2pv)^2+(Bmi2pv)^3;

D2pv=(A12pv/3)^3-((A12pv*B12pv)/6)+(C12pv/2);

E2pv=(B12pv/3)-(A12pv/3)^2;

Delta2pv=(D2pv)^2+(E2pv)^3;

if Delta2pv==0
x12pv=2*nthroot((-1*D2pv),3)-(A12pv/3);
x22pv=-1*nthroot((-1*D2pv),3)-(A12pv/3);
x32pv=-1*nthroot((-1*D2pv),3)-(A12pv/3);

elseif Delta2pv>0
    sqrtdel2pv=sqrt(Delta2pv);
    
    F2pv=nthroot(((-1*D2pv)+sqrtdel2pv),3);
    
    G2pv=nthroot(((-1*D2pv)-sqrtdel2pv),3);
    
    x12pv=F2pv+G2pv-(A12pv/3);
    x22pv=-1*(0.5*(F2pv+G2pv)+(A12pv/3))+((sqrt(3)/2)*(F2pv-G2pv))*1i;
    x32pv=-1*(0.5*(F2pv+G2pv)+(A12pv/3))-((sqrt(3)/2)*(F2pv-G2pv))*1i;
    else
        theta2pv=acos((-1*D2pv)/sqrt((-E2pv)^3));
        x12pv=2*sqrt(-1*E2pv)*cos(theta2pv/3)-(A12pv/3);
        x22pv=2*sqrt(-1*E2pv)*cos((theta2pv/3)+(2/3)*pi)-(A12pv/3);
        x32pv=2*sqrt(-1*E2pv)*cos((theta2pv/3)+(4/3)*pi)-(A12pv/3);
        
end 

roots2pv=[x12pv,x22pv,x32pv];
roots2pv=roots2pv(imag(roots2pv)==0);

if Delta2pv>0
    r2pv=x12pv;
else
    r12pv=max(roots2pv);
    r32pv=min(roots2pv);

    for i = 1:1:Nc
    VAR2pv(i) = 0;
    for j = 1:Nc
        VAR2pv(i) = VAR2pv(i) + y2p(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end
end
% Calculating Fugacities 
for i=1:1:Nc
    lnphi12pv(i) = (r12pv-1)*Bi(i)/Bmi2pv - log(r12pv-Bmi2pv) + ((-Ami2pv)/Bmi2pv/(2*sqrt(2))) * (2*VAR2pv(i)/Ami2pv - Bi(i)/Bmi2pv) * log ((r12pv+(1+sqrt(2))*Bmi2pv)/(r12pv+(1-sqrt(2))*Bmi2pv));
    lnphi32pv(i) = (r32pv-1)*Bi(i)/Bmi2pv - log(r32pv-Bmi2pv) + ((-Ami2pv)/Bmi2pv/(2*sqrt(2))) * (2*VAR2pv(i)/Ami2pv - Bi(i)/Bmi2pv) * log ((r32pv+(1+sqrt(2))*Bmi2pv)/(r32pv+(1-sqrt(2))*Bmi2pv));
end

DelGdimen2pv=0;
for i=1:1:Nc
    DelGdimen2pv = DelGdimen2pv + (y2p(i) * (lnphi12pv(i) - lnphi32pv(i)));
    
    if DelGdimen2pv>0
        r2pv=r32pv;
    else
        r2pv=r12pv;
    end
end
end
 %now we have our roots
 
 for i=1:1:Nc
   VAR2pv(i) = 0;
    for j = 1:Nc
        VAR2pv(i) = VAR2pv(i) + y2p(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end  
 end
 
 for i=1:1:Nc
    lnphi2pv(i) = (r2pv-1)*Bi(i)/Bmi2pv - log(r2pv-Bmi2pv) + ((-Ami2pv)/Bmi2pv/(2*sqrt(2))) * (2*VAR2pv(i)/Ami2pv - Bi(i)/Bmi2pv) * log ((r2pv+(1+sqrt(2))*Bmi2pv)/(r2pv+(1-sqrt(2))*Bmi2pv))
    fugacity2pv(i) = Pres*y2p(i)*exp(lnphi2pv(i))
    phi2pv(i)=exp(lnphi2pv(i))
 end
 
 %************************for Liquid phase *******************************
  
 for i=1:Nc
    Xb2pl(i)=x2p(i)*bi(i);
end

bmix2pl=sum(Xb2pl);

for i=1:Nc
    for j=1:Nc
        cm2pl(i,j)=x2p(i)*x2p(j);
    end
end

for i=1:Nc
    for j=1:Nc
        amcm2pl(i,j)=ami(i,j)*cm2pl(i,j);
    end
end

sumamcm2pl=sum(amcm2pl,1);
amix2pl=sum(sumamcm2pl);

Ami2pl=amix2pl*Pres/RVal^2/Temp^2;
Bmi2pl=bmix2pl*Pres/RVal/Temp;

% Calculating Roots

%cubic PR EOS Z^3+(-1+B)Z^2+(A-3B^2-2B)Z+(-AB+B^2+B^3)=0
    
%coeffcient of z^2
A12pl=-1+Bmi2pl;
%coefficient of z
B12pl=Ami2pl-3*(Bmi2pl)^2-2*Bmi2pl;
%constant
C12pl=-Ami2pl*Bmi2pl+(Bmi2pl)^2+(Bmi2pl)^3;

D2pl=(A12pl/3)^3-((A12pl*B12pl)/6)+(C12pl/2);

E2pl=(B12pl/3)-(A12pl/3)^2;

Delta2pl=(D2pl)^2+(E2pl)^3;

if Delta2pl==0
x12pl=2*nthroot((-1*D2pl),3)-(A12pl/3);
x22pl=-1*nthroot((-1*D2pl),3)-(A12pl/3);
x32pl=-1*nthroot((-1*D2pl),3)-(A12pl/3);

elseif Delta2pl>0
    sqrtdel2pl=sqrt(Delta2pl);
    
    F2pl=nthroot(((-1*D2pl)+sqrtdel2pl),3);
    
    G2pl=nthroot(((-1*D2pl)-sqrtdel2pl),3);
    
    x12pl=F2pl+G2pl-(A12pl/3);
    x22pl=-1*(0.5*(F2pl+G2pl)+(A12pl/3))+((sqrt(3)/2)*(F2pl-G2pl))*1i;
    x32pl=-1*(0.5*(F2pl+G2pl)+(A12pl/3))-((sqrt(3)/2)*(F2pl-G2pl))*1i;
    else
        theta2pl=acos((-1*D2pl)/sqrt((-E2pl)^3));
        x12pl=2*sqrt(-1*E2pl)*cos(theta2pl/3)-(A12pl/3);
        x22pl=2*sqrt(-1*E2pl)*cos((theta2pl/3)+(2/3)*pi)-(A12pl/3);
        x32pl=2*sqrt(-1*E2pl)*cos((theta2pl/3)+(4/3)*pi)-(A12pl/3);
        
end 

roots2pl=[x12pl,x22pl,x32pl];
roots2pl=roots2pl(imag(roots2pl)==0);

if Delta2pl>0
    r2pl=x12pl;
else
    r12pl=max(roots2pl);
    r32pl=min(roots2pl);

    for i = 1:1:Nc
    VAR2pl(i) = 0;
    for j = 1:Nc
        VAR2pl(i) = VAR2pl(i) + x2p(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end
end
% Calculating Fugacities 
for i=1:1:Nc
    lnphi12pl(i) = (r12pl-1)*Bi(i)/Bmi2pl - log(r12pl-Bmi2pl) + ((-Ami2pl)/Bmi2pl/(2*sqrt(2))) * (2*VAR2pl(i)/Ami2pl - Bi(i)/Bmi2pl) * log ((r12pl+(1+sqrt(2))*Bmi2pl)/(r12pl+(1-sqrt(2))*Bmi2pl));
    lnphi32pl(i) = (r32pl-1)*Bi(i)/Bmi2pl - log(r32pl-Bmi2pl) + ((-Ami2pl)/Bmi2pl/(2*sqrt(2))) * (2*VAR2pl(i)/Ami2pl - Bi(i)/Bmi2pl) * log ((r32pl+(1+sqrt(2))*Bmi2pl)/(r32pl+(1-sqrt(2))*Bmi2pl));
end

DelGdimen2pl=0
for i=1:1:Nc
    DelGdimen2pl = DelGdimen2pl + (x2p(i) * (lnphi12pl(i) - lnphi32pl(i)));
    
    if DelGdimen2pl>0
        r2pl=r32pl;
    else
        r2pl=r12pl;
    end
end
end
 %now we have our roots
 
 for i=1:1:Nc
   VAR2pl(i) = 0;
    for j = 1:Nc
        VAR2pl(i) = VAR2pl(i) + x2p(j)*(ami(i,j))*Pres/RVal^2/Temp^2;
    end  
 end
 
 for i=1:1:Nc
    lnphi2pl(i) = (r2pl-1)*Bi(i)/Bmi2pl - log(r2pl-Bmi2pl) + ((-Ami2pl)/Bmi2pl/(2*sqrt(2))) * (2*VAR2pl(i)/Ami2pl - Bi(i)/Bmi2pl) * log ((r2pl+(1+sqrt(2))*Bmi2pl)/(r2pl+(1-sqrt(2))*Bmi2pl))
    fugacity2pl(i) = Pres*x2p(i)*exp(lnphi2pl(i))
    phi2pl(i)=exp(lnphi2pl(i))
 end
 
 for i=1:Nc
     conv2phase(i)=abs((log(x2p(i))+lnphi2pl(i))-(log(y2p(i))+lnphi2pv(i)));
 end
 
 criteria2p=max(conv2phase)
 
 for i=1:Nc
     K(i)=phi2pl(i)/phi2pv(i);
 end
end
disp('number of iterations to reach solution=');
disp(m);
disp('The vapor phase mole number is=');
disp(beta);
betal=1-beta;
disp('The liquid phase mole number is=');
disp(betal);
disp('The vapor phase composition is=');
disp(y2p);
disp('The liquid phase composition is=');
disp(x2p);
end