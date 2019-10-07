%asking user how they want to use the program
PromptOptSel = 'Select Option: (1)Van der Waals (2)Peng-Robinson (3)Cardano Cubic Solver - Enter Number: ';
OptSel = input(PromptOptSel);

if OptSel==1 %Defining values as per VDW equation
    disp ' Program is running Van der Waals EoS';
    %getting inputs for VDW
    PromptPc = 'What is critical Pressure in Pascals? Enter Here: ';
    Pc = input(PromptPc);
    PromptTc = 'What is critical Temperature in Kelvin? Enter Here: ';
    Tc = input(PromptTc);
    PromptPr = 'What is Test Pressure in Pascals? Enter Here: ';
    Pr = input(PromptPr);
    PromptTr = 'What is Test Temperature in Kelvin? Enter Here: ';
    Tr = input(PromptTr);
    %Promptm = 'What is molecular mass in Kg/mol? Enter Here: ';
    %m = input(Promptm);

    %estimating VDW a and b
    a=(27*(8.3144598^2)*Tc^2)/(64*Pc);
    b=(8.3144598*Tc)/(8*Pc);

    %defining cubic parameters constants
    A=(((-Pr*b)+(-8.3144598*Tr))/(Pr));
    B=a/(Pr);
    C=(-a*b)/(Pr);

elseif OptSel==2
    disp 'Program is running Peng-Robinson EoS';
    %getting inputs for Peng-Robinson
    PromptPc = 'What is critical Pressure in Pascals? Enter Here: ';
    Pc = input(PromptPc);
    PromptTc = 'What is critical Temperature in Kelvin? Enter Here: ';
    Tc = input(PromptTc);
    PromptPr = 'What is Test Pressure in Pascals? Enter Here: ';
    Pr = input(PromptPr);
    PromptTr = 'What is Test Temperature in Kelvin? Enter Here: ';
    Tr = input(PromptTr);
    %Promptm = 'What is molecular mass in Kg/mol? Enter Here: ';
    %m = input(Promptm);
    %Mol = input('Enter number of moles: ',Mol);
    PromptOm = 'What is acentric factor? Enter Here as X where X is X.e-3: ';
    Om = input(PromptOm);
    
    %estimating PR parameters
    Tred=Tr/Tc;
    %Om=str2double(Omega);
    Omega=Om*1e-3;
    if Omega>0.49
        Kappa=((0.37964)+(Omega*((1.48503)+(Omega*(-(0.164423)+((0.016666)*Omega))))));
    else
        Kappa=(0.37464+(1.54226*Omega)-(0.26992*((Omega*Omega))));
    end
    Alpha=(1+(Kappa*(1-(Tred^0.5))))^2;
    a=(0.45724*(8.3144598^2)*(Tc^2)*Alpha)/(Pc);
    b=(0.07780*8.3144598*Tc)/(Pc);
    PRA=(a*Pr)/((8.3144598^2)*Tr^2);
    PRB=(b*Pr)/(8.3144598*Tr);


    %defining cubic parameters constants
    A=-1+PRB;
    B=PRA-(3*(PRB^2))-(2*PRB);
    C=-(PRA*PRB)+(PRB^2)+(PRB^3);
    
elseif OptSel==3
    disp 'Running Cubic Solver. Enter the coefficients for cubic equatiuon in the following order: x^3 + Ax^2 + Bx + C = 0 ';
    disp 'Note Matlab doesnt accept inputs with decimals. For e.g: Enter 0.5 as 5e-2 ';
    PromptA = 'Enter A Here: ';
    A = input(PromptA);
    PromptB = 'Enter B Here: ';
    B = input(PromptB);
    PromptC = 'Enter C Here: ';
    C = input(PromptC);
else
    %prompting user to enter 1 2 or 3 instead of any other number
    disp ' Please Enter 1, 2 or 3. Restart program to select again. ';
end
    
    
%Running Cardano solution
D=((A/3)^3)-((A*B)/6)+(C/2);
E=(B/3)-((A/3)^2);

Del=(D^2)+(E^3);

%Solving the cubic equation based on value of Z

if Del==0
    x1=(2*(-D^(1/3)))-(A/3);
    x2=(-1*(-D^(1/3)))-(A/3);
    x3=x2;
elseif Del>0
    F=nthroot(((-D)+((Del^0.5))),3);
    G=nthroot(((-D)-((Del^0.5))),3);
    x1=F+G-(A/3);
    x2=-((0.5*(F+G))+(A/3))+((((3^0.5)/2)*(F-G)))*1i;
    x3=-((0.5*(F+G))+(A/3))-((((3^0.5)/2)*(F-G)))*1i;
elseif Del<0
    Ang=acos(-D/((-E^3)^0.5));
    x1=(2*((-E)^0.5)*cos(Ang/3))-(A/3);
    x2=(2*((-E)^0.5)*cos((Ang/3)+((2*pi)/3)))-(A/3);
    x3=(2*((-E)^0.5)*cos((Ang/3)+((4*pi)/3)))-(A/3);
else
end

%displaying results
if OptSel==2
    disp 'These roots are Z - Compressibility Factors not Volume: '
    disp 'z1: ';
    disp(x1)
    disp 'z2: ';
    disp(x2)
    disp 'z3: ';
    disp(x3)
elseif OptSel==1
    disp 'These are molar volumes (m3/mol): '
    disp 'V1: ';
    disp(x1)
    disp 'V2: ';
    disp(x2)
    disp 'V3: ';
    disp(x3)
else
    disp 'Cubic roots are: '
    disp 'x1: ';
    disp(x1)
    disp 'x2: ';
    disp(x2)
    disp 'x3: ';
    disp(x3)
end

%Asking to continue Program
disp 'Do you want to estiamate Vapor Pressure from PR or end? '
PromptEngProg = 'Enter (1) to End Program, (2) to continue. Enter Here: ';
EngProg = input(PromptEngProg);

if EndProg == 1
    return
else
end

%Estimating Vapor Pressure
if Del>0
    disp 'Entered Pressure is in the vapor only regime. Enter different pressure and restart.'
    return
elseif Del<0
    disp 'Three real roots are found. Ignoring unstable root.'
    roots=[x1,x2,x3];
    Posroots=roots(roots>=0);
    Vl=min(x1,x2,x3);
    Vp=max(x1,x2,x3);
    Vp=(x3*8.3144598*Tr)/Pr;
    disp 'Liquid volume in m3/mol is: '
    disp(Vl)
    disp 'Vapor volume in m3/mol is: '
    disp(Vp)
else
end


    
    



%Author: Imran Tayyab
