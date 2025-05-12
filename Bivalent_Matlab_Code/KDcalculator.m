function [KD] = KDcalculator(RUmax)
RUaim=RUmax/2;

Ltest=1e-6;
KDerror=RUaim*0.01;
OnTime=1e7;
SignalMax=KDcalculatorSIgnal(Ltest, OnTime);
flat=1;

while RUaim>SignalMax+KDerror || RUaim<SignalMax-KDerror || flat==0
    
    if SignalMax>RUmax*0.95
        Ltest=Ltest/10;
    elseif SignalMax<RUmax*0.05
        Ltest=Ltest*10;
    elseif RUaim<SignalMax-KDerror % [Ligand] too high
        Ltest=Ltest*(RUaim/SignalMax);
    elseif RUaim>SignalMax+KDerror % [Ligand] too low
        Ltest=Ltest*(RUaim/SignalMax);
    end
        
        
    [SignalMax, Y, T]=KDcalculatorSIgnal(Ltest, OnTime);
    Ymin=min(find(T>OnTime/2)); %#ok<MXFND>
    Ymax=max(find(T<OnTime)); %#ok<MXFND>
    if Y(Ymin)>(Y(Ymax)*1.0001) || Y(Ymin)<(Y(Ymax)*0.9999)
        OnTime=OnTime*100;
        flat=0;
    else
        flat=1;
    end
    
    if Ltest<1e-22
        KD='covalent like (0)';
        return
    end
end


KD=Ltest;
disp(SignalMax);
end

function [SignalMax, Signal, T] =KDcalculatorSIgnal(Ltest, OnTime)

global L R0 KonTable MWligand Vol Area Ruconv
%% 1) Differential equation calculation
%Set up the options for the ODE solver
options = odeset('RelTol',1e-12,'AbsTol',1e-13);
%OnTime=1000;
OffTime=100;

%% 1.1) %Set up the association phase parameters and calculate the equations

L=Ltest; %#ok<NASGU>                            %The global L used in the differention equations is set to the association value
range = [0 OnTime];                             %The range of the association
y0 = [R0, zeros(1, length(KonTable)-1)];        %The initial condition where Y1=R0 and the others are 0
[T,Y] = ode15s(@Diffs,range,y0, options);       %Using the created DiffsAss,m function for the association phase
TA=T; YA=Y;                                     %The output of the association phase

%% 1.2) Set up the dissosiation phase parameters and calculate the equations

L=0;                                            %The global L used in the differention equations is set to the dissociation value
range=[0 OffTime];                              %The range of the dissociation
y0 = YA(size(YA,1),:);                          %The initial conditions are the end point of the association phase
[T,Y] = ode15s(@Diffs,range,y0,options);        %Useing the created DiffsDiss,m function for the dissociation phase
TD=T; YD=Y;                                     %The output of the dissociation phase

%Giving T and Y the values of the whole experiment
T=[TA; TD+TA(end)];  Y=[YA; YD];                %Add together the values of the Ass and Diss phases

%% 2) Interpret the result
RUmod=(MWligand*Vol)/(Area*Ruconv);
Signal=((R0-Y(:,1))*RUmod); 
SignalMax=max(Signal);
end