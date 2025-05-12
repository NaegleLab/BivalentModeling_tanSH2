function  Plotting(Lplot)
global plottype count SignalMat count2 Kon Kon3

switch plottype
    case 1
        count2 = 1;
        %series=[1,2,8,32,128,512,2028];
        series=[1,2,4,8,16,32,64,128,256,512,1024];
        %series=[1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65356];
        SignalMat = zeros(length(series),1);
        plotting1(Lplot)
        hold on
        count = 1;
        for i=2:(length(series))
            plotting1(Lplot/series(i))
            count = count + 1;
        end

    case 2
        plotting1(Lplot)
    case 3
        plotting1(Lplot,2)
end

end

function []= plotting1(Lplot, plottingtype)

if nargin==1
plottingtype=1;
end

global L OnTime OffTime R0 KonTable MWligand Vol Area Ruconv Kons11 ratioTable count2 SignalMat Kon Kon4 
%OnTime=500000
ratioTable = zeros(1000,1);
%currentRow = 1;
%% 1) Differential equation calculation

%Set up the options for the ODE solver
options = odeset('RelTol',1e-12,'AbsTol',1e-13);

%% 1.1) %Set up the association phase parameters and calculate the equations
global N M 

L=Lplot; %#ok<NASGU>                            %The global L used in the differention equations is set to association value
range = [0 OnTime];                             %The range of the association
if N==1 && M==1
    y0 = [R0, zeros(1, 4)];        %The Initial condition where Y1=R0 and the others are 0
elseif N==2 && M==2
    y0 = [R0, zeros(1, 16)];        %The Initial condition where Y1=R0 and the others are 0 was 13
end

[T,Y] = ode15s(@Diffs,range,y0, options);       %Using the created DiffsAss.m function for the Assiciation phase
TA=T; YA=Y;                                     %The output of the association phase

%% 1.2) Set up the dissociation phase parameters and calculate the equations

L=0;                                            %The global L used in the differention equations is set to dissociation value
range=[0 OffTime];                              %The range of the dissociation
y0 = YA(size(YA,1),:);                          %The initial conditions are the end point of the association phases
[T,Y] = ode15s(@Diffs,range,y0,options);        %Using the created DiffsDiss,m function for the dissociation phase
TD=T; YD=Y;                                     %The output of the dissociation phase

%Giving T and Y the values of the whole experiment
T=[TA; TD+TA(end)];  Y=[YA; YD];                %Add together the values of the Ass and Diss phases

%%
if N==1 && M==1
    multipyer=[0,1,1];
elseif N==2 && M==2
    multipyer=[0,1,1,1,1,2,2,2,2,2,2,2,2,1,1];        %was 0,1,1,1,2,2,1,2,2,2,1,2
end
    
%% 2) Interpret the results
RUmod=(MWligand*Vol)/(Area*Ruconv);
%Signal=((R0-Y(:,1))*RUmod);

Signal=((Y(:,1:end-2)*multipyer.')*10^9); %All signals
%Signal1 = ((Y(:,[1:5,14:end-2])*multipyer([1:5,14:end]).')*RUmod); %All signals except b
%Signal2 = ((Y(:,[14:end-2])*multipyer([14:end]).')*10^9); %only inline and twisted
%Signal3 = ((Y(:,[2])*multipyer([2]).')*10^9); %only a1x0
%Signal4 = ((Y(:,[5])*multipyer([5]).')*10^9); %only x0a2
%Signal = Signal1;
global Multiple 

if plottingtype==1 && Multiple~=1
    %Nonly = Signal3(round(0.5*length(Signal3)));
    %inline = (Y(:,14)*multipyer(14.')*RUmod);
    %twisted = (Y(:,15)*multipyer(15.')*RUmod);
    %Double = inline(round(0.5*length(inline)))+twisted(round(0.5*length(twisted)));
%     Double = Signal2;
%     p(1) = plot(T,Double,'LineWidth',2); %%
%     hold on
    if Kon == 1 | Kon4 ==1
        NSignal = Y(:,[2:3,6:9])*multipyer([2:3,6:9])'*10^9;
        %NSignal = Y(:,[4:5,10:13])*multipyer([4:5,10:13])'*10^9;
        plot(T,Signal,T,NSignal,'LineWidth',2)
        hold on
        plot(T(round(0.5*length(NSignal))),NSignal(round(0.5*length(NSignal))), 'ko')
        ylabel('Signal (nM)'); xlabel('Time (s)') %%
        Signal = NSignal(round(0.5*length(NSignal)));
        SignalMat(count2) = Signal;
        count2 = count2 + 1;
    else
        Double = ((Y(:,[14:end-2])*multipyer([14:end]).')*10^9); %only inline and twisted
        p(1) = plot(T,Double,'LineWidth',2, 'Color', [0 0.4470 0.7410 1-(1/11)*(count2-1)]); %%
        hold on
        p(2) = plot(T(round(0.5*length(T))), Double(round(0.5*length(Double))), 'ko'); %%
        set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %%
        ylabel('Signal (nM)'); xlabel('Time (s)') %%
        SignalMat(count2) = Double(round(0.5*length(Double))); %Double(round(0.5*length(Double)))
        Sig = Double(round(0.5*length(Double)));
        count2 = count2 + 1;
    end
    %plot(T,Signal3, 'LineWidth',2)
    %hold on
    %plot(T, twisted, 'LineWidth',2)
    %hold on
    %plot(T(round(0.5*length(T))), Nonly, 'ko')
    %hold on
    %p(2) = plot(T(round(0.5*length(T))), Double(round(0.5*length(Double))), 'ko'); %%
    %set(get(get(p(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %%
    %plot(T(round(.5*length(T))), Signal(round(.5*length(T))), 'ko')
    %title('Signal')
    %legend('a1a2', 'a1x0', 'a2a1')
    %Signal3(round(0.5*length(Signal3)))
    %Signal = Double(round(0.5*length(Double)))
    %Signal = NSignal(round(0.5*length(NSignal)))
elseif plottingtype==2 && Multiple~=1
    %Lnorm = Y(:,16);
    %ratio = Kons11/Lnorm
    curvemax=max(Y);
    names(1)={'Signal'};
    count = 2; %new
    for i = 3:(length(KonTable))
        if i < 7 || i > 14
            names(count) = {num2str(KonTable{1,i})};
            count = count+1;
        end
    end
%     for i=3:(length(KonTable))
%         names(i-1)={num2str(KonTable{1,i})};
%     end
    %plot(T,Signal,T,Y(:,2:end-2).*multipyer(2:end)*RUmod,'LineWidth',2) %original

    B = multipyer; %new
    for i=1:length(Y(:,2:end-2))-1 %new
        multipyer = [multipyer;B]; %new
    end %new
    %plot(T,Signal,T,Y(:,2:end-2).*multipyer(:,2:end)*RUmod,'LineWidth',2)%alllines
    plot(T,Signal,T,Y(:,[2:3,14]).*multipyer(:,[2:3,14])*10^9,'LineWidth',2)%noB
    %hold on
    %plot(T,Signal,T,Y(:,[14:15])*multipyer([14:end].')*RUmod,'LineWidth',2)%new
%     plot(T,Signal1)
%     hold on
%     plot(T,Signal2)
%     hold on
    %plot(T,Y(:,[2,14:end-2]).*multipyer(:,[2,14:end]).*RUmod) %plots N-only and both tandem
%     inline = (Y(:,14).*multipyer(:,14).*RUmod);
%     twisted = (Y(:,15).*multipyer(:,15).*RUmod);
%     Rtot = R0.*multipyer(:,14).*RUmod;
%     FreeL = Rtot - (inline+twisted);
    %plot(T,Y(:,[14:end-2]).*multipyer(:,[14:end]).*RUmod) %plots both tandem and free ligand?
    %hold on
    %plot(T,FreeL)
    %plot(T, Signal)
    %hold on
    %plot(T(round(0.5*length(T))), Signal(round(0.5*length(Signal))), 'ko')
    %plot(T(round(0.5*length(T))), FreeL(round(0.5*length(FreeL))), 'ko')
    %plot(T(135), Signal(135), 'ko')
    %plot(T,Y(:,[2:5]).*multipyer(:,[2:5]).*RUmod,'LineWidth',2) %noBorDouble
    xlabel('Time (s)')
    ylabel('RU')
    %legend('a1a2', 'a2a1', '[L]')
    %lgd = legend(names); %all legend names
    lgd = legend('Total bound','N-only','Double')
%     lgd.FontSize = 12;
%     lgd.FontName = 'Arial';
%     %Signal(75)
%     Signal(round(0.5*length(Signal)));
    %FreeL(round(0.5*length(FreeL)))
    %Signal1 = Signal1(round(0.5*length(Signal1)));
    %Signal2 = Signal2(round(0.5*length(Signal2)));
    %C = (Y(round(0.5*length(Y)),14)*multipyer(14.')*RUmod)+(Y(round(0.5*length(Y)),15).*multipyer(14.')*RUmod)
    %legend(names(2:5)) %no B or Double
    
%     Assoc = Signal(1:length(TA));
%     Dissoc = Signal(length(TA)+1:end);
    %plot(TD, Dissoc)
    %Signal = NSignal(round(0.5*length(NSignal)))
elseif plottingtype==2 && Multiple==1 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   a=[1,2,3,4,5,6,7,8,9];
%   a(2:border(1)),
%   a(border(1):border(2)),
%   a(border(2):length(KonTable)-1),

    
    names(1)={'Signal'};
    for i=3:(length(KonTable))
        names(i-1)={num2str(KonTable{1,i})};
    end
    legend(names)
    
elseif plottingtype==1 && Multiple==1
    if length(border)==2
        Signal=sum(Y(:,2:border(1)),2)*RUmod+sum(Y(:,border(1)+1:border(2)),2)*RUmod*2 + sum(Y(:,border(2)+1:end-2),2)*RUmod*3;
        plot(T,Signal,'LineWidth',2)
        ylabel('RU'); xlabel('Time (s)')
    elseif length(border)==1
        Signal=sum(Y(:,2:border(1)),2)*RUmod+sum(Y(:,border(1)+1:end-2),2)*RUmod*2;
        plot(T,Signal,'LineWidth',2)
        ylabel('RU'); xlabel('Time (s)')        
    end
end
    


end