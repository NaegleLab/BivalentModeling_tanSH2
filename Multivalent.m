function Multivalent(Nrec,Mlig, BaseL)%symtype
global w
for w = 1:1 %*************
    %%
    % Function description:
    % The function is the main frame of the multivalent kinetic calculation
    % Inputs:   (symtype,Nrec,Mlig,BaseL)
    %           -symtype use 1 in our model (not used currently)
    %           -Nrec: receptor number
    %           -Mlig: ligand number
    %           -BaseL or L0 is the ligand concentration (plotted it's x(0.1-10) fold concentrations)
    
    % Outputs: KonTable and the plots
    % Used functions: All the function
    
    %% 1) The initial parameters
    
    
    %% 1.1) The experiment type choice
    
    %Global parameters for the experiment may used by other functions
    global Kon Koff Kon2 Koff2 Kon3 Koff3 Kon4 Koff4 R0 L0 KMTL OnTime OffTime N M MWreceptor MWligand baseRU Multiple lgd
    global KonTable KoffTable               %The State contant tables as global variables
    
    global LinkerNLig UnitDiameterLig LinkerNRec UnitDiameterRec
    global rerun type plottype
    global count2 SignalMat%*****************
    
    [M,N,rerun,type,plottype,MWreceptor,MWligand,LinkerNLig,UnitDiameterLig,LinkerNRec,...
        UnitDiameterRec,Kon,Koff,Kon2,Koff2,Kon3,Koff3,Kon4,Koff4,OnTime,OffTime,KMTL,L0,baseRU] = parameters;
    
    if isempty(rerun); rerun=1; end
    if isempty(type); type=1; end
    if isempty(plottype); plottype=2; end %2
    
    Multiple=0;
    
    %% 1.2) Input managemant
    
    if nargin==3
        N=Nrec;M=Mlig; L0=BaseL;
    end
    
    %%
    
    
    switch type
        %In case of our model also the default setp
        
        
        %% Trap model
        case 1
            if isempty(N); N=2; end
            if isempty(M); M=2;  end
            
            if isempty(MWreceptor)
                switch N
                    case 1
                        MWreceptor=33000;
                    case 2
                        MWreceptor=41000;
                    case 3
                        MWreceptor=49000;
                end
            end
            
            if isempty(MWligand)
                switch M
                    case 1
                        MWligand=18000;
                    case 2
                        MWligand=21000;
                    case 3
                        MWligand=24000;
                end
            end
            
            %Linker and unit parameters
            if isempty(LinkerNLig); LinkerNLig=23; end
            if isempty(UnitDiameterLig); UnitDiameterLig=11*3; end
            if isempty(LinkerNRec); LinkerNRec=23; end
            if isempty(UnitDiameterRec); UnitDiameterRec=22; end
            
            %Kinetic parameters
            
            if isempty(Kon); Kon=600000; end
            if isempty(Koff); Koff=0.3;  end
            if isempty(Kon2); Kon=600000; end
            if isempty(Koff2); Koff=0.3;  end
            if isempty(Kon3); Kon=600000; end
            if isempty(Koff3); Koff=0.3;  end
            if isempty(Kon4); Kon=600000; end
            if isempty(Koff4); Koff=0.3;  end
            if isempty(OnTime); OnTime=200; end
            if isempty(OffTime); OffTime=240;  end
            if isempty(KMTL); KMTL=1;  end
            if isempty(L0)
                if nargin==3; L0=BaseL;
                else;L0=1e-6;  end
            end
            if isempty(baseRU); baseRU=60; end
            
            
            %% In case of the DNA model---------------
        case 2
            if isempty(N); N=2; end
            if isempty(M); M=2;  end
            if isempty(LinkerNLig); LinkerNLig=23; end
            if isempty(UnitDiameterLig); UnitDiameterLig=11*3; end
            if isempty(LinkerNRec); LinkerNRec=23; end
            if isempty(UnitDiameterRec); UnitDiameterRec=22; end
            
            %Kinetic parameters
            
            if isempty(Kon); Kon=400; end
            if isempty(Koff); Koff=0.32;  end
            if isempty(OnTime); OnTime=90; end
            if isempty(OffTime); OffTime=240;  end
            if isempty(KMTL); KMTL=1;  end
            if isempty(L0); L0=1e-5; end
            if isempty(baseRU); baseRU=2700; end
            
            %The initial RU of the receptor
        case 3
            %Does not need any parameters because they were already given previously as
            %global parameters
    end
    
    %Parameters for RU->concentration convertion-------------------------------
    global Area Vol Ruconv
    Area=1.6*1e-4;
    Vol=1.6e-10;
    Ruconv=1e-12;
    
    R0=baseRU*Ruconv*Area/(MWreceptor*Vol);
    RUmax=baseRU*MWligand/MWreceptor
    
    global MTL FlowRate
    FlowRate=50; %%%%%%%%%%%%%%%%
    flow=(FlowRate/10)/100;
    MTL=0.98*((10^(-0.434*log10(MWligand)-4.059)*10^-4)/(100*10^-9))^(2/3)*(flow/(0.2963*1.6*10^-6))^(1/3); %Determination of diffusion coefficients of peptides and prediction of permeability through a porous membrane
    
    
    
    %% 2) Calculating all the states and the Kon and Koff tables
    global multipyer
    %[KonTable, KoffTable]=RateConstantTable(N, M);
    [KonTable, KoffTable, multipyer]=StatesAndTables(N, M);
    %By using the DiffStates_Kon function it calculates the table values
    %type: 1-our model, 2-DNA model, 3-fullsymmetric,  4-partially symmetric 5-not symmetric
    
    if Multiple==0
    elseif Multiple==1
        if N==2 && M==1
            load('Table21.mat')
        end
        if N==2 && M==2
            load('Table22.mat')
        end
        if N==3 && M==1
            load('Table31.mat')
        end
    else
        save(strcat('T',num2str(N),num2str(M)),'KonTable', 'KoffTable')
    end
    
    
    %% 3) %Create the differential equations and write tham to file
    if rerun~=0
        DiffEquations('Diffs'); %The association phasewith L0 initial ligand CC
        %It also uses Effective concentration script to calculate the rate constants
    end
    if N==1 && M==1
        load('Mix1x1.mat')
    elseif N==2 && M==2
        load('Mix2x2.mat')
    end
    
    %% 4) Plotting
    fig1 = figure;
    switch plottype
        case 1
            title(['Ligand concentration series with L0=',num2str(L0)])
            Plotting(L0)
            %series=[1,2,8,32,128,512,2028];
            series=[1,2,4,8,16,32,64,128,256,512,1024];
            %series=[1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65356];
            names(1)={strcat('L0: ',num2str(L0*10^9,4),' (nM)')};
            for i=2:(length(series))
                names(i)={num2str((L0*10^9)/series(i),4)};
            end
            %lgd = legend(names);
                    % set(gcf,'PaperUnits','inches');
                    % set(gcf,'PaperSize',[8 6]);
                    % %set(gca,'FontSize',12,'FontName','Arial')
                    % a = get(gca,'XTickLabel');
                    % set(gca,'XTickLabel',a,'FontName','Arial','FontSize',10)
                    % xlabel('Time (s)','FontSize',16,'FontName','Arial')
                    % ylabel('Signal (nM)','FontSize',16,'FontName','Arial')
                    % %lgd.FontSize = 10;
                    % %lgd.FontName = 'Arial';
                    % print(fig1,'SerialDilutionNEW.pdf','-dpdf')
            hold off
            %conc = [240000, 120000, 60000, 30000, 15000, 7500, 3750, 1875, 937.5, 468.75, 234.375];
            %conc = [120000, 60000, 30000, 15000, 7500, 3750, 1875, 937.5, 468.75, 234.375, 117.1875];
            %conc = [60000, 30000, 15000, 7500, 3750, 1875, 937.5, 468.75, 234.375, 117.1875, 58.59375];
            %conc = [30000, 15000, 7500, 3750, 1875, 937.5, 468.75, 234.375, 117.1875, 58.59375, 29.296875];
            %conc = [16000, 8000, 4000, 2000, 1000, 500, 250, 125, 62.5, 31.25, 15.625];
            %conc = [8000, 4000, 2000, 1000, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125]; %nM
            %conc = [4000, 2000, 1000, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625]; %nM
            %conc = [2000, 1000, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625, 1.953125];
            %conc = [1000, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625, 1.953125, .9765625];
            %conc = [400, 200, 100, 50, 25, 12.5, 6.25, 3.125, 1.5625, .78125, .390625]
            conc = [250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625, 1.953125, .9765625, .48828125, .24414062];
            %conc = [125, 62.5, 31.25, 15.625, 7.8125, 3.90625, 1.953125, .9765625, .48828125, .24414062, .12207031];
            %conc = [60, 30, 15, 7.5, 3.75, 1.875, .9375, .46875, .234375, .1171875, .05859375];
            %conc = [31.25, 15.625, 7.8125, 3.90625, 1.953125, .9765625, .48828125, .24414062, .12207031, .061, .0305];
            %conc = [15.625, 7.8125, 3.90625, 1.953125, .9765625, .48828125, .24414062, .12207031, .06103516, .03051758, .01525879];
            %conc = [60000, 30000, 15000, 7500, 3750, 1875, 937.5, 468.75, 234.375, 117.1875, 58.59375, 29.296875,...
            %    14.6484375, 7.32421875, 3.66210938, 1.83105469, .9155];
            x = conc;
            %signal = flip(signal);
            y = SignalMat'./(R0*10^9); %using L = free ligand (SignalMat = nM, R0 = M)
            %halfMax = max(y)/2;
            halfMax = 0.5;
            [curve, goodness, output] = fit(x',y','smoothingspline');
            %plot(curve)
            %legend('hide')
            %hold on
            %axis([0 125 0 1])
            %curve2 = zeros(1e5,1);
            span = linspace(conc(1)/1e5,conc(1),1e5);
%             for z = 1 : 1e5
%                 if curve(z) <= halfMax
%                     curve2(z) = 0;
%                 else
%                     curve2(z) = 1;
%                 end
%             end
            for z = span
                if curve(z) > halfMax
                    break
                end
            end
            %xHigh = min(find(curve2));
            xHigh = z;
            if xHigh == conc(1)
                bad = 1+1
            end
            %xLow = xHigh-1;
            xLow = xHigh - (conc(1)/1e5);
            Kd = (xHigh + xLow)/2
            KA = 1/Kd
           % plot(Kd, (curve(xHigh)+curve(xLow))/2, 'k*')
           % hold on
           % plot([0:1:10]', 0.5*ones(1,11), 'k--')
           % hold on
           % plot(10.9305*ones(1,6), [0:.1:.5]', 'k--')
           % set(gcf,'PaperUnits','inches');
           %  set(gcf,'PaperSize',[8 6]);
           %  set(gca,'FontSize',12,'FontName','Arial')
           %  a = get(gca,'XTickLabel');
           %  set(gca,'XTickLabel',a,'FontName','Arial','FontSize',10)
           %  xlabel('Concentration (nM)', 'FontSize',16,'FontName','Arial')
           %  ylabel('Fraction Bound', 'FontSize',16,'FontName','Arial')
           %  %lgd.FontSize = 12;
           %  %lgd.FontName = 'Arial';
           %  print(fig1,'ConcEffectCurve.pdf','-dpdf')
        case 2
            Plotting(L0)
            hold off
        case 3
            Plotting(L0)
            %set(gcf,'PaperUnits','inches');
            %set(gcf,'PaperSize',[4 3]);
            %set(gca,'FontSize',12,'FontName','Arial')
            %a = get(gca,'XTickLabel');
            %set(gca,'XTickLabel',a,'FontName','Arial','FontSize',10)
            %xlabel('Time (s)', 'FontSize',16,'FontName','Arial')
            %ylabel('Signal (nM)', 'FontSize',16,'FontName','Arial')
            %lgd.FontSize = 12;
            %lgd.FontName = 'Arial';
            %print(fig1,'ExModel.eps','-depsc')
            hold off
        case 4
            sttext=['Ligand conc:=',num2str(L0), 'times 1000 to 0.001','With MTL'];
            suptitle(sttext)
            %KMTL=1;
            for i=1:9
                Koff=1e-5*10^i;
                Kon=Koff/(500*1e-9);
                subplot(3,3,i)
                Plotting(L0)
                title(['Kon=',num2str(Kon),', Koff=',num2str(Koff)])
            end
            legend('L0','L0*3','L0*10','L0*30','L0*100','L0*1000','L0/3','L0/10','L0/30','L0/100','L0/1000')
    end
    
    %KDcalculator(RUmax)
end %******************
end