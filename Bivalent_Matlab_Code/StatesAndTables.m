function [TableOn, TableOff, multipyer] = StatesAndTables(N,M)
%% Showing all the possible ligand names
% Longer for transparency

if N==3
    if M==3
        lig={'a1','a2','a3','b1','b2','b3','c1','c2','c3','x0','y0','z0'};
    elseif M==2
        lig={'a1','a2','b1','b2','c1','c2','x0','y0','z0'};
    elseif M==1
        lig={'a1','b1','c1','x0','y0','z0'};
    end
elseif N==2
    if M==3
        lig={'a1','a2','a3','b1','b2','b3','x0','y0'};
    elseif M==2
        lig={'a1','a2','b1','b2','x0','y0'};
    elseif M==1
        lig={'a1','b1','x0','y0'};
    end
elseif N==1
    if M==3
        lig={'a1','a2','a3','x0'};
    elseif M==2
        lig={'a1','a2','x0'};
    elseif M==1
        lig={'a1','x0'};
    end
end

%% Creating all the unique states (taking into consideration identities and symmetries)

states={};
counter=1;

for i = lig
    for j = lig
        if prod(j{1}==i{1})==0
           %j{1}~=i{1},  j{1},i{1}
            for k = lig
                if prod(k{1}==j{1})==0
                    if prod(k{1}==i{1})==0
                        if N==3
                            states(counter,1)={[i{1},j{1},k{1}]};
                        elseif N==2
                            states(counter,1)={[i{1},j{1}]};
                        elseif N==1
                            states(counter,1)={[i{1}]};
                        end                        
                    counter=counter+1;
                    end
                end

            end

        end
    end
end

for i=1:length(states)
    states{i,1}(states{i}=='y')='x';
    states{i,1}(states{i}=='z')='x';    
end

Allstates=unique(states);
if N==1 && M==1
    Allstates(1,1)=lig(1);Allstates(2,1)=lig(2);
end
states=Allstates;

states(1,2)=states(1,1);
for i =2:length(Allstates)
    states(i,2) = {StateMapper(N, M, Allstates{i,1}, states(1:i-1,2))};
end

Allstates=states;
Uniquestates=unique(Allstates(:,2));

%%
counter =1;
for step=[3,2,1,0]
    
for i =1:length(Uniquestates)
    if sum(Uniquestates{i}=='x')==step
        Uniquestatesorderd(counter)=Uniquestates(i);
        counter=counter+1;
    end    
end
end
Uniquestates=Uniquestatesorderd;
for i = 1:length(Uniquestates)
    if sum(Uniquestates{i}=='c')==1
        multipyer(i)=3;
    elseif sum(Uniquestates{i}=='b')==1
        multipyer(i)=2;
    elseif sum(Uniquestates{i}=='a')>=1
        multipyer(i)=1;
    else
        multipyer(i)=0;
    end
end
%%
TableOn=cell(length(Uniquestates)+1);   %Creates the empty cell matrix
TableOn(1,2:end)=Uniquestates;          %Gives the row names
TableOn(2:end,1)=Uniquestates;          %Gives the column names
TableOff=TableOn;                       %Creates the empty Kofftable

statemap=containers.Map;
rankmap=containers.Map;
for i=1:length(Allstates)
    statemap(Allstates{i,1})=Allstates{i,2};
end
for i=1:length(Uniquestates)
    rankmap(Uniquestates{i})=i+1;
end

for i = 1:length(Allstates)
    associated=Allstates{i,1};
    posass=rankmap(statemap(associated));
    for j = 1:length(Allstates)
        disociated=Allstates{j,1};
        posdiss=rankmap(statemap(disociated));
        if length(disociated(disociated~=associated))==2
            if prod(disociated(disociated~=associated)=='x0')
                [Kon, Koff]=findKon(associated, disociated, M, N);
                TableOn(posdiss,posass)={Kon};
                TableOff(posass,posdiss)={Koff};
            end
        end
    end
end
%%
end

function [Kon, Koff] = findKon(ass, diss, M, N)

switch max([sum(ass=='a'),sum(ass=='b'),sum(ass=='c')]) % Number of different ligands on the receptor
    case 3
            Kon=findKonValency([str2num(ass(2)),str2num(ass(4)),str2num(ass(6))], [str2num(diss(2)),str2num(diss(4)),str2num(diss(6))]);
            
    case 2
        if max([sum(diss=='a'),sum(diss=='b'),sum(diss=='c')])==2
            Kon='1*Kons00';
           
        else
            if N==2
                Kon=findKonValency([str2num(ass(2)),str2num(ass(4))], [str2num(diss(2)),str2num(diss(4))]); 
                
            else
                if sum(ass=='a')==2
                    ass(max([find(ass=='c'),find(ass=='b')])+1)=0;
                    diss(max([find(ass=='c'),find(ass=='b')])+1)=0;
                elseif sum(ass=='b')==2
                    ass(max([find(ass=='c'),find(ass=='a')])+1)=0;
                    diss(max([find(ass=='c'),find(ass=='a')])+1)=0;
                elseif sum(ass=='c')==2
                    ass(max([find(ass=='a'),find(ass=='b')])+1)=0;
                    diss(max([find(ass=='a'),find(ass=='b')])+1)=0;
                end
                Kon=findKonValency([str2num(ass(2)),str2num(ass(4)),str2num(ass(6))], [str2num(diss(2)),str2num(diss(4)),str2num(diss(6))]); 
                
            end
            
        end
    case 1
        switch N
            case 3
                if sum(diss=='x')==3
                    if prod(find(ass~=diss)==[3,4])==1
                        if M==3
                            if ass(4)=='1' || ass(4)=='3' 
                                Kon='2*Kons00'; 
                            elseif ass(4)=='2'
                                Kon='1*Kons00';
                            else
                                Kon='wut'
                            end
                        elseif M==1
                            Kon='1*Kons00';
                        else
                           Kon='2*Kons00'; 
                        end
                    else
                        Kon='2*Kons00';
                    end
                    
                elseif M==1
                    if prod([diss(2),diss(6)]==['0','0'])==1
                        Kon='2*Kons00';
                    else
                        Kon='1*Kons00';
                    end
                elseif M==2 %The M==1 above is for N3M1 correction is this ok?
                    Kon='1*Kons00';
                elseif M==3
                    Kon='1*Kons00';
                end
            case 2
                if sum(diss=='x')==2
                    Kon='2*Kons00';
                else
                    Kon='1*Kons00';
                end
            case 1
                if M==3
                    if ass(2)=='2'
                        Kon='1*Kons00';
                    else
                        Kon='2*Kons00';
                    end
                else
                Kon=strcat(num2str(M),'*Kons00');
                end
        end
        

end

        %% Koff
        
        if N==3
            if prod(ass([2,6])~='0')==1 % non of the ends are empty
                if M==1
                    if prod(find(diss~=ass)==[3,4])==1
                        Koff='1*Koff';
                    else
                        Koff='2*Koff';
                    end
                elseif M==2
                    if prod(find(ass~=diss)~=[3,4])==1
                        if ass(2)~=ass(6) && ass(1)==ass(5) && ass(4)=='0'
                            Koff='2*Koff';
                        else
                            Koff='1*Koff';
                        end
                    else
                       Koff='1*Koff'; 
                    end
                elseif M==3
                     if prod(find(ass~=diss)~=[3,4])==1
                        if prod([ass(2),ass(6)]=='22')==1 || prod([ass(2),ass(6)]=='31')==1 || prod([ass(2),ass(6)]=='13')==1
                            Koff='2*Koff';
                        else
                            Koff='1*Koff';
                        end
                    else
                       Koff='1*Koff'; 
                    end 
                end
            else
                Koff='1*Koff';
            end
            
        elseif N==2
            if M==2
                if prod(diss([2,4])=='0')==1
                    Koff='1*Koff';
                elseif ass(2)~=ass(4)
                    Koff='2*Koff';
                else
                    Koff='1*Koff'; 
                end
            elseif M==3
                if prod(diss([2,4])=='0')==2
                    Koff='1*Koff';
                elseif prod(ass([2,4])=='13')==1 || prod(ass([2,4])=='31')==1 || prod(ass([2,4])=='22')==1
                    Koff='2*Koff';
                else
                    Koff='1*Koff';
                end
            elseif M==1
                if prod(diss([2,4])=='0')==1
                    Koff='1*Koff';
                else
                    Koff='2*Koff';
                end
            end
        else 
            Koff='1*Koff';
        end
end

function [MappedState] = StateMapper(N, M, state, UniqueStates)


        %"a" should be the only duplicated state
        if sum(state=='b')==2 || sum(state=='b')==3
            
            state(state=='b')='g';
            state(state=='a')='b';
            state(state=='g')='a';
            
        end
        if sum(state=='c')==2 || sum(state=='c')==3
            state(state=='c')='g';
            state(state=='a')='c';
            state(state=='g')='a';
        end
        %there must be always an "a" state
        if sum(state=='b')==1 && sum(state=='a')==0
            state(state=='b')='a';
        elseif sum(state=='c')==1 && sum(state=='a')==0
            state(state=='c')='a';
        end
        %"b" instead of "c" 
        if sum(state=='c')==1 && sum(state=='b')==0
            state(state=='c')='b';
        end

        %for "abc" states the order must be "a", "b", "c"
        if sum(state=='a')==1 && sum(state=='b')==1 && sum(state=='c')==1
            state(1)='a';
            state(3)='b';
            state(5)='c';
        end

        %for "ab" states the order must be "a", "b"
        if sum(state=='a')==1 && sum(state=='b')==1 && sum(state=='c')==0
            if find(state=='a')>find(state=='b')
                state(state=='a')='g';
                state(state=='b')='a';
                state(state=='g')='b';
            end
        end
        %%%

%%
if N==3
        if sum(state=='a')==1 && sum(state=='b')==1 && sum(state=='c')==1
            newstate=['a',state(6),'b',state(4),'c',state(2)];
            newstate([2,4,6])=[num2str(M+1-str2num(newstate(2))),num2str(M+1-str2num(newstate(4))),num2str(M+1-str2num(newstate(6)))];
        elseif sum(state=='a')==1 && sum(state=='b')==1 && sum(state=='c')==0
            newstate=[state(5:6),state(3:4),state(1:2)];
            newstate([2,4,6])=[num2str(M+1-str2num(newstate(2))),num2str(M+1-str2num(newstate(4))),num2str(M+1-str2num(newstate(6)))];
            newstate(newstate==num2str(M+1))='0';
                newstate(newstate=='a')='g';
                newstate(newstate=='b')='a';
                newstate(newstate=='g')='b';
        else
            newstate=[state(5:6),state(3:4),state(1:2)];
            newstate([2,4,6])=[num2str(M+1-str2num(newstate(2))),num2str(M+1-str2num(newstate(4))),num2str(M+1-str2num(newstate(6)))];
            newstate(newstate==num2str(M+1))='0';
        end
        %[states{i},' ',newstate]
        for j=1:length(UniqueStates)
            if prod(UniqueStates{j}==newstate)==1
                state=UniqueStates{j};
            end  
        end
elseif N==2
        if sum(state=='a')==1 && sum(state=='b')==1
            newstate=[state(3:4),state(1:2)];
            newstate([2,4])=[num2str(M+1-str2num(newstate(2))),num2str(M+1-str2num(newstate(4)))];
            newstate(newstate==num2str(M+1))='0';
                newstate(newstate=='a')='g';
                newstate(newstate=='b')='a';
                newstate(newstate=='g')='b';
        else
            newstate=[state(3:4),state(1:2)];
            newstate([2,4])=[num2str(M+1-str2num(newstate(2))),num2str(M+1-str2num(newstate(4)))];
            newstate(newstate==num2str(M+1))='0';
        end
        %[states{i},' ',newstate]
        for j=1:length(UniqueStates)
            if prod(UniqueStates{j}==newstate)==1
                state=UniqueStates{j};
            end  
        end
elseif N==1        
            newstate=state;
            newstate(2)=num2str(M+1-str2num(newstate(2)));
            newstate(newstate==num2str(M+1))='0';
        for j=1:length(UniqueStates)
            if prod(UniqueStates{j}==newstate)==1
                state=UniqueStates{j};
            end  
        end
end
%%        
        MappedState=state;
end

function Kon=findKonValency(ass, diss)
    
    new_ligand=     max(ass-diss);                                          %The value of the associating ligand
    min_ligand=     min(diss(diss~=0));                                     %the value of the smallest previous ligand
    max_ligand=     max(diss);                                              %the value of the largest previous ligand
    
    new_location=   find(ass==new_ligand);                                  %location of the associating ligand
    min_location=   find(ass==min_ligand);                                  %location od the previous smallest ligand
    max_location=   find(ass==max_ligand);                                  %location of the largest previous ligand
    
    
    %----------------------------------------------------------------------The ligand values
    Kon_ligandX=0;                                                          %variable for bint ends
    if new_ligand>min_ligand && new_ligand<max_ligand                       %Check for both ends
        Kon_ligand=(10*(max_ligand-new_ligand))+(new_ligand-min_ligand);    %creates a two number digit e.g. (KonX)12(X11)
        Kon_ligandX=1;                  
    elseif new_ligand>max_ligand                                            %If the new ligand is larger then max gets the value
        Kon_ligand=new_ligand-max_ligand;
    elseif new_ligand<min_ligand                                            %If the new ligand is smaller then min gets the value
        Kon_ligand=min_ligand-new_ligand;    
    end
    
    %----------------------------------------------------------------------
    Kon_locationX=0; %itt a ketto kozott kell lennie
    if min_location>max_location
        reverse=min_location;
        min_location=max_location;
        max_location=reverse;
        reverse=1;
    else 
        reverse=0;
    end
    
    %new_location, min_location, max_location, ass, diss
    if new_location>min_location && new_location<max_location      %????????????? & &&         %Check for bint location
        Kon_location=(10*(max_location-new_location))+(new_location-min_location);  %bi values (KonX11X)23
        Kon_locationX=1;
    elseif new_location>max_location                                        %Location value if the location is large then the previous
        Kon_location=new_location-max_location;
    elseif new_location<min_location                                        %Location number if the location is smaller then the previous
        Kon_location=min_location-new_location;    
    end  
    
    %----------------------------------------------------------------------The symmetry value calculation
                                    %in case of partial symmetry
            
            if Kon_ligandX==1 && Kon_locationX==1                           %in case of double XX%reverse or not
                if reverse==1
                    Kon_type='r'; 
                else
                    Kon_type='s';
                end 

            elseif Kon_ligandX==1 || Kon_locationX==1                       %in case of simple X
                Kon_type='s';                                               %there is no reverse form

            else % Kon_ligandX==0 && Kon_locationX==0                       %in case of no X
                if new_ligand<min_ligand && new_location<min_location %????????????? & &     %if the ligand and location are both smaller = simple
                        
                    if reverse==1
                        Kon_type='r'; 
                    else
                        Kon_type='s';
                    end 

                elseif new_ligand>max_ligand && new_location>max_location %??????????? & && & %if the ligand and location are both larger = simple
                    if reverse==1
                        Kon_type='r'; 
                    else
                        Kon_type='s';
                    end 
                else                                                        %if either the ligand or location is reverse = reverse
                    if reverse==1
                        Kon_type='s'; 
                    else
                        Kon_type='r';
                    end 
                end
            end
            
   
    if Kon_locationX==1
        Kon_location=['X',num2str(Kon_location)];
    else
        Kon_location=num2str(Kon_location);
    end
    
    if Kon_ligandX==1
        Kon_ligand=['X',num2str(Kon_ligand)];
    else
        Kon_ligand=num2str(Kon_ligand);
    end
    
    Kon=['1*Kon',Kon_type, Kon_ligand, Kon_location];

end