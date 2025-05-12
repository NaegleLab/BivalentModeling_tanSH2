function [] = DiffEquations(name)
global KonTable KoffTable            %Uses the global tables
global N M

delete Diffs.m
dir
diff = fopen(strcat(name,'.m'),'w');    %Creates function for Diffs
B=KonTable;
C=KoffTable;

%fprintf(diff, 'import mlreportgen.dom.*\n');
%fprintf(diff, 'd = Document(''test'', ''html'');\n');
%fprintf(diff, 'ol = UnorderedList();\n');
fprintf(diff, 'function dydt = %s(t,y) \n\n', name);    %Initiate the function
fprintf(diff, 'global Kon Koff Kon2 Koff2 Kon3 Koff3 Kon4 Koff4 L R0 KMTL w N M MTL Kons11 ratioTable currentRow\n');   %Start with the global var Kon and Koff 
fprintf(diff, 'km=MTL*KMTL;\n'); 
fprintf(diff, 'KonA11=Kon(w);\n');
fprintf(diff, 'KonA12=Kon2(w);\n');
fprintf(diff, 'KonA21=Kon3;\n');
fprintf(diff, 'KonA22=Kon4;\n');
fprintf(diff, 'KoffA11=Koff;\n');
fprintf(diff, 'KoffA12=Koff2;\n');
fprintf(diff, 'KoffA21=Koff3;\n');
fprintf(diff, 'KoffA22=Koff4;\n');

Konlist=KonTable(2:end,2:end); %%remove headers?
Konlist=Konlist(~cellfun(@isempty, Konlist)); %%removes empty cells
for i=1:length(Konlist)
    Konlist(i)={Konlist{i}(3:end)}; %%moves rows into new list
end
Konlist=unique(Konlist); %%keeps unique rows and sorts

for i=1:length(Konlist)
    
            if prod(Konlist{i}(end-1:end)=='00')
            fprintf(diff, '%s=Kon; \n', Konlist{i});    
            else
            Leff = EffC_Calculator(Konlist{i});
            fprintf(diff, '%s=%d; \n', Konlist{i}, (Leff));% Kon->(EffC_Calculator(Konlist{i}))
            %**************************************************************************************************
            %**Adding RO to the Kon represents the chance of binding another receptor on the surface***
            %**************************************************************************************************
            end
end

%fprintf(diff,'ratio = Kons11/y(16)\n')
%fprintf(diff,'ratioTable(currentRow) = ratio;\n')
%fprintf(diff,'currentRow = currentRow + 1;\n')
%fprintf(diff,'append(ol,ratio);\n')

diffs11=['dydt = zeros(4,1);\n', ...
'dydt(1)=+1*Koff*y(2)-1*Kons00*y(1)*y(4)-1*Kon2*y(1)*y(5)+1*Koff2*y(3);\n', ...
'dydt(4)=km*L-km*y(4)+1*Koff*y(2)-1*Kons00*y(1)*y(4);\n', ...
'dydt(5)=km*L-km*y(5)+1*Koff2*y(3)-1*Kon2*y(1)*y(5);\n', ... %Switched L2 --> L
'\n', ...
'dydt(2)=+1*Kons00*y(1)*y(4)-1*Koff*y(2);\n', ...
'dydt(3)=+1*Kon2*y(1)*y(5)-1*Koff2*y(3);\n'];

diffs22=['dydt = zeros(17,1);\n', ...
'dydt(1)=y(2)*KoffA11+y(3)*KoffA12+y(4)*KoffA21+y(5)*KoffA22-y(1)*y(16)*KonA11-y(1)*y(16)*KonA12-y(1)*y(16)*KonA21-y(1)*y(16)*KonA22;\n', ...
'dydt(16)=km*L-km*y(16)+y(2)*KoffA11+y(3)*KoffA12+y(4)*KoffA21+y(5)*KoffA22+y(6)*KoffA12+y(6)*KoffA11+y(7)*KoffA22+y(7)*KoffA11+y(8)*KoffA12+y(8)*KoffA11+y(9)*KoffA21+y(9)*KoffA12+y(10)*KoffA21+y(10)*KoffA12+y(11)*KoffA22+y(11)*KoffA21+y(12)*KoffA22+y(12)*KoffA11+y(13)*KoffA22+y(13)*KoffA21-y(1)*y(16)*KonA11-y(1)*y(16)*KonA12-y(1)*y(16)*KonA21-y(1)*y(16)*KonA22-y(2)*y(16)*KonA12-y(2)*y(16)*KonA22-y(3)*y(16)*KonA11-y(3)*y(16)*KonA21-y(4)*y(16)*KonA12-y(4)*y(16)*KonA22-y(5)*y(16)*KonA11-y(5)*y(16)*KonA21;\n', ...
'dydt(17)=y(16)*0.01;\n', ...
'dydt(2)=y(1)*y(16)*KonA11-y(16)*y(2)*KonA12-y(16)*y(2)*KonA22-y(2)*Kons11*KonA22+y(6)*KoffA12+y(7)*KoffA22+y(8)*KoffA12+y(12)*KoffA22+y(14)*KoffA22-y(2)*KoffA11;\n', ...
'dydt(3)=y(1)*y(16)*KonA12-y(16)*y(3)*KonA11-y(16)*y(3)*KonA21-y(3)*Konr11*KonA21+y(6)*KoffA11+y(8)*KoffA11+y(9)*KoffA21+y(10)*KoffA21+y(15)*KoffA21-y(3)*KoffA12;\n', ...
'dydt(4)=y(1)*y(16)*KonA21-y(16)*y(4)*KonA12-y(16)*y(4)*KonA22-y(4)*Konr11*KonA12+y(9)*KoffA12+y(10)*KoffA12+y(11)*KoffA22+y(13)*KoffA22+y(15)*KoffA12-y(4)*KoffA21;\n', ...
'dydt(5)=y(1)*y(16)*KonA22-y(16)*y(5)*KonA11-y(16)*y(5)*KonA21-y(5)*Kons11*KonA11+y(7)*KoffA11+y(11)*KoffA21+y(12)*KoffA11+y(13)*KoffA21+y(14)*KoffA11-y(5)*KoffA22;\n', ...
'dydt(6)=y(2)*y(16)*KonA12-y(6)*KoffA12-y(6)*KoffA11;\n', ...
'dydt(7)=y(2)*y(16)*KonA22-y(7)*KoffA22-y(7)*KoffA11;\n', ...
'dydt(8)=y(3)*y(16)*KonA11-y(8)*KoffA12-y(8)*KoffA11;\n', ...
'dydt(9)=y(3)*y(16)*KonA21-y(9)*KoffA21-y(9)*KoffA12;\n', ...
'dydt(10)=y(4)*y(16)*KonA12-y(10)*KoffA21-y(10)*KoffA12;\n', ...
'dydt(11)=y(4)*y(16)*KonA22-y(11)*KoffA22-y(11)*KoffA21;\n', ...
'dydt(12)=y(5)*y(16)*KonA11-y(12)*KoffA22-y(12)*KoffA11;\n', ...
'dydt(13)=y(5)*y(16)*KonA21-y(13)*KoffA22-y(13)*KoffA21;\n', ...
'dydt(14)=y(2)*Kons11*KonA22+y(5)*Kons11*KonA11-y(14)*KoffA22-y(14)*KoffA11;\n', ...
'dydt(15)=y(3)*Konr11*KonA21+y(4)*Konr11*KonA12-y(15)*KoffA21-y(15)*KoffA12;\n'];


if N==1 &&M==1
  fprintf(diff,diffs11);
elseif N==2 && M==2
  fprintf(diff,diffs22);
end
  fprintf(diff,'end'); %Close the function with end  
  fclose(diff);        %Close the .m function to be usable
end