function EffectiveConcentration = EffC_Calculator(KonLocal)

resolution=1;
global w LinkerNLig UnitDiameterLig LinkerNRec UnitDiameterRec D type LinkerTypeLig LinkerTypeRec %****************

LinkerTypeLig='flexible'; LinkerTypeRec='flexible'; notalpha=0;
LinkerN1=LinkerNLig; UnitDiameter1=UnitDiameterLig;
LinkerN2=LinkerNRec; UnitDiameter2=UnitDiameterRec;
linkerlpLig = 30; %30;%702.21; %NEW %200 %234.07 %185.418 %702.21
linkerlpRec = 30; %576.8745; %NEW %1000 %192.2915 %180.9115 %576.8745


    
%% Parameterizing the matrix and linker sizes
modifier=1;
if sum(KonLocal=='X')==0
        if KonLocal(4)=='s' %input is Kons11
            reverse=0;
            
        elseif KonLocal(4)=='r' %input is Konr11
            reverse=1;
            %modifier=0.03;%!!!!!!!!!!!!
        end
    UnitN1=str2double(KonLocal(5))+1;
    UnitN2=str2double(KonLocal(6))+1;
elseif sum(KonLocal=='X')==2
    UnitN1=min(str2double(KonLocal(6)),str2double(KonLocal(7)))+1;
    UnitN2=min(str2double(KonLocal(9)),str2double(KonLocal(10)))+1;
        if KonLocal(4)=='s'
            
            modifier=5/((1+(min(LinkerNLig, LinkerNRec)+1)/5)*(1+abs(LinkerNLig - LinkerNRec)/10));
            reverse=0;
        elseif KonLocal(4)=='r'
            modifier=1/((1+(min(LinkerNLig, LinkerNRec)+1)/10)*(1+abs(LinkerNLig - LinkerNRec)/10));
            reverse=1;
        end
elseif find(KonLocal=='X')==5
    UnitN1=min(str2double(KonLocal(6)),str2double(KonLocal(7)))+1;
    UnitN2=str2double(KonLocal(8))+1;
    modifier=1;
    reverse=0;
elseif find(KonLocal=='X')==6
    UnitN1=str2double(KonLocal(5))+1;
    UnitN2=min(str2double(KonLocal(7)),str2double(KonLocal(8)))+1;
    modifier=1;
    reverse=0;
end

%% The receptor PDF calculation
switch LinkerTypeRec
    case 'flexible'
    PDF2=PDFMatrix(LinkerN2, UnitN2, UnitDiameter2, type, [linkerlpRec,1], 'receptor', resolution); %[200,1]
    case 'alpha'
    PDF2=PDFMatrix(LinkerN2, UnitN2, UnitDiameter2+notalpha*3, type, [linkerlpRec,0.55], 'receptor', resolution); %[1000,0.55]    
end

%% THe ligand PDF calculation
switch LinkerTypeLig
    case 'flexible'
        PDF1=PDFMatrix(LinkerN1, UnitN1, UnitDiameter1, type, [linkerlpLig,1], 'ligand', resolution); %[200,1]
    case 'alpha'
        PDF1=PDFMatrix(LinkerN1-notalpha, UnitN1, UnitDiameter1+notalpha*3, type, [linkerlpLig,0.55], 'ligand', resolution); %[1000,0.55]
end

if length(PDF1)>length(PDF2)
    PDFsmaller=zeros(size(PDF1));
    border=(length(PDF1)-length(PDF2))/2;
    PDFsmaller(border+1:length(PDF1)-border,border+1:length(PDF1)-border,border+1:length(PDF1)-border)=PDF2;
    PDF2=PDFsmaller;
elseif length(PDF1)<length(PDF2)
    PDFsmaller=zeros(size(PDF2));
    border=(length(PDF2)-length(PDF1))/2;
    PDFsmaller(border+1:length(PDF2)-border,border+1:length(PDF2)-border,border+1:length(PDF2)-border)=PDF1;
    PDF1=PDFsmaller;
end
%% Calculation of the [Leff]

if reverse ==0
    shift = ceil((UnitDiameter2 + UnitDiameter1)/2);
elseif reverse ==1
    shift= UnitDiameter2 + UnitDiameter1;
end
%shift = 0 %!!!!!!!!

PDF1 = [PDF1;zeros(shift, size(PDF1,2),size(PDF1,3))];
size(PDF1);
PDF2 = [zeros(shift, size(PDF2,2),size(PDF2,3));PDF2];

clear D DistanceMatrix
EffectiveConcentration=sum(sum(sum(PDF1.*PDF2)))/((6*(10^-4))/resolution^3)*modifier;

end

function PDF_Matrix=PDFMatrix(LinkerN, UnitN, UnitDiameter, type, linkerlp, recorlig, resolution)
%The function representing a section of the Matrix similarly as the whole
%Input: LinkerN: linker unit number in piece; UnitN: The number of units of Receptor and Ligand
%       UnitDiameter: The diameter of the Receptor and Ligand
%       type: the general type of reaction, resolution: by default it is 1A

%% Initial steps

global DistanceMatrix D

if nargin==4; resolution=1; end         %The argument setup
if type==1; LinkerU=3; end              %the type setup
if UnitN==1; PDF_Matrix=1;return; end

%% Calculating the main variables, matrices

%Linker and size calculation-----------------------------------------------
LinkerL=ceil(LinkerN*LinkerU);                %The linker length calculation
MatrixWidth=ceil(max(LinkerL,UnitDiameter)); %The width of the matrix
Size=UnitN*UnitDiameter+(UnitN-1)*LinkerL+2*MatrixWidth; %The size of the matrix

%Distance matrix calculation-----------------------------------------------
    n=ceil(Size*resolution);            %the matrix radius determined by the resolution and the Ligand length
    D=zeros((2*n+1),(2*n+1),(2*n+1));   %The zero matrix with the middle 1
    D(n+1,n+1,n+1)=1;                   %The 1 in the middle starts the count of the distance
    D=bwdist(D);                        %Creates the distance matrix
    D=double(D);                        %Converts the matrix to numerical matrix
    D=D/resolution;                     %It should start as around zero
    

DistanceMatrix=D(1:end,n+1-MatrixWidth:n+1+MatrixWidth,n+1-MatrixWidth:n+1+MatrixWidth); %!!!!!!!!!!!!!!!!!!!! Need to reduce to incresase speed
%Creates a Matrix of a relevant smaller proportion of the whole matrix to reduce the computation

%PDF calculation-----------------------------------------------------------
PDF=EffectiveCC(DistanceMatrix,LinkerL, linkerlp);

%Unit01Matrix calculation for addUnit()------------------------------------
switch recorlig
    case 'receptor'
        Unit01Matrix=D>=UnitDiameter&D<(UnitDiameter+1);
        Unit01Matrix=Unit01Matrix(n+1-UnitDiameter:n+1+UnitDiameter,n+1-UnitDiameter:n+1+UnitDiameter,n+1-UnitDiameter:n+1+UnitDiameter);
        Unit01Matrix=Unit01Matrix.*ones(size(Unit01Matrix));
        x=2;
%E(ceil(length(E)*0.8):end,:,:)=0;  %Would be some sort of angle correction !!!!!!!!!!!
    case 'ligand'
        Unit01Matrix=D(find(D(:,n+1,n+1)<=UnitDiameter),find(D(:,n+1,n+1)<=UnitDiameter),find(D(:,n+1,n+1)<=UnitDiameter)); %#ok<FNDSB>
        Unit01Matrix=EffectiveCC(Unit01Matrix,UnitDiameter, [9,1]);
end
%LinkerPFDMatrix calcualtion for addLinker()-------------------------------
LinkerPFDMatrix=D(find(D(:,n+1,n+1)<=LinkerL),find(D(:,n+1,n+1)<=LinkerL),find(D(:,n+1,n+1)<=LinkerL)); %#ok<FNDSB>
LinkerPFDMatrix=EffectiveCC(LinkerPFDMatrix,LinkerL, linkerlp);
x=2;

%% Initiate an at least 2Unit containing matrix

%double addUnit to add two units to the end of the linker------------------
%PDF=addUnit(PDF, UnitDiameter, Unit01Matrix);
PDF=addUnit(PDF, UnitDiameter, Unit01Matrix);

%% Repairs the Linker and Unit adding proccess till needed
for i=1:(UnitN-2)
PDF=addLinker(PDF, LinkerL ,LinkerPFDMatrix);
PDF=addUnit(PDF, UnitDiameter, Unit01Matrix);
end

%% The final calculation
PDF_Matrix=PDF2Dto3D(PDF((n+1):end,MatrixWidth+1,MatrixWidth+1), D);
PDF_Matrix=PDF_Matrix/sum(PDF_Matrix(:));


end

function [PDF_added] = addUnit(PDF, UnitDiameter, Unit01Matrix)
%PDF in any form, X axis must be the longest, E is a 0-1 matrix
%E(ceil(length(E)*0.8):end,:,:)=0;  %Would be some sort of angle correction !!!!!!!!!!!
Unit01Matrix(ceil(length(Unit01Matrix)*0.8):end,:,:)=0;

    PDFSize=size(PDF);                                                          %Get the size of the PDF
    PDFmiddleX=ceil(PDFSize(1)/2);                                              %Get the center of the PDF 
    PDFmiddleYZ=ceil(PDFSize(2)/2);                                             %get the center of the ZY axes
    PDF_added=zeros(1,PDFSize(1)-PDFmiddleX-UnitDiameter);                      %The output 2D vector

        PositionYZ=((PDFmiddleYZ)-UnitDiameter):((PDFmiddleYZ)+UnitDiameter);   %After the ZY coordinates
    for i= 0:(length(PDF_added)-1)                                              %Goes trough the PDF added matrix
        newUnit01Matrix = Unit01Matrix;
        newUnit01Matrix(ceil(length(newUnit01Matrix)*0.8):end,:,:)=0;
        PositionX=((PDFmiddleX+i)-UnitDiameter):((PDFmiddleX+i)+UnitDiameter);  %The position of the X axis: the matrix with the center of the PDF added point
        PDF_added(i+1)=sum(sum(sum(newUnit01Matrix.*PDF(PositionX,PositionYZ,PositionYZ))));  %Add the value~ to the correct position
        % The value~ is the E 0-1(in R dustance)matrix multiplied by the local matrix of interrest in the PDF
    end

    PDF_added=PDF2Dto3D(PDF_added);
    
end

function PDF_added = addLinker(PDF, LinkerL ,LinkerPFDMatrix)
    %
    
    PDFSize=size(PDF);                                                          %Get the size of the PDF
    PDFmiddleX=ceil(PDFSize(1)/2);                                              %Get the center of the PDF 
    PDFmiddleYZ=ceil(PDFSize(2)/2);                                             %get the center of the ZY axes
    PDF_added=zeros(1,PDFSize(1)-PDFmiddleX-LinkerL);                           %The output 2D vector

        PositionYZ=((PDFmiddleYZ)-LinkerL):((PDFmiddleYZ)+LinkerL);
    for i= 0:(length(PDF_added)-1)
        PositionX=((PDFmiddleX+i)-LinkerL):((PDFmiddleX+i)+LinkerL);
        PDF_added(i+1)=sum(sum(sum(LinkerPFDMatrix.*PDF(PositionX,PositionYZ,PositionYZ))));
    end
    

    PDF_added=PDF2Dto3D(PDF_added);

end

function PDF = PDF2Dto3D(PDF2D, D)
PDF2D=PDF2D/sum(PDF2D);
global DistanceMatrix;
if nargin==2
    localDistanceMatrix=D;
else
    localDistanceMatrix=DistanceMatrix;
end

    if size(PDF2D,1)>1
       PDF2D=PDF2D'; 
    end
    Location_relevant=find(PDF2D>1e-15);                                    %Distance from 0 in A
    PDF2D_relevant=PDF2D(Location_relevant);                                %PDF of the relevant location (PDF)
    
    Location_from_0=Location_relevant-min(Location_relevant);               %Distance from the start of the relevant loc (A)
    %PolyCoeffs = polyfit(Location_from_0,PDF2D_relevant,5);                 %Fitting (from shifted to 0)
    curve = fit(Location_from_0',PDF2D_relevant','smoothingspline');
    

    PDF=zeros(size(localDistanceMatrix));
    %figure 
    %plot(PDF2D_relevant)
    PDF(min(Location_relevant)<=localDistanceMatrix & max(Location_relevant)>=localDistanceMatrix)=...      %The affected part of the PDF
        curve((localDistanceMatrix(min(Location_relevant)<=localDistanceMatrix & max(Location_relevant)>=localDistanceMatrix)-min(Location_relevant)));
        %polyval(PolyCoeffs,(localDistanceMatrix(min(Location_relevant)<=localDistanceMatrix & max(Location_relevant)>=localDistanceMatrix)-min(Location_relevant)));
    
    %figure 
    %plot(PDF(:,ceil(size(localDistanceMatrix,2)/2),ceil(size(localDistanceMatrix,2)/2)));
        %The new PDF gets the fitted values
    PDF(PDF<0)=0; %Sometimes yields negative values close to 0 and the values must be positive
    
end

function [p] = EffectiveCC(x, L, parameters)
%Eff values representing volumes

%lp=9;               %The lp parameter in A units
lp=parameters(1);
L=L*parameters(2);

eff = @(r) (((1./(1-r.^2/L^2).^(4.5)).*exp(-(9*L/(8*lp))*(1./(1-r.^2/L^2)))));
%the eff is the given distribution for the peptide PDF

border=(x>-L&x<L);  %Logically FALSE outside the range and TRUE inside
x=x.*border;        %A correction factor: if param>L then 0
p=eff(x);           %The calculation of the vector to the PDF values
p=p.*border;        %The output is cropped also

end

