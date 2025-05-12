function [M, N, rerun, type, plottype, MWreceptor, MWligand, LinkerNLig, UnitDiameterLig, LinkerNRec,...
    UnitDiameterRec, Kon, Koff, Kon2, Koff2, Kon3, Koff3, Kon4, Koff4, OnTime, OffTime, KMTL, L0, baseRU] = parameters
    M = 2;
    N = 2;
    rerun = 1;
    type = 1;
    plottype = 1; %plottype =1 will calculate 
    MWreceptor = 865.94 + 807.9;%736.82 + 833.85;%865.94 + 807.9; %average of both pTyr binding motifs
    MWligand = 11587.07+11764.29;%11587.07+11764.29;%10352.9 + 10704.32;10594.2+10736.19 %average diameter of SH2 domains in PTPN11
    LinkerNLig = 5.75; 
    UnitDiameterLig = 36;%35
    LinkerNRec = 31;%10; %31
    UnitDiameterRec = 19;
    Kon = 50000;%251889.169; 
    Koff = 1; %1
    Kon2 = 1;
    Koff2 = 1; %1
    Kon3 = 1;%1907013.997;
    Koff3 = 1; %1
    Kon4 = 5000000;%426712.182;
    Koff4 = 1; %1
    OnTime = 30; %500;
    OffTime = 30; %500
    KMTL = 1;
    L0 = .25e-6; %4e-6 %1e-6 %initial Ligand concentration. This value should be higher than the expected affinity in a sweep
    baseRU = 25; %25
end