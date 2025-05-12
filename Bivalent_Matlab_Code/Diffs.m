function dydt = Diffs(t,y) 

global Kon Koff Kon2 Koff2 Kon3 Koff3 Kon4 Koff4 L R0 KMTL w N M MTL Kons11 ratioTable currentRow
km=MTL*KMTL;
KonA11=Kon(w);
KonA12=Kon2(w);
KonA21=Kon3;
KonA22=Kon4;
KoffA11=Koff;
KoffA12=Koff2;
KoffA21=Koff3;
KoffA22=Koff4;
Konr11=7.465401e-04; 
Kons00=Kon; 
Kons11=9.711442e-04; 
dydt = zeros(17,1);
dydt(1)=y(2)*KoffA11+y(3)*KoffA12+y(4)*KoffA21+y(5)*KoffA22-y(1)*y(16)*KonA11-y(1)*y(16)*KonA12-y(1)*y(16)*KonA21-y(1)*y(16)*KonA22;
dydt(16)=km*L-km*y(16)+y(2)*KoffA11+y(3)*KoffA12+y(4)*KoffA21+y(5)*KoffA22+y(6)*KoffA12+y(6)*KoffA11+y(7)*KoffA22+y(7)*KoffA11+y(8)*KoffA12+y(8)*KoffA11+y(9)*KoffA21+y(9)*KoffA12+y(10)*KoffA21+y(10)*KoffA12+y(11)*KoffA22+y(11)*KoffA21+y(12)*KoffA22+y(12)*KoffA11+y(13)*KoffA22+y(13)*KoffA21-y(1)*y(16)*KonA11-y(1)*y(16)*KonA12-y(1)*y(16)*KonA21-y(1)*y(16)*KonA22-y(2)*y(16)*KonA12-y(2)*y(16)*KonA22-y(3)*y(16)*KonA11-y(3)*y(16)*KonA21-y(4)*y(16)*KonA12-y(4)*y(16)*KonA22-y(5)*y(16)*KonA11-y(5)*y(16)*KonA21;
dydt(17)=y(16)*0.01;
dydt(2)=y(1)*y(16)*KonA11-y(16)*y(2)*KonA12-y(16)*y(2)*KonA22-y(2)*Kons11*KonA22+y(6)*KoffA12+y(7)*KoffA22+y(8)*KoffA12+y(12)*KoffA22+y(14)*KoffA22-y(2)*KoffA11;
dydt(3)=y(1)*y(16)*KonA12-y(16)*y(3)*KonA11-y(16)*y(3)*KonA21-y(3)*Konr11*KonA21+y(6)*KoffA11+y(8)*KoffA11+y(9)*KoffA21+y(10)*KoffA21+y(15)*KoffA21-y(3)*KoffA12;
dydt(4)=y(1)*y(16)*KonA21-y(16)*y(4)*KonA12-y(16)*y(4)*KonA22-y(4)*Konr11*KonA12+y(9)*KoffA12+y(10)*KoffA12+y(11)*KoffA22+y(13)*KoffA22+y(15)*KoffA12-y(4)*KoffA21;
dydt(5)=y(1)*y(16)*KonA22-y(16)*y(5)*KonA11-y(16)*y(5)*KonA21-y(5)*Kons11*KonA11+y(7)*KoffA11+y(11)*KoffA21+y(12)*KoffA11+y(13)*KoffA21+y(14)*KoffA11-y(5)*KoffA22;
dydt(6)=y(2)*y(16)*KonA12-y(6)*KoffA12-y(6)*KoffA11;
dydt(7)=y(2)*y(16)*KonA22-y(7)*KoffA22-y(7)*KoffA11;
dydt(8)=y(3)*y(16)*KonA11-y(8)*KoffA12-y(8)*KoffA11;
dydt(9)=y(3)*y(16)*KonA21-y(9)*KoffA21-y(9)*KoffA12;
dydt(10)=y(4)*y(16)*KonA12-y(10)*KoffA21-y(10)*KoffA12;
dydt(11)=y(4)*y(16)*KonA22-y(11)*KoffA22-y(11)*KoffA21;
dydt(12)=y(5)*y(16)*KonA11-y(12)*KoffA22-y(12)*KoffA11;
dydt(13)=y(5)*y(16)*KonA21-y(13)*KoffA22-y(13)*KoffA21;
dydt(14)=y(2)*Kons11*KonA22+y(5)*Kons11*KonA11-y(14)*KoffA22-y(14)*KoffA11;
dydt(15)=y(3)*Konr11*KonA21+y(4)*Konr11*KonA12-y(15)*KoffA21-y(15)*KoffA12;
end