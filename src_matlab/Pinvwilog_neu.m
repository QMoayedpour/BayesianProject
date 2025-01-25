function fd = pinvwilog_neu(Qcol,detQ1log,nu,Scol,detSlog)

% dichte einer invertiert wishart-verteilung 
%input: Qcol ... argument (matrix INVERTIERT!, wird als spalte übergeben)
%       mehrere Argumente: anzhal der Spalten entspricht ANzahl der Argumente
%       detQ1log ... log(det(Q1))
%       nu ... shape parameter skalar:identisch für alle argumente
%                              zeilenvektor: varriert mit dem argument
%       Scol  ...  scale matrix, wird als spalte übergeben
%       detSlog ...  log(det(Scol)) (zeilenvektor)
% ergegbis: funktionswert skalar oder Zeile

D = size(Qcol,1);
r = -.5 + (.25 + 2*D)^.5; %anzahl der Zeilen von Q wird aus der dimension von Qcol rekontruiert
index=0.5*[[1:r].*[2:r+1]]';
if (size(nu,2)==1); nu=nu(1,size(Scol,2));end
%trQS=trace(S*Q1);
% berechnung aus den Spaltenvektoren
trQS=2*sum(Qcol.*Scol,1)-sum(Qcol(index,:).*Scol(index,:),1);
fd=nu.*detSlog+(nu+(r+1)/2)*detQ1log-trQS-r*(r-1)/4*log(pi);
index2=[1:r]';
fd=fd-sum(gammaln(nu(ones(r,1),:)+0.5-0.5*index2(:,ones(1,size(nu,2)))),1);