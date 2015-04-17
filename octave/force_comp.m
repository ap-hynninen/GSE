clear;
forceEW=load('../forceEW.txt');
%forceGSEr=load('../forceGSEr.txt');
forceGSEk=load('../forceGSEk.txt');
forceLES=load('../forceLES.txt');
forceSPME=load('../forceSPME.txt');
%diffGSEr = forceEW-forceGSEr;
diffGSEk = forceEW-forceGSEk;
diffLES = forceEW-forceLES;
diffSPME = forceEW-forceSPME;

plot(sum(diffLES.^2,2),'r');
hold on;
plot(sum(diffGSEk.^2,2),'k');
%plot(sum(diffGSEr.^2,2));
plot(sum(diffSPME.^2,2),'g');
hold off;