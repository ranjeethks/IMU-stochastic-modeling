% Ranjeeth KS, University of Calgary
close all;
clc;
b1=9.5238095238095238095238095238095e-5;
b2=20.5238095238095238095238095238095e-5;
b3=27.5238095238095238095238095238095e-5;

 t=-15000:1/200:15000-(1/200);
   r1 = exp(-b1*abs(t));
    t=-15000:1/200:15000-(1/200);
    t1=-15000:1/200:15000-(1/200);
   r2 = (1+b2*abs(t1)).*exp(-b2*abs(t));
    t2=-15000:1/200:15000-(1/200);
   t1=-15000:1/200:15000-(1/200);
    t3=-15000:1/200:15000-(1/200);
   t=-15000:1/200:15000-(1/200);
   r3 = (1+b3*abs(t1)+(1/3)*b3*b3*(t2.*t3)).*exp(-b3*abs(t));
   plot(r1);hold on;plot(r2,'r');hold on;plot(r3,'g');
