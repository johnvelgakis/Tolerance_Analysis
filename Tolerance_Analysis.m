%Uniform Distribution:

% define nominal and tolerance parameters
Rnom=1e3;
dR=0.05;
%loop to generate 1000 resistors
for ind=1:1000
% determine individual component value R(ind)=Rnom*(1+2*(rand-0.5)*dR);
end
%plot a histogram
hist(R);

%Gaussian (normal) Probability Distribution:

% define nominal and tolerance parameters
Rnom=1e3;
dR=0.05;
%loop to generate 1000 resistors
for ind=1:1000
% determine individual component value R(ind)=Rnom*(1+dR/3*randn);
end
%plot a histogram
hist(R);

%Potential Divider:

   % A2 The Matlab program code for a uniform probability distribution Monte Carlo simulation
% to determine the gain calculation of a potential divider circuit
% define nominal and tolerance parameters
% R1=3k? and R2=1k?.
Rnom1=3e3; Rnom2=1e3; dR=0.01;
%loop to generate 1000 resistor pairs R1, R2
for ind=1:1000
% determine individual component value G(ind)=Rnom2*(1+2*(rand-0.5)*dR)/( Rnom1*(1+2*(rand-0.5)*dR)+ Rnom2*(1+2*(rand-0.5)*dR));
end
%plot a histogram
hist(G);
Gnom=Rnom2/(Rnom1+ Rnom2); Gmean = mean(G);
% standard deviation
s = std(G);
% tolerance TG
TG=s/Gnom;
% A2a ? Supplement Matlab program code of routine A2
% to determine the yield of a potential divider circuit % define tolerance parameters
dGnom=0.02;
%loop to determine the yield (how many circuits meet the specification)
n=0;
for ind=1:1000
% determine how many meet the specification dG= abs((G(ind)-Gnom)/Gnom);
if (dG <= dGnom) n=n+1;
end
end
yield=n/1000;

%Low pass filter:

% A3 ? The Matlab program code for a Gaussian distribution Monte Carlo simulation
% to determine the corner frequency fc=1kHz ?5% where fc=1/(2? RC) of a potential divider
% define nominal and tolerance parameters
% component values R=1.6k? and C=100nF ?10%. Rnom=1.6e3;
Cnom=100e-9;
dR=0.01;
dC=0.1;
%loop to generate 1000 resistor pairs R1, R2
ntot=1000;
for ind=1:ntot
% determine individual component value
f(ind)=1/( 2*pi*(Rnom*(1+dR/3*randn)*Cnom*(1+dC/3*randn)));
end
%plot a histogram
hist(f);
fnom=1/(2*pi*Rnom*Cnom); fmean = mean(f);
% standard deviation
s = std(f);
% tolerance Tfc
Tfc=s/fnom;
% A3a ? Supplement Matlab program code of routine A3
% to determine the yield of a potential divider circuit % define tolerance parameters
dfnom=0.05;
%loop to determine the yield (how many circuits meet the specification)
n=0;
for ind=1:ntot
% determine how many meet the specification df= abs((f(ind)-fnom)/fnom);
if (df <= dfnom) n=n+1;
end
end
yield=n/ntot;

%2-bit R-2R ladder DAC:

 % A4 – The Matlab program code for a Gaussian distribution Monte Carlo simulation
% to calculate the output voltage of a DAC for various input codes with uniform distribution % of resistors of various tolerances
% define nominal and tolerance parameters
R1=20e3;
R2=10e3;
R3=20e3;
R4=20e3;
Vref=5;
LSB=1.25;
LSBhalf=LSB/2;
% (a0, a1) is the digital input
dR=0.0;
%Vout2(1)=0;
for ind2=1:1500 dR=dR+0.001;
for ind1=1:3 a0=mod(ind1,2); a1=round((ind1-a0)/2);
if a1 == 2; a1=0; end
d1=Vref*a1;
d0=Vref*a0;
Voutnom=(d1*(1/R1 + R2/(R1*R3) + R2/(R1*R4) ) + d0/R3 )/(1/R1 +1/R3 + 1/R4 + R2/(R1*R3) + R2/(R1*R4) );
for ind=1:1000
% determine individual component value
Vout(ind)=(d1*(1/(R1*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand- 0.5)*dR)*R3*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R4*(1+2*(rand- 0.5)*dR)))+d0/(R3*(1+2*(rand-0.5)*dR)))/(1/(R1*(1+2*(rand-0.5)*dR))+1/(R3*(1+2*(rand- 0.5)*dR))+1/(R4*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R3*(1+2*(rand- 0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R4*(1+2*(rand-0.5)*dR)));
end
Voutmean = mean(Vout);
% standard deviation
s = std(Vout);
% tolerance TVout
TVout=s/Voutnom; Vout2(ind1)=Voutmean;
end
%plot a histogram %hist(Vout2);
A(1)=0; for i=1:3
A(i+1)=Vout2(i);
end

for j=1:3 DNL(j)=abs((A(j+1)-A(j)-LSB)/LSB);
end
maxxDNL=max(DNL); DNLmax(ind2)=maxxDNL; DNLmin=min(DNL);
if maxxDNL > LSBhalf
return % Stop the ind2 for loop else
savedR=dR;
end end
%plot a histogram %hist(DNL);
% A5 – The Matlab program code for a Gaussian distribution Monte Carlo simulation % to calculate the Gain error at full-scale,
% for various input codes with a uniform distribution of resistors of various tolerances % define nominal and tolerance parameters
R1=20e3; R2=10e3; R3=20e3; R4=20e3; dR=0.0; Vref=5;
% start up the counter for the gain yield loop
n=0;
%for ind1=1:1
% (a0, a1) is the digital input % a0=mod(ind1,2);
% a1=round((ind1-a0)/2); % if a1 == 2;
% a1=0;
% end
for ind2=1:1500
dR=dR+0.001; a0=1;
a1=1;
d1=Vref*a1;
d0=Vref*a0;
Voutnom=(d1*(1/R1 + R2/(R1*R3) + R2/(R1*R4) ) + d0/R3 )/(1/R1 +1/R3 + 1/R4 + R2/(R1*R3) + R2/(R1*R4) ); VoutFS=3.75;
LSB=1.25;
LSBtenth=0.1*LSB;
Gainnom=(Voutnom-VoutFS)/LSB;
for ind=1:1000
% determine individual component value
Voutind=(d1*(1/(R1*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand- 0.5)*dR)*R3*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R4*(1+2*(rand- 0.5)*dR)))+d0/(R3*(1+2*(rand-0.5)*dR)))/(1/(R1*(1+2*(rand-0.5)*dR))+1/(R3*(1+2*(rand-0.5)*dR))+1/(R4*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R3*(1+2*(rand- 0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R4*(1+2*(rand-0.5)*dR))); Gain(ind)=(Voutind-VoutFS)/LSB;
end
Gainmean = mean(Gain);
% standard deviation
s = std(Gain);
% tolerance TGain
TGain=s/Gainnom;
% Vout2(ind1)=Voutmean; % end
%plot a histogram
hist(Gain);
%find the maximum and minimum values of Gain error for each resistor tolerance value
maxxGain=max(Gain); maxGain(ind2)=maxxGain; minGain=min(Gain);
if maxxGain > LSBtenth
return % Stop the ind2 for loop
else
savedR=dR; saveindex2=ind2; end
end

% A6 – The Matlab program code for a Gaussian distribution Monte Carlo simulation
% to calculate the DNL for various input codes with uniform distribution of resistors and % various tolerances.
% define nominal and tolerance parameters
R1=20e3; R2=10e3; R3=20e3; R4=20e3; dR=0.01; Vref=5;
% start up the counter for the gain yield loop
n=0;
%for ind1=1:1
% (a0, a1) is the digital input a0=1;
a1=1;
% a0=mod(ind1,2);
% a1=round((ind1-a0)/2);
% if a1 == 2;
% a1=0;
% end
d1=Vref*a1;
d0=Vref*a0;
Voutnom=(d1*(1/R1 + R2/(R1*R3) + R2/(R1*R4) ) + d0/R3 )/(1/R1 +1/R3 + 1/R4 + R2/(R1*R3) + R2/(R1*R4) ); VoutFS=3.75;
LSB=1.25;
for ind=1:1000

% determine individual component value
Voutind=(d1*(1/(R1*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand- 0.5)*dR)*R3*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R4*(1+2*(rand- 0.5)*dR)))+d0/(R3*(1+2*(rand-0.5)*dR)))/(1/(R1*(1+2*(rand-0.5)*dR))+1/(R3*(1+2*(rand- 0.5)*dR))+1/(R4*(1+2*(rand-0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R3*(1+2*(rand- 0.5)*dR))+(R2*(1+2*(rand-0.5)*dR))/(R1*(1+2*(rand-0.5)*dR)*R4*(1+2*(rand-0.5)*dR))); Vout(ind)=(Voutind-VoutFS)/LSB;
end
Voutmean = mean(Vout);
% standard deviation
s = std(Vout);
% tolerance TVout
TVout=s/Voutnom;
% Vout2(ind1)=Voutmean; % end
%plot a histogram
hist(Vout);
%find the maximum and minimum values of Gain error for each resistor tolerance value maxVout=max(Vout);
minVout=min(Vout);

