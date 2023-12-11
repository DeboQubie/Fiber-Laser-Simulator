clc 
close all 
clear all

bit_rate = 2.5*(10^9);

b = [0 1 0 1 1 1 1 0 1 1 1 0];
L = length(b);

t = 0:(.01/bit_rate):(L/bit_rate)-(.01/bit_rate);

j=1;
p=[];

for i=1:length(t)
    if t(i)<=(j/bit_rate)
        p(i) = b(j);
    else
        j = j+1;
        p(i) = b(j);
    end
end

plot(t,p);
xlabel("Time in sec");
ylabel("Normalized optical power");
title("Digital NRZ signal power");
