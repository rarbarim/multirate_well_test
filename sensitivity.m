function Pw_sens = sensitivity(q,t2,t3,x)
rw = 0.25; % m
h = 100; % m
k = 10^-12; % m2
Pi = 1.5*10^7 ; % m2
ceff = 2*10^-9; % Pa^-1
S = 0;
visc = 0.02; % Pa.S
dens = 900; % kg/m
por = 0.25;
gamma = exp(0.5722);

sens1 = [rw;h;k;Pi;ceff;S;visc;dens;por];
sens = [sens1 sens1 sens1];
sens(x,:) = [0.8*sens(x,1) sens(x,2) 1.2*sens(x,3)];

S = [ 0 0 0];
Pi= [1.5*10^7 1.5*10^7 1.5*10^7];
if x == 6
    S = [1 2 3];
elseif x == 4
    Pi = [0.8*1.5*10^7 1.5*10^7 1.2*1.5*10^7];
end
 
pw_sens =zeros(length(q),length(t3));
for i= 1 : 3
    A(i) = (sens(7,i))/(4*pi.*sens(3,i).*sens(2,i));
    B(i) = (4.*sens(3,i))/(gamma.*sens(9,i).*sens(7,i).*sens(5,i).*sens(1,i).^2); 
for j = 1 : length(t3)
    if j <= t2(1)
        pw_sens(i,j) =  - A(i)*q(1)*log(S(i)+B(i)*t3(j));
    elseif j <= t2(2) && j > t2(1)
        pw_sens(i,j) = - A(i)*q(1)*log(S(i)+B(i)*t3(j)) -A(i)*(q(2)-q(1))*log(S(i)+B(i)*(t3(j)-t2(1)));
    elseif j <= t2(3) && j > t2(2)
        pw_sens(i,j) = - A(i)*q(1)*log(S(i)+B(i)*t3(j)) -A(i)*(q(2)-q(1))*log(S(i)+B(i)*(t3(j)-t2(1)))-A(i)*(q(3)-q(2))*log(S(i)+B(i)*(t3(j)-t2(2)));
    elseif j <= t2(4) && j > t2(3)
        pw_sens(i,j) = - A(i)*q(1)*log(S(i)+B(i)*t3(j)) -A(i)*(q(2)-q(1))*log(S(i)+B(i)*(t3(j)-t2(1)))-A(i)*(q(3)-q(2))*log(S(i)+B(i)*(t3(j)-t2(2)))-A(i)*(q(4)-q(3))*log(S(i)+B(i)*(t3(j)-t2(3)));    
    elseif j <= t2(5) && j > t2(4)
        pw_sens(i,j) = - A(i)*q(1)*log(S(i)+B(i)*t3(j)) -A(i)*(q(2)-q(1))*log(S(i)+B(i)*(t3(j)-t2(1)))-A(i)*(q(3)-q(2))*log(S(i)+B(i)*(t3(j)-t2(2)))-A(i)*(q(4)-q(3))*log(S(i)+B(i)*(t3(j)-t2(3)))-A(i)*(q(5)-q(4))*log(S(i)+B(i)*(t3(j)-t2(4))); 
    end
end
Pw_sens(i,1:length(t3)+1) = [ Pi(i) (Pi(i)+pw_sens(i,:))];

end
return 