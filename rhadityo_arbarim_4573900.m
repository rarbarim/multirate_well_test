% Multirate Well Testing
% R B Arbarim : 4573900

clear all

disp('Multi Rate Test Simulator by RB Arbarim - 4573900')
disp(' ')
for i = 1 : 5
       fprintf('Q%s\n',num2str(i))
       q(i) = str2double(input('enter rate in m3/s (recommended value not more than 0.02 m3/s) = ','s'));
       t(i) = str2double(input('how long does it last in days = ','s'));
end

% q = [0.0200 0 0.0040 0.0090 0.0100];
% t = [3     6     2     4     2];
conv = 86400; % conversion from day to second
t1 = conv.*t; % s
t2 = cumsum(t1);
t3 = 1:sum(t1);
t4 = [ 0 t3./conv];

rw = 0.25; % m
h = 100; % m
k = 10^-12; % m2
Pi = 1.5*10^7; % m2
ceff = 2*10^-9; % Pa^-1
S = 0;
visc = 0.02; % Pa.S
dens = 900; % kg/m
por = 0.25;
Bo = 1; % Rm3/STm3
gamma = exp(0.5722);

A = (visc*Bo)/(4*pi*k*h);
B = (4*k)/(gamma*por*visc*ceff*rw^2);

a = cell(length(q),1);
figure(1)
for i = 1 : length(q)
    hold on
        if i == 1
            Q(i,1:sum(t1)) = q(i);
            plot((1:t2(i))./conv,Q(i,i:t2(i)),'linewidth',3);
        else
            Q(i,(t2(i-1)):sum(t1)) = q(i);
            plot((t2(i-1):t2(i))./conv,Q(i,t2(i-1):t2(i)),'linewidth',3 );
        end
    a{i} = ['Q' num2str(i) ' = ' num2str(q(i)) ' m3/s'];
end

legend(a)
xlabel('time (days)')
ylabel('Rate (m3/s)')
hold off

for i = 1 : length(t3)
    if i <= t2(1)
        pw(i) =  - A*q(1)*log(S+B*t3(i));
    elseif i <= t2(2) && i > t2(1)
        pw(i) = - A*q(1)*log(S+B*t3(i)) -A*(q(2)-q(1))*log(S+B*(t3(i)-t2(1)));
    elseif i <= t2(3) && i > t2(2)
        pw(i) = - A*q(1)*log(S+B*t3(i)) -A*(q(2)-q(1))*log(S+B*(t3(i)-t2(1)))-A*(q(3)-q(2))*log(S+B*(t3(i)-t2(2)));
    elseif i <= t2(4) && i > t2(3)
        pw(i) = - A*q(1)*log(S+B*t3(i)) -A*(q(2)-q(1))*log(S+B*(t3(i)-t2(1)))-A*(q(3)-q(2))*log(S+B*(t3(i)-t2(2)))-A*(q(4)-q(3))*log(S+B*(t3(i)-t2(3)));    
    elseif i <= t2(5) && i > t2(4)
        pw(i) = - A*q(1)*log(S+B*t3(i)) -A*(q(2)-q(1))*log(S+B*(t3(i)-t2(1)))-A*(q(3)-q(2))*log(S+B*(t3(i)-t2(2)))-A*(q(4)-q(3))*log(S+B*(t3(i)-t2(3)))-A*(q(5)-q(4))*log(S+B*(t3(i)-t2(4))); 
    end
end

figure(2)
Pw = [ Pi (Pi+pw)];
plot(t4,Pw)
xlabel('time (days)')
ylabel('Wellbore Pressure (Pa)')

disp(' ')
disp('do sensitivity plus and minus 20% of its initial value');
disp('press 1 to do sensitivity in wellbore radius (rw)')
disp('press 2 to do sensitivity in reservoir thickness (h)')
disp('press 3 to do sensitivity in permeability (k)')
disp('press 4 to do sensitivity in initial pressure (Pi)')
disp('press 5 to do sensitivity in compressibility (Ceff)')
disp('press 6 to do sensitivity in skin (S)')
disp('press 7 to do sensitivity in viscosity (mu)')
disp('press 8 to do sensitivity in density (rho)')
disp('press 9 to do sensitivity in porosity')
x = str2double(input('which parameter (1 to 9) do you want to do sensitivity? = ','s'));
disp('see figure 3 for the sensitivity result')

Pw_sens = sensitivity(q,t2,t3,x);
b = cell(9,3);
b(1,1:3) = {'0.8 rw' , 'rw' , '1.2 rw'};
b(2,1:3) = {'0.8 h' , 'h' , '1.2 h'};
b(3,1:3) = {'0.8 k', 'k' , '1.2 k'};
b(4,1:3) = {'0.8 Pi' , 'Pi' , '1.2 Pi'};
b(5,1:3) = {'0.8 ceff' , 'ceff',  '1.2 ceff'};
b(6,1:3) = {'0.8 S',  'S',  '1.2 S'};
b(7,1:3) = {'0.8 mu',  'mu',  '1.2 mu'};
b(8,1:3) = {'0.8 rho' , 'rho',  '1.2 rho'};
b(9,1:3) = {'0.8 por' , 'por' , '1.2 por'};

figure(3)
hold on
for i = 1 : 3
    plot(t4,Pw_sens(i,:))
end
hold off

xlabel('time (days)')
ylabel('Wellbore Pressure (Pa)')
legend(b(x,:))