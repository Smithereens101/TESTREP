clear;
clc;
close all;

%Hertzian Indentation Stress Field EquationsDavid Louapre,* and Kristin Breder
%Int. J. Appl. Ceram. Technol., 12 [5] 1071–1079 (2015)DOI:10.1111/ijac.12317
%https://ceramics.onlinelibrary.wiley.com/doi/epdf/10.1111/ijac.12317?saml_referrer
%nu - poisson's ratio of the halspace,
%nu_p - nu prime - poisson's ratio of the sphere
%a - half of the contact area radiuss
%E - Young's modulus of the halfspace
%E_p - Young's modulus of the sphere
%ro - radius of the spherical indenter
%P - normal force

%syms r z nu nu_p E_p ro P

N = 500;
p0 = 3.4795;
s_0 = 3/2*p0;
a = 1;
nu = 0.25;
rlist = linspace(1e-4,1.5,N);
zlist = linspace(1e-4,1.5,N);


i = 0;
j = 0;
stress = zeros(3,3);
s1 = zeros(1,N);
s2 = zeros(1,N);
s3 = zeros(1,N);
s13= zeros(1,N);

for r = rlist
    i = i+1;
       stress =  alternative_hertz_stress(r, 1e-3, a, nu, p0)/s_0;
       s1(i) = 0.5*(stress(1,1)+stress(3,3))+sqrt((0.5*(stress(1,1)-stress(3,3)))^2+stress(3,1)^2);
       s2(i) = stress(2,2); 
       s3(i) = 0.5*(stress(1,1)+stress(3,3))-sqrt((0.5*(stress(1,1)-stress(3,3)))^2+stress(3,1)^2);
       s13(i)= 0.5*(s1(i)-s3(i));
end

%% Figure 2

figure
hold on

plot(rlist/a,s1,'-',LineWidth=0.8)
plot(rlist/a,s2,'--',LineWidth=0.8)
plot(rlist/a,s3,'-.',LineWidth=0.8)
xlabel('r/a',FontSize=14)
ylabel('$\sigma$/$\sigma _0$',interpreter = 'latex',FontSize=14)
legend('\sigma _1','\sigma _2','\sigma _3',Location='southeast')
set(gca,"FontWeight","bold")




