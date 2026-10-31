function [Smodel,I1model,I2model,I3model,I4model, Rmodel, Dmodel] = si4rdsrk4(time,s0,i10,i20,i30,i40,r0,d0,Ptotal,p1,p2,p3,p4,beta1,beta2,beta3,beta4,r1,r2,r3,r4,d1,d2,d3,d4,gamma,p_vax)
Szero = s0;
I1zero = i10;
I2zero = i20;
I3zero = i30;
I4zero = i40;
Rzero = r0;
Dzero = d0;
T = length(time);
N = T;
M = 500;
sigma1 = 0.05;
% mu = 0.000021;
mu = 0;
dt = T/N; R = 1; Dt = R*dt; L = N/R;
global paramestimate
%======================================================
Pastparamestimate = paramestimate;
past_p1      = Pastparamestimate(1,:);
past_p2      = Pastparamestimate(2,:);
past_p3      = Pastparamestimate(3,:);
past_p4      = Pastparamestimate(4,:);
past_beta1   = Pastparamestimate(5,:);
past_beta2   = Pastparamestimate(6,:);
past_beta3   = Pastparamestimate(7,:);
past_beta4   = Pastparamestimate(8,:);
past_r1      = Pastparamestimate(9,:);
past_r2      = Pastparamestimate(10,:);
past_r3      = Pastparamestimate(11,:);
past_r4      = Pastparamestimate(12,:);
past_delta1  = Pastparamestimate(13,:);
past_delta2  = Pastparamestimate(14,:);
past_delta3  = Pastparamestimate(15,:);
past_delta4  = Pastparamestimate(16,:);
past_gamma   = Pastparamestimate(17,:);
past_pvax   = Pastparamestimate(18,:);

past_beta1(past_beta1 == 0) = beta1;
past_beta2(past_beta2 == 0) = beta2;
past_beta3(past_beta3 == 0) = beta3;
past_beta4(past_beta4 == 0) = beta4;
past_r1(past_r1 == 0) = r1;
past_r2(past_r2 == 0) = r2;
past_r3(past_r3 == 0) = r3;
past_r4(past_r4 == 0) = r4;
past_delta1(past_delta1 == 0) = d1;
past_delta2(past_delta2 == 0) = d2;
past_delta3(past_delta3 == 0) = d3;
past_delta4(past_delta4 == 0) = d4;
past_p1(past_p1 == 0) = p1;
past_p2(past_p2 == 0) = p2;
past_p3(past_p3 == 0) = p3;
past_p4(past_p4 == 0) = p4;
past_gamma(past_gamma == 0) = gamma;
past_pvax(past_pvax == 0) = p_vax;

past_beta1 = past_beta1(:,end-T:end);
past_beta2 = past_beta2(:,end-T:end);
past_beta3 = past_beta3(:,end-T:end);
past_beta4 = past_beta4(:,end-T:end);
past_r1 = past_r1(:,end-T:end);
past_r2 = past_r2(:,end-T:end);
past_r3 = past_r3(:,end-T:end);
past_r4 = past_r4(:,end-T:end);
past_delta1 = past_delta1(:,end-T:end);
past_delta2 = past_delta2(:,end-T:end);
past_delta3 = past_delta3(:,end-T:end);
past_delta4 = past_delta4(:,end-T:end);
past_p1 = past_p1(:,end-T:end);
past_p2 = past_p2(:,end-T:end);
past_p3 = past_p3(:,end-T:end);
past_p4 = past_p4(:,end-T:end);
past_gamma = past_gamma(:,end-T:end);
past_pvax = past_pvax(:,end-T:end);
%======================================================
% 4-stage SRK Tableus form %Refer Equation 4.55
A = [0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0]; %A
B_1 = [0 0 0 0;-0.7242916356 0 0 0;0.4237353406 -0.1994437050 0 0;-1.578475506 0.840100343 1.738375163 0]; %B^(1)
B_2 = [0 0 0 0;2.702000410 0 0 0;1.757261649 0 0 0;-2.918524118 0 0 0]; %B^(2)
G_1 = [-.7800788474;0.07363768240;1.486520013;0.2199211524]; %GAMMA 1
G_2 = [1.693950844;1.636107882;-3.024009558;-0.3060491602]; %GAMMA 2
% alpha =[1/6; 1/3; 1/3; 1/6]; %ALPHA
%======================================================
% plot(Datay(1,:),Datay(2,:),'ro');hold on %plot data
% plot(Datax(1,:),Datax(2,:),'bo');hold on %plot data
%======================================================
% plot(Datay(1,:),Datay(2,:),'ro');hold on %plot data
% plot(Datax(1,:),Datax(2,:),'bo');hold on %plot data
%======================================================
%SRK-4 Method
Ssrk4_2 = zeros(M,N);
I1srk4_2 = zeros(M,N);
I2srk4_2 = zeros(M,N);
I3srk4_2 = zeros(M,N);
I4srk4_2 = zeros(M,N);
Rsrk4_2 = zeros(M,N);
Dsrk4_2 = zeros(M,N);
Stemp = Szero.*ones(M,1);
I1temp = I1zero.*ones(M,1);
I2temp = I2zero.*ones(M,1);
I3temp = I3zero.*ones(M,1);
I4temp = I4zero.*ones(M,1);
Rtemp = Rzero.*ones(M,1);
Dtemp = Dzero.*ones(M,1);
SdW1 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
SdW2 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I1dW1 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I1dW2 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I2dW1 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I2dW2 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I3dW1 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I3dW2 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I4dW1 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
I4dW2 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
RdW1 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
RdW2 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
DdW1 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments
DdW2 = sqrt(dt).*randn(M,N);       %Numerical Brownian increments

for j = 1:L
    p_vax = past_pvax(j)/(past_p1(j) + past_p2(j) + past_p3(j) + past_p4(j) + past_pvax(j));
    p1 = past_p1(j)/(past_p1(j) + past_p2(j) + past_p3(j) + past_p4(j) + past_pvax(j));
    p2 = past_p2(j)/(past_p1(j) + past_p2(j) + past_p3(j) + past_p4(j) + past_pvax(j));
    p3 = past_p3(j)/(past_p1(j) + past_p2(j) + past_p3(j) + past_p4(j) + past_pvax(j));
    p4 = past_p4(j)/(past_p1(j) + past_p2(j) + past_p3(j) + past_p4(j) + past_pvax(j));
    
    beta1 = past_beta1(j);
    beta2 = past_beta2(j);
    beta3 = past_beta3(j);
    beta4 = past_beta4(j);
    r1 = past_r1(j);
    r2 = past_r2(j);
    r3 = past_r3(j);
    r4 = past_r4(j);
    d1 = past_delta1(j); 
    d2 = past_delta2(j); 
    d3 = past_delta3(j); 
    d4 = past_delta4(j); 
    gamma = past_gamma(j);
    
    SWinc1 = sum(SdW1(:,R*(j-1)+1:R*j),2);
    SWinc2 = sum(SdW2(:,R*(j-1)+1:R*j),2);
    I1Winc1 = sum(I1dW1(:,R*(j-1)+1:R*j),2);
    I1Winc2 = sum(I1dW2(:,R*(j-1)+1:R*j),2);
    I2Winc1 = sum(I2dW1(:,R*(j-1)+1:R*j),2);
    I2Winc2 = sum(I2dW2(:,R*(j-1)+1:R*j),2);
    I3Winc1 = sum(I3dW1(:,R*(j-1)+1:R*j),2);
    I3Winc2 = sum(I3dW2(:,R*(j-1)+1:R*j),2);
    I4Winc1 = sum(I4dW1(:,R*(j-1)+1:R*j),2);
    I4Winc2 = sum(I4dW2(:,R*(j-1)+1:R*j),2);
    RWinc1 = sum(RdW1(:,R*(j-1)+1:R*j),2);
    RWinc2 = sum(RdW2(:,R*(j-1)+1:R*j),2);
    DWinc1 = sum(DdW1(:,R*(j-1)+1:R*j),2);
    DWinc2 = sum(DdW2(:,R*(j-1)+1:R*j),2);
    
    SJ10 = (0.5*(Dt^(1/2))*(SWinc1+(1/sqrt(3)).*SWinc2));
    I1J10 = (0.5*(Dt^(1/2))*(I1Winc1+(1/sqrt(3)).*I1Winc2));
    I2J10 = (0.5*(Dt^(1/2))*(I2Winc1+(1/sqrt(3)).*I2Winc2));
    I3J10 = (0.5*(Dt^(1/2))*(I3Winc1+(1/sqrt(3)).*I3Winc2));
    I4J10 = (0.5*(Dt^(1/2))*(I4Winc1+(1/sqrt(3)).*I4Winc2));
    RJ10 = (0.5*(Dt^(1/2))*(RWinc1+(1/sqrt(3)).*RWinc2));
    DJ10 = (0.5*(Dt^(1/2))*(DWinc1+(1/sqrt(3)).*DWinc2));
    
    S1 = Stemp;
    I11 = I1temp;
    I21 = I2temp;
    I31 = I3temp;
    I41 = I4temp;
    R1 = Rtemp;
    
    %Refer Equation 4.56
    f1_S = mu*(S1+I11+I21+I31+141+R1)-S1./Ptotal.*p_vax.*(beta1.*I11+beta2.*I21+beta3.*I31+beta4.*I41)+gamma.*R1 - mu*S1;
    g1_S = -sigma1.*I11.*S1.*p_vax./Ptotal;
    f1_I1 = p1.*S1./Ptotal.*(beta1.*I11+beta2.*I21+beta3.*I31+beta4.*I41)-(r1+d1+mu).*I11;
    g1_I1 = sigma1.*p1.*I11.*S1./Ptotal;
    f1_I2 = p2.*S1./Ptotal.*(beta1.*I11+beta2.*I21+beta3.*I31+beta4.*I41)-(r2+d2+mu).*I21;
    g1_I2 = sigma1.*p2.*I11.*S1./Ptotal;
    f1_I3 = p3.*S1./Ptotal.*(beta1.*I11+beta2.*I21+beta3.*I31+beta4.*I41)-(r3+d3+mu).*I31;
    g1_I3 = sigma1.*p3.*I11.*S1./Ptotal;
    f1_I4 = p4.*S1./Ptotal.*(beta1.*I11+beta2.*I21+beta3.*I31+beta4.*I41)-(r4+d4+mu).*I41;
    g1_I4 = sigma1.*p4.*I11.*S1./Ptotal;
    f1_R = r1.*I11+r2.*I21+r3.*I31+r4.*I41 - mu*R1;
    g1_R = 0;
    f1_D = d1.*I11+d2.*I21+d3.*I31+d4.*I41;
    g1_D = 0;
    
    S2 = Stemp + Dt*(A(2)*f1_S) + (B_1(2).*SWinc1 + B_2(2).*SJ10).*g1_S;
    I12 = I1temp + Dt*(A(2)*f1_I1) + (B_1(2).*I1Winc1 + B_2(2).*I1J10).*g1_I1;
    I22 = I2temp + Dt*(A(2)*f1_I2) + (B_1(2).*I2Winc1 + B_2(2).*I2J10).*g1_I2;
    I32 = I3temp + Dt*(A(2)*f1_I3) + (B_1(2).*I3Winc1 + B_2(2).*I3J10).*g1_I3;
    I42 = I4temp + Dt*(A(2)*f1_I4) + (B_1(2).*I4Winc1 + B_2(2).*I4J10).*g1_I4;
    R2 = Rtemp + Dt*(A(2)*f1_R) + (B_1(2).*RWinc1 + B_2(2).*RJ10).*g1_R;
    
    %Refer Equation 4.56
    f2_S = mu*(S2+I12+I22+I32+142+R2)-S2./Ptotal.*p_vax.*(beta1.*I12+beta2.*I22+beta3.*I32+beta4.*I42)+gamma.*R2 - mu*S2;
    g2_S = -sigma1.*I12.*S2.*p_vax./Ptotal;
    f2_I1 = p1.*S2./Ptotal.*(beta1.*I12+beta2.*I22+beta3.*I32+beta4.*I42)-(r1+d1+mu).*I12;
    g2_I1 = sigma1.*p1.*I12.*S2./Ptotal;
    f2_I2 = p2.*S2./Ptotal.*(beta1.*I12+beta2.*I22+beta3.*I32+beta4.*I42)-(r2+d2+mu).*I22;
    g2_I2 = sigma1.*p2.*I12.*S2./Ptotal;
    f2_I3 = p3.*S2./Ptotal.*(beta1.*I12+beta2.*I22+beta3.*I32+beta4.*I42)-(r3+d3+mu).*I32;
    g2_I3 = sigma1.*p3.*I12.*S2./Ptotal;
    f2_I4 = p4.*S2./Ptotal.*(beta1.*I12+beta2.*I22+beta3.*I32+beta4.*I42)-(r4+d4+mu).*I42;
    g2_I4 = sigma1.*p4.*I12.*S2./Ptotal;
    f2_R = r1.*I12+r2.*I22+r3.*I32+r4.*I42 - mu*R2;
    g2_R = 0;
    f2_D = d1.*I12+d2.*I22+d3.*I32+d4.*I42;
    g2_D = 0;
    S3 = Stemp + Dt*(A(7)*f2_S) + (B_1(3).*SWinc1 + B_2(3).*SJ10).*g2_S + B_1(7)*SJ10.*g2_S;
    I13 = I1temp + Dt*(A(7)*f2_I1) + (B_1(3).*I1Winc1 + B_2(3).*I1J10).*g2_I1 + B_1(7)*I1J10.*g2_I1;
    I23 = I2temp + Dt*(A(7)*f2_I2) + (B_1(3).*I2Winc1 + B_2(3).*I2J10).*g2_I2 + B_1(7)*I2J10.*g2_I2;
    I33 = I3temp + Dt*(A(7)*f2_I3) + (B_1(3).*I3Winc1 + B_2(3).*I3J10).*g2_I3 + B_1(7)*I3J10.*g2_I3;
    I43 = I4temp + Dt*(A(7)*f2_I4) + (B_1(3).*I4Winc1 + B_2(3).*I4J10).*g2_I4 + B_1(7)*I4J10.*g2_I4;
    R3 = Rtemp + Dt*(A(7)*f2_R) + (B_1(3).*RWinc1 + B_2(3).*RJ10).*g2_R + B_1(7)*RJ10.*g2_R;
    
    %Refer Equation 4.56
    f3_S = mu*(S3+I13+I23+I33+143+R3)-S3./Ptotal.*p_vax.*(beta1.*I13+beta2.*I23+beta3.*I33+beta4.*I43)+gamma.*R3-mu*S3;
    g3_S = -sigma1.*I13.*p_vax.*S3./Ptotal;
    f3_I1 = p1.*S3./Ptotal.*(beta1.*I13+beta2.*I23+beta3.*I33+beta4.*I43)-(r1+d1+mu).*I13;
    g3_I1 = sigma1.*p1.*I13.*S3./Ptotal;
    f3_I2 = p2.*S3./Ptotal.*(beta1.*I13+beta2.*I23+beta3.*I33+beta4.*I43)-(r2+d2+mu).*I23;
    g3_I2 = sigma1.*p2.*I13.*S3./Ptotal;
    f3_I3 = p3.*S3./Ptotal.*(beta1.*I13+beta2.*I23+beta3.*I33+beta4.*I43)-(r3+d3+mu).*I33;
    g3_I3 = sigma1.*p3.*I13.*S3./Ptotal;
    f3_I4 = p4.*S3./Ptotal.*(beta1.*I13+beta2.*I23+beta3.*I33+beta4.*I43)-(r4+d4+mu).*I43;
    g3_I4 = sigma1.*p4.*I13.*S3./Ptotal;
    f3_R = r1.*I13+r2.*I23+r3.*I33+r4.*I43-mu*R3;
    g3_R = 0;
    f3_D = d1.*I13+d2.*I23+d3.*I33+d4.*I43;
    g3_D = 0;
    S4 = Stemp + Dt*(f3_S) + (B_1(4).*SWinc1 + B_2(4).*SJ10).*g3_S + (B_1(8).*SWinc1.*g3_S + B_1(12).*SWinc1.*g3_S);
    I14 = I1temp + Dt*(f3_I1) + (B_1(4).*I1Winc1 + B_2(4).*I1J10).*g3_I1 + (B_1(8).*I1Winc1.*g3_S + B_1(12).*I1Winc1.*g3_I1);
    I24 = I2temp + Dt*(f3_I2) + (B_1(4).*I2Winc1 + B_2(4).*I2J10).*g3_I2 + (B_1(8).*I2Winc1.*g3_I2 + B_1(12).*I2Winc1.*g3_I2);
    I34 = I3temp + Dt*(f3_I3) + (B_1(4).*I3Winc1 + B_2(4).*I3J10).*g3_I3 + (B_1(8).*I3Winc1.*g3_I3 + B_1(12).*I3Winc1.*g3_I3);
    I44 = I4temp + Dt*(f3_I4) + (B_1(4).*I4Winc1 + B_2(4).*I4J10).*g3_I4 + (B_1(8).*I4Winc1.*g3_I4 + B_1(12).*I4Winc1.*g3_I4);
    R4 = Rtemp + Dt*(f3_R) + (B_1(4).*RWinc1 + B_2(4).*RJ10).*g3_R + (B_1(8).*RWinc1.*g3_R + B_1(12).*RWinc1.*g3_R);
    
    %Refer Equation 4.56
    f4_S = mu*(S4+I14+I24+I34+144+R4)-S4./Ptotal.*p_vax.*(beta1.*I14+beta2.*I24+beta3.*I34+beta4.*I44)+gamma.*R4-mu*S4;
    g4_S = -sigma1.*I14.*p_vax.*S4./Ptotal;
    f4_I1 = p1.*S4./Ptotal.*(beta1.*I14+beta2.*I24+beta3.*I34+beta4.*I44)-(r1+d1+mu).*I14;
    g4_I1 = sigma1.*p1.*I14.*S4./Ptotal;
    f4_I2 = p2.*S4./Ptotal.*(beta1.*I14+beta2.*I24+beta3.*I34+beta4.*I44)-(r2+d2+mu).*I24;
    g4_I2 = sigma1.*p2.*I14.*S4./Ptotal;
    f4_I3 = p3.*S4./Ptotal.*(beta1.*I14+beta2.*I24+beta3.*I34+beta4.*I44)-(r3+d3+mu).*I34;
    g4_I3 = sigma1.*p3.*I14.*S4./Ptotal;
    f4_I4 = p4.*S4./Ptotal.*(beta1.*I14+beta2.*I24+beta3.*I34+beta4.*I44)-(r4+d4+mu).*I44;
    g4_I4 = sigma1.*p4.*I14.*S4./Ptotal;
    f4_R = r1.*I14+r2.*I24+r3.*I34+r4.*I44 -mu*R4;
    g4_R = 0;
    f4_D = d1.*I14+d2.*I24+d3.*I34+d4.*I44;
    g4_D = 0;
        
    Stemp = Stemp + Dt.*((1/6)*f1_S + (1/3).*f2_S + (1/3).*f3_S + (1/6).*f4_S) + (G_1(1).*SWinc1 + G_2(1).*SJ10).*g1_S + (G_1(2).*SWinc1 + G_2(2).*SJ10).*g2_S + (G_1(3).*SWinc1 + G_2(3).*SJ10).*g3_S + (G_1(4).*SWinc1 + G_2(4).*SJ10).*g4_S;
    I1temp = I1temp + Dt.*((1/6)*f1_I1 + (1/3).*f2_I1 + (1/3).*f3_I1 + (1/6).*f4_I1) + (G_1(1).*I1Winc1 + G_2(1).*I1J10).*g1_I1 + (G_1(2).*I1Winc1 + G_2(2).*I1J10).*g2_I1 + (G_1(3).*I1Winc1 + G_2(3).*I1J10).*g3_I1 + (G_1(4).*I1Winc1 + G_2(4).*I1J10).*g4_I1;
    I2temp = I2temp + Dt.*((1/6)*f1_I2 + (1/3).*f2_I2 + (1/3).*f3_I2 + (1/6).*f4_I2) + (G_1(1).*I2Winc1 + G_2(1).*I2J10).*g1_I2 + (G_1(2).*I2Winc1 + G_2(2).*I2J10).*g2_I2 + (G_1(3).*I2Winc1 + G_2(3).*I2J10).*g3_I2 + (G_1(4).*I2Winc1 + G_2(4).*I2J10).*g4_I2;
    I3temp = I3temp + Dt.*((1/6)*f1_I3 + (1/3).*f2_I3 + (1/3).*f3_I3 + (1/6).*f4_I3) + (G_1(1).*I3Winc1 + G_2(1).*I3J10).*g1_I3 + (G_1(2).*I3Winc1 + G_2(2).*I3J10).*g2_I3 + (G_1(3).*I3Winc1 + G_2(3).*I3J10).*g3_I3 + (G_1(4).*I3Winc1 + G_2(4).*I3J10).*g4_I3;
    I4temp = I4temp + Dt.*((1/6)*f1_I4 + (1/3).*f2_I4 + (1/3).*f3_I4 + (1/6).*f4_I4) + (G_1(1).*I4Winc1 + G_2(1).*I4J10).*g1_I4 + (G_1(2).*I4Winc1 + G_2(2).*I4J10).*g2_I4 + (G_1(3).*I4Winc1 + G_2(3).*I4J10).*g3_I4 + (G_1(4).*I4Winc1 + G_2(4).*I4J10).*g4_I4;
    Rtemp = Rtemp + Dt.*((1/6)*f1_R + (1/3).*f2_R + (1/3).*f3_R + (1/6).*f4_R) + (G_1(1).*RWinc1 + G_2(1).*RJ10).*g1_R + (G_1(2).*RWinc1 + G_2(2).*RJ10).*g2_R + (G_1(3).*RWinc1 + G_2(3).*RJ10).*g3_R + (G_1(4).*RWinc1 + G_2(4).*RJ10).*g4_R;
    Dtemp = Dtemp + Dt.*((1/6)*f1_D + (1/3).*f2_D + (1/3).*f3_D + (1/6).*f4_D) + (G_1(1).*DWinc1 + G_2(1).*DJ10).*g1_D + (G_1(2).*DWinc1 + G_2(2).*DJ10).*g2_D + (G_1(3).*DWinc1 + G_2(3).*DJ10).*g3_D + (G_1(4).*DWinc1 + G_2(4).*DJ10).*g4_D;

    %Refer Equation 4.56
    Ssrk4_2(:,j) = Stemp.';
    I1srk4_2(:,j) = I1temp.';
    I2srk4_2(:,j) = I2temp.';
    I3srk4_2(:,j) = I3temp.';
    I4srk4_2(:,j) = I4temp.';
    Rsrk4_2(:,j) = Rtemp.';
    Dsrk4_2(:,j) = Dtemp.';
end
% RMSE_srk4x = rmse(Xsrk4_2,Datax);
% RMSE_srk4y = rmse(Ysrk4_2,Datay);
Save_srk4 = mean(Ssrk4_2,1);
I1ave_srk4 = mean(I1srk4_2,1);
I2ave_srk4 = mean(I2srk4_2,1);
I3ave_srk4 = mean(I3srk4_2,1);
I4ave_srk4 = mean(I4srk4_2,1);
Rave_srk4 = mean(Rsrk4_2,1);
Dave_srk4 = mean(Dsrk4_2,1);
% Smodel = Save_srk4;
Smodel = Save_srk4;
I1model = I1ave_srk4;
I2model = I2ave_srk4;
I3model = I3ave_srk4;
I4model = I4ave_srk4;
Rmodel = Rave_srk4;
Dmodel = Dave_srk4;
end
%======================================================