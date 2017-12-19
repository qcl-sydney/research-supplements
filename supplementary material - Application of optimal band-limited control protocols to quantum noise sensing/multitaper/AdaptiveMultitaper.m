%%  weightsfinal 
% produces multitaper estimate using adaptive weighting

%% output: 

%o1=S_est, colums of this array are estimates of the spectrum produced at 
%each iteration of the adaptive weighting procedure

%o2=D_array, a length(shifts) x length(K) x (number of iterations) sized
%array that contains the weights produced by each iteration of the adaptive
%weighting procedure. Each entry D_Array_{i,j,k} is the weight of the
%estimate for Slepian K(j) modulated by carrier shifts(i)at iteration k.

%o3=E_array, a length(shifts) x length(K) matrix, where entry E_array_{ij} is 
% the estimate produced by the Slepian of order K(j) modulated by 
% carrier shifts(i)



function [o1,o2,o3]=AdaptiveMultitaper()



%% input parameters
load('data07272016.mat','pBright');   
%pBright is a matrix with dimensions length(shifts) x length(K), and
% elements p_{w_0,k}= P(bright) for Slepian k modulated by ssb 
%with carrier w_0             
              
K=[1 3 5 7];  %list of Slepian orders used

N=2240;  %number of points in Slepians

W=7/2240; %dimensionless width of band

delT=10^(-6); %seconds.  length of increments. 

Wtot=14*10^4; %Hz. Range of reconstruction.

Ninc=1000;  %number of steps in integrations to determine biases

shifts=(pi/delT)*(10/2240)*(0:.5:4);    %Hz. list of carrier frequencies



%% derived quantities   
Winc=Wtot*(0:Ninc)/Ninc; %list of Ninc+1 frequencies from 0 to Wtot, 
                         %spaced by Wtot/Ninc;
                         
inc=(Wtot/Ninc);   %size of frequency spacing in Winc                      

HB_start=round((shifts/Wtot)*Ninc)+1; %list of elements in Winc where the 
                                   %half-bands begin for each carrier

HBN=round(((2*W*pi/delT)/Wtot)*Ninc); %number of frequencies in Winc enclosed 
                                    %by the half-band

HB_end=HB_start+HBN; %list of elements in Winc where the 
                     %half-bands end for each carrier



%% function that produces interpolated spectrum

%This function takes an estimate of the spectrum or the derivative of the spectrum 
%at the centers of the
%half-bands and produces a linear interpolation on the 
%frequencies Winc. The interpolation is used in the numerical integration
%to determine the biases.
%At frequencies greater than the largest center half-band,
%the interpolation is zero.

function int=Sinterp(est)
wz=Winc(max(HB_start)+ceil(HBN/2)):inc:Wtot;
z=zeros(1,length(wz));
    
int=interp1(cat(2,min(shifts),shifts+(pi/delT)*(W/2),wz),cat(2,est(1),transpose(est),z),...
    Winc,'linear');
end
    

                           
     

%% ssb modulated DPSWFs   
% get ssb modulated DPSWFs for all carriers and Slepian orders

% output is a cell "S_array" 
% entries of cell correspond to the carrier frequencies in "shifts" 
% each entry is a length(Winc) x length(K) matrix. The columns of these
% matrices are the DPSWFs for the Slepians K(1), K(2),... at the 
% frequencies in Winc

S_array=arrayfun(@(w0)(UssbPi(N,W,w0,Wtot,delT,K,Ninc)),shifts,...
    'UniformOutput',false);

%% C "normalization"(area under Slepians in half-band)
%this is an length(shifts) x length(K) matrix, where entry C_array_{ij} is 
%the area of the FF produced by the Slepian of order K(j) in the half-band 
% of carrier shifts(i) 
[X,Y]=meshgrid(1:length(K),1:length(shifts));
C_array=arrayfun(@(k,c)(trapz(Winc(HB_start(c):HB_end(c)),(((4/pi)*sin(delT*(Winc(HB_start(c):HB_end(c))+.1)/2).^2)./...
    ((Winc(HB_start(c):HB_end(c))+.1).^2)).*abs(transpose(S_array{1,c}(HB_start(c):HB_end(c),k))).^2)),X,Y);


%% Estimates
% this is a length(shifts) x length(K) matrix, where entry E_array_{ij} is 
% the estimate produced by the Slepian of order K(j) modulated by 
% carrier shifts(i) 
E_array=(1-pBright)./C_array;
S0=E_array(:,1);


%% Estimated broad band bias
% this function produces an length(shifts) x length(K) matrix, where entry 
% B_array_{ij} is  the estimated broadband bias by the Slepian of order 
% K(j) modulated by carrier shifts(i) based the spectral estimate "esti"
function b=B_array(esti)
est=Sinterp(esti);


% %estimated upperbound of broadband bias
% bi=(cat(1,zeros(1,length(K)),arrayfun(@(k,c)(trapz(Winc(1:HB_start(c)),(4/pi)*(sin(delT*(Winc(1:HB_start(c))+.1)/2).^2)./...
%     ((Winc(1:HB_start(c))+.1).^2).*(abs(transpose(S_array{1,c}((1:HB_start(c)),k))).^2))),...
%     X(2:length(shifts),:),Y(2:length(shifts),:)))+...
%     arrayfun(@(k,c)(trapz(Winc(HB_end(c):Ninc),(4/pi)*(sin(delT*(Winc(HB_end(c):Ninc)+.1)/2).^2)./...
%     ((Winc(HB_end(c):Ninc)+.1).^2).*(abs(transpose(S_array{1,c}(HB_end(c):Ninc,k))).^2))),...
%     X,Y))*trapz(Winc,est);

%estimated broadband bias
bi=cat(1,zeros(1,length(K)),arrayfun(@(k,c)(trapz(Winc(1:HB_start(c)),(4/pi)*(sin(delT*(Winc(1:HB_start(c))+.1)/2).^2)./...
    ((Winc(1:HB_start(c))+.1).^2).*est(1:HB_start(c)).*(abs(transpose(S_array{1,c}((1:HB_start(c)),k))).^2))),...
    X(2:length(shifts),:),Y(2:length(shifts),:)))+...
    arrayfun(@(k,c)(trapz(Winc(HB_end(c):Ninc),(4/pi)*(sin(delT*(Winc(HB_end(c):Ninc)+.1)/2).^2)./...
    ((Winc(HB_end(c):Ninc)+.1).^2).*est(HB_end(c):Ninc).*(abs(transpose(S_array{1,c}(HB_end(c):Ninc,k))).^2))),...
    X,Y);

b=bi./C_array;
end

%% Estimated local bias 
% This function takes an estimate of the spectrum and takes the finite 
% difference to estimate S'(\omega) on Winc
function d=D1(est)
wd=transpose(shifts+(pi/delT)*(W/2));
est_l=length(est);
di=(est(2:est_l)-est(1:(est_l-1)))./(wd(2:est_l)-wd(1:(est_l-1)));
d=cat(1,di, di(est_l-1));    
end



% this function produces an length(shifts) x length(K) matrix, where entry 
% L_array_{ij} is  the estimated local bias by the Slepian of order 
% K(j) modulated by carrier shifts(i) based the spectral estimate "esti"
function l=L_array(esti)
d1_est=Sinterp(D1(esti));
li=arrayfun(@(k,c)(trapz(Winc(HB_start(c):HB_end(c)),(((4/pi)*sin(delT*(Winc(HB_start(c):HB_end(c))+.1)/2).^2)./...
    ((Winc(HB_start(c):HB_end(c))+.1).^2)).*(abs(transpose(S_array{1,c}(HB_start(c):HB_end(c),k))).^2).*...
    d1_est(HB_start(c):HB_end(c)).*(Winc(HB_start(c):HB_end(c))-(shifts(c)+(pi/delT)*(W/2))))),X,Y);
l=li./C_array;
end 



%% Estimated infidelity
% this function returns an length(shifts) x length(K) matrix, where entry 
% I_array_{ij} is  the estimated infidelity based the spectral estimate 
% "esti" when the Slepian pulse of order K(j) modulated by 
% carrier shifts(i) is applied to the qubit

function l=I_array(esti)
est=Sinterp(esti);
l=arrayfun(@(k,c)(trapz(Winc,(4/pi)*(sin(delT*(Winc+.1)/2).^2)./...
    ((Winc+.1).^2).*est.*(abs(transpose(S_array{1,c}(:,k))).^2))),X,Y);
end 



%% Adaptive weighting
it=10; %number of interations
S_est=zeros(length(shifts),it+1); %array to store spectral estimates at
                                  %each iteration
S_est(:,1)=S0; % initial spectral estimate (the estimate produced by the 
               % order 1 Slepian)
D_array=zeros(length(shifts),length(K),it+1); % array to store weighting 
                                              % coefficients at each
                                              % iteration

for i=1:it
Si=S_est(:,i);
Si_int=Sinterp(Si);

Sc=arrayfun(@(c)(Si_int(round(HB_start(c)+HBN/2))),Y); % estimate of 
                                                          % spectrum at the
                                                          % center
                                                          % half-bands

D_array(:,:,i)=(I_array(Si)./(C_array.* (Sc+B_array(Si)+L_array(Si))));
%solve for weights that minimize RMS error

Snew=sum(D_array(:,:,i).*E_array,2)./sum(D_array(:,:,i),2);
%compute new estimate of spectrum

S_est(:,i+1)=Snew;    
end 

o1=S_est;
o2=D_array;
o3=E_array;


end