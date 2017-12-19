%% UssbPi
%returns DPSWFs in frequency domain of Slepians with single sideband
%modulation

%pulse is normalized so that the integral of its absolute value in the time
%domain is Pi 

%% input
%N=number of segments
%W=half-width of band
%K=[k1,...,kn] row vector containing Slepian orders that were used in the multitaper
%w0=carrier frequency
%Wtot=frequency range of data
%Ninc= number of increments used in integration for biases
%delT=time increments

%% output
%WF is a 2Ninc+1 x length(K) sized array. Each column is a DPSWF corresponding
%to the Slepian ki\in K. Rows are the DPSWF computed at frequencies 
%(pi/delT)*(-Ninc:Ninc)/Ninc 



function WF=UssbPi(N,W,w0,Wtot,delT,K,Ninc)

nk=length(K);


%Get dpss sequences
M=dpss(N,N*W);
MK=M(:,K+1);

%Hilbert transform
H=hilbert(MK);


%shift by carrier frequency
c=arrayfun(@(n)(cos(w0*n*delT)),transpose(0:(N-1)));
s=arrayfun(@(n)(sin(w0*n*delT)),transpose(0:(N-1)));

shiftc=repmat(c,1,nk);
shifts=repmat(s,1,nk);


%amplitude modulation 

%The real part of the output of
%the "hilbert" function is the unmodulated sequence. The imaginary part is
%the Hilbert transform of the unmodulated sequence.
vssb=real(H).*shiftc-imag(H).*shifts;


% normalize so that the area underneath all modulated Slepians is Pi.
norms=arrayfun(@(k)(trapz((0:(N-1))*delT,abs(transpose(vssb(:,k))))),1:nk);
norm_array=repmat(norms,N,1);
vssb=vssb.*(pi*(norm_array.^(-1)));    

%% wave functions
[X,Y]=meshgrid(0:(N-1),(Wtot)*(0:Ninc)/Ninc);

%symmetric fourier transform
DFTc=arrayfun(@(x,y)(exp(1i*y*delT*(x-(N-1)/2))),X,Y);


%% output
WF=cell2mat(arrayfun(@(k)(DFTc*vssb(:,k)),1:nk,'UniformOutput',false));

end