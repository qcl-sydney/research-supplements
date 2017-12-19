%% UssbFisher
%integrated value of DPSS filters (for each k)  in each segment

%% input
%N=number of segments
%W=half-width of band
%K=[k1,...,kn] row vector containing Slepian orders that were used in the multitaper
%w0=carrier frequency
%delT= time increments
%oms= start and endpoints of segments in half-band starting at w0

%% output
%out= array 
%dimensions are (# segments per half-band)x nk 
%each entry is the integrated value of the filter corresponding to the DPSS of
%order k in the segment


function out=UssbFisher(N,W,w0,delT,K,oms)

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

%(matlab is strange. The real part of the output of
%the "hilbert" function is the unmodulated sequence. The imaginary part is
%the Hilbert transform of the unmodulated sequence.)
vssb=real(H).*shiftc-imag(H).*shifts;


% normalize so that the area underneath all modulated Slepians is Pi.
norms=arrayfun(@(k)(trapz((0:(N-1))*delT,abs(transpose(vssb(:,k))))),1:nk);
norm_array=repmat(norms,N,1);
vssb=vssb.*(pi*(norm_array.^(-1)));    

% integrated values
seg_avs=zeros((length(oms)-1),nk);
for i=1:(length(oms)-1)
inc=(oms(i+1)-oms(i))/100;
oms_i=oms(i):inc:oms(i+1);
[X,Y]=meshgrid(0:(N-1),oms_i);
%symmetric fourier transform
DFTc=arrayfun(@(x,y)(exp(1i*y*delT*(x-(N-1)/2))),X,Y);
%WF=cell2mat(arrayfun(@(k)(trapz(oms_i,abs(DFTc*vssb(:,k)).^2))/(oms(i+1)-oms(i)),1:nk,'UniformOutput',false));
WF=cell2mat(arrayfun(@(k)(trapz(transpose(oms_i),4*(abs(DFTc*vssb(:,k)).^2).*...
    (sin(delT*transpose(oms_i)/2).^2)./((transpose(oms_i).^2+.01)))),1:nk,...
    'UniformOutput',false));
seg_avs(i,:)=WF;
end

out=seg_avs;

end









