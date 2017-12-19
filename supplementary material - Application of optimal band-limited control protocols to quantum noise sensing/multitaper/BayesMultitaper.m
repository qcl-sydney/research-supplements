function [out1, out2] = BayesMultitaper()

%% output
%out1=spectral estimates in each segment
%out2=center frequencies of each segment


%% Input
%load probability bright from experiment
load('data07272016.mat','pBright'); 




%% Parameters
%list of Slepian orders used
K=[1 3 5 7]; 

%number of ions
num_ions=10000;

%number of points in Slepians
N=2240;  

%dimensionless width of band
W=7/2240; 

%seconds.  length of increments.
delT=10^(-6);  

%Hz. list of carrier frequencies
shifts=(pi/delT)*(10/2240)*(0:.5:4);    

% centers of half-bands, estimation points of the multitaper reconstruction
om_est=shifts+pi*W/delT;  


%% Determine segments
%these are the segments at which the spectrum is estimated

% number segments per half-band. This cannot be made too large or else the
% reconstruction matrix will be ill-conditioned
segs_num=4;

%end points of half-bands
shifts_end=shifts+2*pi*W/delT;

%length of segments in half-bands
segs_length=(2*pi*W/delT)/segs_num;

%start and end points of segments
segs=cell2mat(arrayfun(@(c)(shifts(c):segs_length:shifts_end(c)),transpose(1:length(shifts)),'UniformOutput',false));

%centers of segments
segs_centers=cell2mat(arrayfun(@(c)(segs(c,1:segs_num)+segs_length/2),transpose(1:length(shifts)),'UniformOutput',false));


%% Prior
% prior mean
S_est=(10^(-4))*[0.3614,    0.3124,    0.3108,    0.3023,    0.3254,...
                 0.2432,    0.0300,    0.0230,    0.0202];
%interpolation of prior mean at the center of each segment
function s=S_prior(oms)
    s=interp1(cat(2,min(min(segs_centers)),om_est,max(max(segs_centers))),cat(2,S_est(1),S_est,S_est(length(shifts))),...
    oms,'linear');
end

%creates array of S_prior with rows corresponding to bands and columns 
%corresponding to segments
prior_array=cell2mat(arrayfun(@(c)(S_prior(segs_centers(c,:))),...
    transpose(1:length(shifts)),'UniformOutput',false));

%variance of prior 
var0=(1.5*10^(-5))^2;




%% Average value of FFs in the segments
% average values of FFs are needed in inversion matrix

% This is a 1D cell array (dimension is # of carrier frequencies x 1)
% Each entry in the cell array is a matrix with dimension # segments x # 
% Slepians overlapping segment
FF_array=arrayfun(@(c)(UssbBayes(N,W,shifts(c),delT,K,segs(c,:))),...
    transpose(1:length(shifts)),'UniformOutput',false);




%% Create linear inversion

%initialize output array
%dimensions are # carriers x # of segements per half-band
%each entry is the spectral estimate given by inversion
S_bayes=zeros(length(shifts),segs_num);


% elements of inversion matrix
function e=elt(i,l,p)   
    e=sum((segs_length*num_ions*var0/pi^2)*...
        arrayfun(@(k)((pBright(i,k)*(1-pBright(i,k)))^(-1)*...
        FF_array{i,1}(l,k)*FF_array{i,1}(p,k)),1:length(K)));
end


%produce estimates for the segments in each half-band
%each iteration of this loop corresponds to a half-band
for i=1:length(shifts)

% vector depending on experimental measurements 
% (the inverse of the inversion matrix is applied to this vector)    
vec=transpose(prior_array(i,:))+(num_ions*var0/pi)*FF_array{i,1}*...
    transpose(pBright(i,:).^(-1));

% inversion matrix
[X,Y]=meshgrid(1:segs_num);
mat=cell2mat(arrayfun(@(x,y)elt(i,x,y),X,Y,'UniformOutput',false))+...
    eye(segs_num);

% estimates for segments in half-band "i" are stored in S_bayes
S_bayes(i,:)=mat\vec;
       
end

%% Weighting 
% the segment estimates are weighted by the Fisher info they carry about 
% the spectrum



% put segments, shifts and estimates into a more useful form

% segs_all is an array where each row corresponds to a segment and the 
% columns ordered from 1 to 5 are 
%        1. the start frequency of the segment
%        2. the end frequency of the segment
%        3. the half-band containing the segment
%        4. estimate of spectrum in the segment
%        5. center frequency of the segment
segs_all=zeros(segs_num*length(shifts),5);
for i=1:length(shifts)
    for j=1:segs_num
        segs_all((i-1)*4+j,1)=segs(i,j);
        segs_all((i-1)*4+j,2)=segs(i,j+1);
        segs_all((i-1)*4+j,3)=i;
        segs_all((i-1)*4+j,4)=S_bayes(i,j);
        segs_all((i-1)*4+j,5)=segs_centers(i,j);
        
    end
end


% For each segment, find the segments that it overlaps

% overlaps is a cell array with each row corresponding to a segment. 

% The first column is a list of the "down overlaps" of each segment, the 
% overlapping segments centered at smaller frequencies 

% The second column is a list of the "up overlaps" of each segment, the 
% overlapping segments centered at larger frequencies 
 
overlaps=cell(segs_num*length(shifts),2);
for i=1:segs_num*length(shifts)
    up_overlaps=[];
    down_overlaps=[];
    for j=1:segs_num*length(shifts)
        if (segs_all(j,1)>segs_all(i,1))&&(segs_all(j,1)<segs_all(i,2))
           up_overlaps=cat(2,up_overlaps,j);
        end
        
        if (segs_all(j,2)>segs_all(i,1))&&(segs_all(j,2)<segs_all(i,2))
           down_overlaps=cat(2,down_overlaps,j);
        end
    end

    overlaps{i,1}=down_overlaps;
    overlaps{i,2}=up_overlaps;    
end

% S_final is an column vector where each entry is the final weighted
% spectral estimate in each segment
S_final=zeros(segs_num*length(shifts),1);

% this loop goes through all the segments and computes the Fisher info for
% the segment and all overlapping segments
for i=1:segs_num*length(shifts)
    
    if isempty(overlaps{i,1})==0
        down_list=zeros(length(overlaps{i,1}),2);
        %compute the Fisher info for each of the down overlaps
        for j=1:length(overlaps{i,1})
            segj=overlaps{i,1}(j);
            % numerator in expression for Fisher info for each Slepian in
            % the segment
            area=UssbFisher(N,W,shifts(segs_all(segj,3)),delT,K,[segs_all(i,1),segs_all(segj,2)]);
            % the weighting coefficients of the segment, i.e. Fisher info 
            % for each Slepian in the segment
            weights=(area.^2)./(pBright(segs_all(segj,3)).*(1-pBright(segs_all(segj,3))));
            % sum the fisher info for each of the Slepians
            down_list(j,1)=sum(weights);
             % unnormalized weighted estimate (the estimate in the segment
             % times the weight that depends on the Fisher info. In the
             % next step this unnormalized estimates will be summed
             % together and "normalized" by dividing by the sum of all
             % weights
            down_list(j,2)=sum(weights)*segs_all(segj,4);    
        end
    end
    
    
    if isempty(overlaps{i,2})==0
        up_list=zeros(length(overlaps{i,2}),2);
         %compute the Fisher info for each of the up overlaps
        for j=1:length(overlaps{i,2})
            segj=overlaps{i,2}(j);
            area=UssbFisher(N,W,shifts(segs_all(segj,3)),delT,K,[segs_all(segj,1),segs_all(i,2)]);
            weights=(area.^2)./(pBright(segs_all(segj,3)).*(1-pBright(segs_all(segj,3))));
            up_list(j,1)=sum(weights);
            up_list(j,2)=sum(weights)*segs_all(segj,4);    
        end
    end 
 
 % weight depending on Fisher info of segment itself   
 self_area=UssbFisher(N,W,shifts(segs_all(i,3)),delT,K,[segs_all(i,1),segs_all(i,2)]);
 self_weights=(self_area.^2)./(pBright(segs_all(i,3)).*(1-pBright(segs_all(i,3))));
    
 if (isempty(overlaps{i,1})==1)&&(isempty(overlaps{i,2})==1)
     % if there are no overlapping segments estimate of segment stays the
     % same
    S_final(i)=segs_all(i,4);
 
 elseif (isempty(overlaps{i,1})==0)&&(isempty(overlaps{i,2})==0)
     % if there are up overlaps & down overlaps, the final weighted 
     % estimate is the sum of the estimates in each of the overlapping
     % segments weighted by the Fisher info and divided by all the weights
     % to "normalize"
    S_final(i)=(self_weights*segs_all(i,4)+sum(up_list(:,2))+sum(down_list(:,2)))/...
    (self_weights+sum(up_list(:,1))+sum(down_list(:,1)));
 
 elseif (isempty(overlaps{i,1})==0)&&(isempty(overlaps{i,2})==1) 
     % final weighted estimate if there are only down overlaps
      S_final(i)=(self_weights*segs_all(i,4)+sum(down_list(:,2)))/...
    (self_weights+sum(down_list(:,1)));

 else
     % final weighted estimate if there are only up overlaps
     S_final(i)=(self_weights*segs_all(i,4)+sum(up_list(:,2)))/...
    (self_weights+sum(up_list(:,1)));
 end

end    

% sort the segments in order of increasing frequency
[B, I]=sort(segs_all(:,5));

out1=S_final(I); %final spectral estimate in each segment
out2=B; %final spectral estimate in each segment

plot(out2,out1);

end


























