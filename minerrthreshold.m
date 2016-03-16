function [bw,T]=minerrthreshold(img)
% [bw,T]=minerrthreshold(img,lvl)
% MINERRTHRESHOLD segments a grayscale images & identifies an 'optimal threshold'
% as defined as the threshold with the minimum error thresholding 
% all the math is from the paper:
% 
% Kittler, J., and J. Illingworth. 
% â€œMinimum error thresholding.â€? 
% Pattern Recognition 19, no. 1 (1986): 41-47.
%
% Algorithm summary:
% the basic idea is that they defined a score funciton for misclassification. 
% All I do here is define this function as a function ofthe threshold and 
% use fminsearch to find this threshold. 
%
% The modification by Cho et al improves the estimate of the variance to
% take into avount the fact that it is trauncated
% 
% inputs: 
%   lvl - how many gray levels are there (defaults 4096)
% outputs: 
%   bw - segmented image
%   

% move img to be only integer
% keep it of class double for convenience
switch class(img)
    case 'uint8'
        lvl=256;
    case 'uint16'
        lvl=2^16;
    case 'double'
        lvl=length(unique(img(:)));
        lvl=min(lvl,2^16);
        img=gray2ind(img,lvl);

end

img=double(img);


% create the histogram for the image
[H,bins]=hist(img(:),lvl);
H=H./sum(H);

% Define all accesory functions  
Pb=@(T) sum(H(bins<T));
Pf=@(T) sum(H(bins>=T));
mub=@(T) 1/Pb(T)*sum(H(bins<T).*bins(bins<T));
muf=@(T) 1/Pf(T)*sum(H(bins>=T).*bins(bins>=T));
sigb=@(T) sqrt(1/Pb(T)*sum(H(bins<T).*(bins(bins<T)-mub(T)).^2));
sigf=@(T) sqrt(1/Pf(T)*sum(H(bins>=T).*(bins(bins>=T)-muf(T)).^2));

% define the function to minimize 
J=@(T) 1+2*(Pb(T)*log(sigb(T))+Pf(T)*log(sigf(T)))-2*(Pb(T)*log(Pb(T))+Pf(T)*log(Pf(T)));

% initial guess
T0=graythresh(img/lvl)*lvl;

% brute force checking all levels
[T,bla,flag]=fminsearch(J,T0);
if flag==1
    bw=img>T;
    T=T/lvl;
else
    bw =zeros(size(img));
    T = 1;
end