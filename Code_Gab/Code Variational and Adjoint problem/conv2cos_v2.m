function [vwshort,vwfull,vwextended] = conv2cos_v2(v,w)
% computes the convolution for cosine series V and W 
% Works for same size V and W
% works in all dimensions
% vwpfull is just the full cosine convolution
% vwpextended is are all modes, not just the positive ones
% vwpshort has size(v)

N = size(v)-1;
dim=length(N);

p0=2; % takes care of p=0 case

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  sfull{j}=p0*N(j)+1:2*p0*N(j)+1;
end

vv=v(sflip{1:dim});
ww=w(sflip{1:dim});


[~,vwextended]=convtensor(vv,ww);
vwextended=real(vwextended);

vwfull = vwextended(sfull{1:dim});
vwshort = vwfull(s{1:dim});

return
