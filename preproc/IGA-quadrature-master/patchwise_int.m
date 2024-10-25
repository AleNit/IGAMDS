
function [IPu,IPv,nqpu,nqpv] = patchwise_int(U,V,CP,over,order,regularity)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

addpath 'IGA-quadrature-master';

[pu,wu]=IGA_quadrature_master(unique(U),over,order,regularity);
nu = length(pu);
ue = zeros(nu,1);
for ii = 1:length(pu)
    ue(ii) = findspan(pu(ii),U,length(CP(:,1,1)));  
end
IPu = [pu,wu,ue];

[pv,wv]=IGA_quadrature_master(unique(V),over,order,regularity);
nv = length(pv);
ve = zeros(nv,1);
for ii = 1:length(pv)
    ve(ii) = findspan(pv(ii),V,length(CP(1,:,1)));  
end
IPv = [pv,wv,ve];

% compute array of quadrature point count per element
uU = unique(IPu(:,3));
nqpu = histc(IPu(:,3),uU);

uV = unique(IPv(:,3));
nqpv = histc(IPv(:,3),uV);

end

