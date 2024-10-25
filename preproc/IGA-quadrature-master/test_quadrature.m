
%Knots=test_quadrature();

function Knots=test_quadrature()
%esempio
pars=[0,2,3,4,0,1];
Knots.k=[0 0.0500000000000000 0.100000000000000 0.150000000000000 0.200000000000000 0.250000000000000 0.300000000000000 0.350000000000000 0.400000000000000 0.450000000000000 0.500000000000000 0.550000000000000 0.600000000000000 0.650000000000000 0.700000000000000 0.750000000000000 0.800000000000000 0.850000000000000 0.900000000000000 0.950000000000000 1;3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 25];
%--------------------
over =pars(1);%false;%
npg = pars(2);
nv  = pars(3);
order = pars(4);
smooth =pars(5);
quadrature_master= pars(6);
Knots=quadrature(Knots,quadrature_master,over, order, smooth,npg);


function Knots=quadrature(Knots,quadrature_master,over, order, smooth,npg)
if quadrature_master
    addpath 'IGA-quadrature-master';
    for i=1:length(Knots)
        Knots(i).Ggauss=[];
        Knots(i).w=[];
        Knots(i).xi=[];
        [Knots(i)]=IGA_quadrature_master(Knots(i),over, order, smooth);
    end
else
    for i=1:length(Knots)
        Knots(i).xi=[];
        Knots(i).w=[];
        Knots(i).Ggauss=[];
        Knots(i).Ggauss(1:length(Knots(i).k(1,:))-1)=npg;
    end
end


