function Psi = force_inference(Cdat,Fdat,Psi,lmt)
%FORCE_INFERENCE Infer Voronoi parameters by minimizing error function
%   Cdat - Cellular information
%   Fdat - Face information
%   Psi - Voronoi parameters (power, weight, site)
%   lmt - Limitation of minimizing

cellnum = length(Cdat)-1;

check_real = isreal(error_function(Psi,Fdat));
while ~check_real
    Psi = Psi.*random('Normal',1,0.001,cellnum,5);
    check_real = isreal(error_function(Psi,Fdat));
end
pp = Psi(:,1);
Psi(:,1) = pp/mean(pp);

%% minimize error function
ef = @(psi)error_function(psi,Fdat);
options1 = optimoptions(@fmincon,'Display','iter','MaxFunctionEvaluations',2e5,'ObjectiveLimit',lmt);
Aeq = [ones(1,cellnum),zeros(1,cellnum*4)];
beq = cellnum;
lb = zeros(5*cellnum,1);
ub = [5*ones(cellnum,1);100*ones(cellnum,1);300*ones(cellnum*3,1)];
Psi = fmincon(ef,Psi,[],[],Aeq,beq,lb,ub,[],options1);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%
function error = error_function(psi,Fdat)

psi = [0,0,0,0,0;psi];
error = 0;
pxn = 0;
for ii=1:length(Fdat)
    px = Fdat(ii).Pixels;
    c1 = Fdat(ii).Cells(1);
    c2 = Fdat(ii).Cells(2);
    p1 = psi(c1,1);
    p2 = psi(c2,1);
    q1 = psi(c1,3:5);
    q2 = psi(c2,3:5);
    t1 = psi(c1,2);
    t2 = psi(c2,2);
    rho = (p1*q1-p2*q2)/(p1-p2);
    r = sqrt(p1*p2*(q1-q2)*(q1-q2)'/(p1-p2)^2+(p1*t1^2-p2*t2^2)/(p1-p2));
    ww = px-ones(size(px,1),1)*rho;
    dd = sqrt(sum(ww.^2,2));
    error = error+sum((dd-r).^2);
    pxn = pxn+size(px,1);
end

error = error/pxn;

end


