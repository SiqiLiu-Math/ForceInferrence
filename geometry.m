function [Cdat,Fdat,Edat] = geometry(Cdat,Fdat,Edat,Psi,Mu)
%GEOMETRY geometry and mechanical results
%   Cdat - Cellular information
%   Fdat - Face (membrane) information
%   Edat - Edge information
%   Psi - Voronoi parameters (power, site, weight)
%   Mu - Constant for zero mode

%% Cell parameters
Psi = [0,0,0,0,0;Psi];
for ii=1:length(Cdat)
    Cdat(ii).Power = Psi(ii,1);
    Cdat(ii).Weight = Psi(ii,2);
    Cdat(ii).Site = Psi(ii,3:5);
    Cdat(ii).Pressure  = Psi(ii,1)-Psi(ii,1)^2*Mu;
end

%% Membrane force and geometry
for ii=1:length(Fdat)
    c1 = Fdat(ii).Cells(1);
    c2 = Fdat(ii).Cells(2);
    q1 = Cdat(c1).Site;
    q2 = Cdat(c2).Site;
    p1 = Cdat(c1).Power;
    p2 = Cdat(c2).Power;
    t1 = Cdat(c1).Weight;
    t2 = Cdat(c2).Weight;
    px = Fdat(ii).Pixels;
    rho = (p1*q1-p2*q2)/(p1-p2);
    r = sqrt(p1*p2*(q1-q2)*(q1-q2)'/(p1-p2)^2+(p1*t1^2-p2*t2^2)/(p1-p2));
    dp = px-ones(size(px,1),1)*rho;
    er = mean((sqrt(sum(dp.^2,2))-r).^2);
    Fdat(ii).Radius = r;
    Fdat(ii).Centroid = rho;
    Fdat(ii).Error = er;
    Fdat(ii).Tension = abs(p1-p2)*r/2*(1-Mu*(p1+p2));
end

%% Line force and geometry
for ii=1:length(Edat)
    f1 = Edat(ii).Faces(1);
    f2 = Edat(ii).Faces(2);
    o1 = Fdat(f1).Centroid;
    o2 = Fdat(f2).Centroid;
    r1 = Fdat(f1).Radius;
    r2 = Fdat(f2).Radius;
    c1 = Edat(ii).Cells(1);
    c2 = Edat(ii).Cells(2);
    c3 = Edat(ii).Cells(3);
    q1 = Cdat(c1).Site;
    q2 = Cdat(c2).Site;
    q3 = Cdat(c3).Site;
    p1 = Cdat(c1).Power;
    p2 = Cdat(c2).Power;
    p3 = Cdat(c3).Power;
    d = (r1^2-r2^2)/norm(o1-o2)+norm(o1-o2);
    rho = (o2-o1)/norm(o1-o2)*d/2+o1;
    r = sqrt(r1^2-(d/2)^2);
    n = (o1-o2)/norm(o1-o2);
    px = Edat(ii).Pixels;
    rr = mean(px,1)-rho;
    r0 = cross(rr,n)/norm(cross(rr,n))*r+rho;
    qq1 = (q1-r0)*p1+r0;
    qq2 = (q2-r0)*p2+r0;
    qq3 = (q3-r0)*p3+r0;
    s = cross(qq1-qq2,qq1-qq3);
    dp = px-ones(size(px,1),1)*rho;
    er = mean((sqrt(sum(dp.^2,2))-r).^2);
    Edat(ii).Radius = r;
    Edat(ii).Center = rho;
    Edat(ii).Axis = n;
    Edat(ii).Tension = norm(s)/2*Mu;
    Edat(ii).Error = er;
end

%% Stress tensor of cell
apicalface = length(Cdat(1).Faces);
apicaledge = length(Cdat(1).Edges);
for ii=2:length(Cdat)
    v = Cdat(ii).Volume;
    p = Cdat(ii).Pressure;
    sigma = -p*eye(3);
    for ff=Cdat(ii).Faces
        rho = Fdat(ff).Centroid;
        ten = Fdat(ff).Tension;
        px  = Fdat(ff).Pixels;
        for jj=1:size(px,1)
            n = (px(jj,:)-rho)/norm(px(jj,:)-rho);
            sigma = sigma+ten/2/v*(eye(3)-n'*n);
            if ff<apicalface+1
                sigma = sigma+ten/2/v*(eye(3)-n'*n);
            end
        end
    end
    for ee=Cdat(ii).Edges
        rho = Edat(ee).Center;
        ten = Edat(ee).Tension;
        px  = Edat(ee).Pixels;
        eta = Edat(ee).Axis;
        for jj=1:size(px,1)
            n = (px(jj,:)-rho)/norm(px(jj,:)-rho);
            t = cross(eta,n);
            sigma = sigma+ten/3/v*(t'*t);
            if ee<apicaledge+1
                sigma = sigma+ten/6/v*(t'*t);
            end
        end
    end
    Cdat(ii).StressTensor = sigma;
end


end

