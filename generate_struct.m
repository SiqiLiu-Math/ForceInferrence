function [Cdat,Fdat,Edat,Psi,L] = generate_struct(im)
%GENERATE_STRUCT Extract geometry and topology inf from image
%   Cdat - Cellular information
%   Fdat - Face (membrane) information
%   Edat - Edge information
%   Psi - Voronoi parameters (power, site, weight)
%   L - Label matrix of image



%% read image
info = imfinfo(im);
x = info(1).Height;
y = info(1).Width;
z = length(info);
L = zeros(x,y,z);
for ii=1:z
    L(:,:,ii) = imread(im,ii);
end
L = L(2:end-1,2:end-1,2:end-1);
dim = size(L);


%% generate Cdat
Cdat = regionprops(L,'PixelIdxList');
for ii=length(Cdat):-1:1
    Cdat(ii).Index = ii;
    if isempty(Cdat(ii).PixelIdxList)
        Cdat(ii) = [];
    end
end
cellnum = length(Cdat)-1;
Lt = (L>1);
Ld = imdilate(Lt,ones(3,3,3));
Cdat(1).PixelIdxList = find(Ld-Lt);
for ii=1:cellnum+1
    [x1,y1,z1] = ind2sub(dim,Cdat(ii).PixelIdxList);
    ll = length(Cdat(ii).PixelIdxList);
    x2 = x1 + ones(ll,1)*[0,0,0,0,1,1,1,1];
    y2 = y1 + ones(ll,1)*[0,0,1,1,0,0,1,1];
    z2 = z1 + ones(ll,1)*[0,1,0,1,0,1,0,1];
    id2 = sub2ind(dim+1,x2,y2,z2);
    Cdat(ii).Id = unique(id2);
    Cdat(ii).Center = [mean(x1),mean(y1),mean(z1)];
    Cdat(ii).Volume = ll;
    Cdat(ii).AdjacentCells = [];
    Cdat(ii).Faces = [];
    Cdat(ii).Edges = [];
    L(Cdat(ii).PixelIdxList) = ii;
end

%% generate Fdat
Fdat = struct('Cells',[],'Pixels',[],'Id',[],'Edges',[]);
kk = 0;
ad = zeros(cellnum+1);
for ii=1:cellnum
    for jj=ii+1:cellnum+1
        id = intersect(Cdat(ii).Id,Cdat(jj).Id,'sorted');
        if ~isempty(id)
            kk = kk+1;
            [x1,y1,z1] = ind2sub(dim+1,id);
            Fdat(kk).Cells = [ii,jj];
            Fdat(kk).Id = id;
            Fdat(kk).Index = [Cdat(ii).Index,Cdat(jj).Index];
            Fdat(kk).Pixels = [x1,y1,z1];
            Cdat(ii).AdjacentCells = [Cdat(ii).AdjacentCells, jj];
            Cdat(jj).AdjacentCells = [Cdat(jj).AdjacentCells, ii];
            Cdat(ii).Faces = [Cdat(ii).Faces, kk];
            Cdat(jj).Faces = [Cdat(jj).Faces, kk];
            ad(ii,jj) = kk;
        end
    end
end
facenum = kk;


%% generate Edat
Edat = struct('Id',[],'Faces',[],'Cells',[]);
kk = 0;
for f1 = 1:facenum
    c1 = Fdat(f1).Cells(1);
    c2 = Fdat(f1).Cells(2);
    adc = intersect(Cdat(c1).AdjacentCells,Cdat(c2).AdjacentCells);
    adc = adc(adc > max(c1,c2));
    for c3 = adc
        f2 = ad(c1,c3);
        f3 = ad(c2,c3);
        id = intersect(Fdat(f1).Id,Fdat(f2).Id,'sorted');
        if ~isempty(id)
            kk = kk+1;
            [x1,y1,z1] = ind2sub(dim+1,id);
            Edat(kk).Id=id;
            Edat(kk).Pixels = [x1,y1,z1];
            Edat(kk).Faces = [f1,f2,f3];
            Edat(kk).Cells = [c1,c2,c3];
            Fdat(f1).Edges = [Fdat(f1).Edges, kk];
            Fdat(f2).Edges = [Fdat(f2).Edges, kk];
            Fdat(f3).Edges = [Fdat(f3).Edges, kk];
            Cdat(c1).Edges = [Cdat(c1).Edges, kk];
            Cdat(c2).Edges = [Cdat(c2).Edges, kk];
            Cdat(c3).Edges = [Cdat(c3).Edges, kk];
        end
    end
end


%% Generate Psi
Psi = zeros(cellnum,5);
for ii=1:cellnum
    Psi(ii,1)   = random('Normal',2,0.1);
    Psi(ii,3:5) = Cdat(ii+1).Center;
    Psi(ii,2)   = random('Normal',40,1);
end
Cdat = rmfield(Cdat,{'PixelIdxList','Id'});

end
