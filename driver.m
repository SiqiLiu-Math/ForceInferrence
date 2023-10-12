clear

path1 = 'Astec/';
path2 = 'results/';

%% read images and generate data
image_list = dir([path1,'*.tif']);
for nn = 1:length(image_list)
    [Cdat,Fdat,Edat,Psi,L] = generate_struct([path1,image_list(nn).name]);
    save([path2,image_list(nn).name(1:end-4)],'Cdat','Fdat','Edat','Psi','L');
end

%% infer Voronoi parameters
result_list = dir([path2,'*.mat']);
lmt = 5;
for nn=1:length(result_list)
    load([path2,result_list(nn).name]);
    Psi = force_inference(Cdat,Fdat,Psi,lmt);
    save([path2,result_list(nn).name],'Cdat','Fdat','Edat','Psi','L');
end

lmt = 2;   % improve the results
for nn=1:length(result_list)
    load([path2,result_list(nn).name]);
    Psi = force_inference(Cdat,Fdat,Psi,lmt);
    save([path2,result_list(nn).name],'Cdat','Fdat','Edat','Psi','L');
end

%% compute mechanical results
for nn=1:length(result_list)
    load([path2,result_list(nn).name]);
    [Cdat,Fdat,Edat] = geometry(Cdat,Fdat,Edat,Psi,0);
    save([path2,result_list(nn).name],'Cdat','Fdat','Edat','Psi','L');
end

