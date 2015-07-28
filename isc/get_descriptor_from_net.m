function [dsc_raw,M] = get_descriptor_from_net(in_ray,in_ring,desc,ar)
Nrays = size(in_ray,1);
nscales =size(in_ring,1);
ndesc = size(desc,2);
dsc_raw = zeros(ndesc,nscales,Nrays);

count = 1;
I = [];
J = [];
V = [];

for sc =1:nscales,
    wt = find(in_ring(sc,:)>.001);
    arz = ar(wt)'.* in_ring(sc,wt);
    prds = repmat(arz,[Nrays,1]).*in_ray(:,wt);
    for ray = [1:Nrays],
        prd  = prds(ray,:);
        keep = find(prd>.01*max(prd));
        if isempty(keep),[~,keep] = max(prd);end
        nrm  = max(sum(prd(keep)),eps);
        wntr = wt(keep);
        wgtr = prd(keep);
        dsc_raw(:,sc,ray) = (wgtr*desc(wntr,:))/nrm;
        
        wntr = wntr';
        vals = wgtr ./ nrm;
        vals = vals';
        
        I = [I;count * ones(length(wntr),1)];
        J = [J;wntr];
        V = [V;double(vals)]; 
        count = count + 1;
    end
end

M = sparse(I,J,V,nscales * Nrays, length(ar));

