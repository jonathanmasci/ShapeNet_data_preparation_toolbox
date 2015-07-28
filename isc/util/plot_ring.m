function plot_ring(shape,rr);

vertex = [shape.X';shape.Y';shape.Z'];


face = shape.TRIV';

[v1,v2] = compute_levelset_mesh(vertex,face,shape.D,rr,[]);

X  =[v1(1,:)];
Y = [v1(2,:)];
Z = [v1(3,:)];
dsts  = dist2([X;Y;Z]',[X;Y;Z]');
dsts = dsts  + 1e6*eye(size(dsts));
[srt,idx] = sort(dsts,2);
i1 = idx(:,1);

for k=1:length(i1)
    h = plot3(X([k,i1(k)]),Y([k,i1(k)]),Z([k,i1(k)]),'linewidth',2,'color',[0,0,0]);
    nn  = i1;
    
    %% trailing nns
    nns = idx(k,:);
    st = 2;
    cont = 1;
    while cont,
        idnn = nns(st);
        if min(dsts(idnn,nns(1:st-1)))<dsts(k,idnn)
            st = st+1;
            if st==length(i1)
                cont = 0;
            end
        else
            cont = 0;
        end
    end
    i2 = idnn;
    hold on;
    h = plot3(X([k,i2]),Y([k,i2]),Z([k,i2]),'linewidth',2,'color',[0,0,0]);
end


    function [v1,v2] = compute_levelset_mesh(vertex,face,f,tau,options);
        
        % compute_levelset_mesh - compute level set curve
        %
        %   [v1,v2] = compute_levelset_mesh(vertex,face,f,tau,options);
        %
        %   v1(:,i),v2(:,i) is a segment of a levelt set of point f(x)=tau on the
        %   mesh.
        %
        %   Copyrigh (c) 2007 Gabriel Peyre
        
        
        if length(tau)>1
            v1 = [];
            v2 = [];
            for i=1:length(tau)
                [w1,w2] = compute_levelset_mesh(vertex,face,f,tau(i),options);
                v1 = [v1, w1];
                v2 = [v2, w2];
            end
            return;
        end
        
        f = f-tau;
        
        s = sum( f(face)>0 );
        I = find( s<3 & s>0 );
        m = length(I);
        
        % find intersection points
        t = f(face(:,I))./( f(face(:,I))-f(face([2 3 1],I)) );
        t(t<0 | t>1) = Inf;
        [tmp,J1] = min(t,[],1); t1 = t(J1+(0:m-1)*3);
        t(J1+(0:m-1)*3) = Inf;
        [tmp,J2] = min(t,[],1); t2 = t(J2+(0:m-1)*3);
        A1 = face(J1+(I-1)*3);
        A2 = face(mod(J1,3)+1+(I-1)*3);
        
        B1 = face(J2+(I-1)*3);
        B2 = face(mod(J2,3)+1+(I-1)*3);
        
        v1 = vertex(:,A1).*repmat(1-t1,[3 1]) + vertex(:,A2).*repmat(t1,[3 1]);
        v2 = vertex(:,B1).*repmat(1-t2,[3 1]) + vertex(:,B2).*repmat(t2,[3 1]);
