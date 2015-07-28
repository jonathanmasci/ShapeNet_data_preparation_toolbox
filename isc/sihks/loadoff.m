function shape = loadoff(filename)

shape = [];

f = fopen(filename, 'rt');

n = '';
while isempty(n)
    fgetl(f);
    n = sscanf(fgetl(f), '%d %d %d');
end
    
nv = n(1);
nt = n(2);
data = fscanf(f, '%f');

if(length(data) == nv*3 + nt*3)
    numsInTri = 3;
else 
    if(length(data) == nv*3 + nt*4)
        numsInTri = 4;
    else
        error('file format not supported');
    end
end


shape.TRIV = reshape(data(end-numsInTri*nt+1:end), [4 nt])';
if(numsInTri ==4)
    shape.TRIV = shape.TRIV(:,2:4);
end
if(shape.TRIV(1) == 0)
    shape.TRIV = shape.TRIV + 1;
end


data = data(1:end-numsInTri*nt);

data = reshape(data, [length(data)/nv nv]);

shape.X = data(1,:)';
shape.Y = data(2,:)';
shape.Z = data(3,:)';

fclose(f);

