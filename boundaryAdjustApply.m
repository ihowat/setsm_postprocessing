function boundaryAdjustApply(fileName)

m=matfile(fileName);

if any(strcmp(fields(m),'adjusted'))
    if m.adjusted==1
        fprintf('adjusted flag true, skipping\n')
        return
    end
end

if ~any(strcmp(fields(m),'dz0'))
    fprintf('dz0 doesnt exist, skipping\n')
    return
end

sz=  size(m,'z');

m.Properties.Writable = true;
dz0 = imresize(m.dz0,sz);
m.z = m.z - dz0;
m.adjusted = true;

