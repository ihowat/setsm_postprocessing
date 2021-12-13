function redoBoundaryAdjust(fileName,tileNeighborIndexFile)

load(tileNeighborIndexFile,'fileNames','nN')

i = find(strcmp(fileNames,fileName));
if isempty(i)
    error('filename not found in tileNeighborIndexFile fileNames')
end

neighborFiles = cell(8,1);
neighborFiles(~isnan(nN(i,:))) = fileNames(nN(i,~isnan(nN(i,:))));

j=1;
for j=1:length(neighborFiles)
    
    if isempty(neighborFiles{j})
        continue
    end

    m=matfile(neighborFiles{j});
    m.Properties.Writable = true;
    sz=size(m,'z');
    
    % NEED TO ADD undoBuffer
    
    if any(strcmpi(fields(m),'adjusted')) && any(strcmpi(fields(m),'dz0'))
        if m.adjusted == 1
            fprintf('undoing adjustment for %s\n',neighborFiles{j})
            dz0 = imresize(m.dz0,sz);
            m.z = m.z + dz0;
            m.adjusted=false;
        else
            fprintf('adjusted flag is zero, no change made\n')
        end
    else
        fprintf('no adjusted and/or dz0 field, no change made\n')
    end
end


boundaryAdjustCalc(fileName,neighborFiles)

j=1;
for j=1:length(neighborFiles)
    
    if isempty(neighborFiles{j})
        continue
    end
    
    k = find(strcmp(fileNames,neighborFiles{j}));
    if isempty(k)
        error('%s not found in tileNeighborIndexFile fileNames',neighborFiles{j})
    end

    neighborFiles1 = cell(8,1);
    neighborFiles1(~isnan(nN(k,:))) = fileNames(nN(k,~isnan(nN(k,:))));
    
    boundaryAdjustCalc(neighborFiles{j},neighborFiles1);
end
    
    
    