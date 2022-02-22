function undoTile(fileName)

m=matfile(fileName);
m.Properties.Writable = true;
s=whos(m);
sz= s(strcmp({s.name},'z')).size;
% undoBuffer

%undo dz0 addjustment 
if any(strcmpi(fields(m),'adjusted')) && any(strcmpi(fields(m),'dz0'))
    if m.adjusted == 1
        fprintf('undoing adjustment\n')
        dz0 = imresize(m.dz0,sz);
        m.z = m.z + dz0;
        m.adjusted=false;
    else
        fprintf('adjusted flag is zero, no change made\n')
    end
else
     fprintf('no adjusted and/or dz0 field, no change made\n')
end

%undo dzfit
if any(strcmpi(fields(m),'dzfit')) 
    %if m.fitted == 1
    if ~any(strcmpi(fields(m),'dzfitApplied')) || m.dzfitApplied == 1
        fprintf('undoing dzfit\n')
        dzfit = imresize(m.dzfit,sz);
        m.z = m.z + dzfit;
        m.dzfitApplied=false;
    else
        fprintf('dzfitApplied flag is zero, no change made\n')
    end
        
%     else
%         fprinft('fitted flag is zero, no change made\n')
%     end
else
     fprintf('no dzfit field, no change made\n')
end

    


    
        
