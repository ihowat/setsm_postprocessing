function undoBoundaryAdjust(fileNames,varargin)


if ~iscell(fileNames)
    fileNames = {fileNames};
end


%% dz0 apply loop - this could be run simultaneously/parallel
i=1;
for i =1:length(fileNames)
    fprintf('removing offset from %s ....',fileNames{i})
    
     m0=matfile(fileNames{i});
%      
%     if any(strcmp(fields(m0),'adjusted'))
%         if ~m0.adjusted 
%             fprintf('adjusted flag already false, skipping\n')
%             continue
%         end
%     end
%     
%     if ~any(strcmp(fields(m0),'dz0'))
%             fprintf('dz0 doesnt exist, skipping\n')
%             continue
%     end
    
   m0.Properties.Writable = true;
   m0.z = m0.z + m0.dzfit;
   m0.adjusted = false;
   fprintf('done\n')
end

        