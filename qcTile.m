function qcTile(tileFile,meta)
% qcTile: quality control REMA/ArcticDEM tile mosaic
% qcTile(tileFile,meta) where tileFile contains the file path/name of the
% tile and meta is the meta data structure. The tileFile must end with
% *dem.mat. If a file named *dem_shade.tif exists, that will be used for
% the qc display and, if not, a hillshade will be created. Hillshades must
% exist for the individual strips, however.

% open the tileFile into a matfile structure
m=matfile(tileFile);


% check if shade tif exists and read if it does
shadeFile = strrep(tileFile,'dem.mat','dem_shade.tif');
if strcmp(tileFile,shadeFile); shadeFile =[]; end
if exist(shadeFile,'file')
    disp('reading geotiff')
    h = readGeotiff(shadeFile);
    x = h.x;
    y = h.y;
    h = h.z;
    
else
    % no shade so create one
    h=hillshade(m.z,m.x,m.y);
    x = m.x;
    y = m.y;
end

% plot the hillshade in the current figure window
imagesc(x,y,uint8(h));
axis xy equal;
colormap gray
set(gca,'clim', [150 200])

hold on

deleteFlag = false;
k=1;
while k
    
    s=input('Draw ROI? (y/n)\n','s');
    
    if strcmpi(s,'y')
        
        deleteFlag = true;
        
        if ~exist('f','var')
            %% match strip files to database
            f=m.f; % extract image filename vector from from matfile
            f(isnan(m.rmse))=[]; % remove unused data
            [~,IA]= intersect(meta.f,f); % find index of tile files in meta
            meta = structfun(@(x) ( x(IA,:) ), meta, 'UniformOutput', false); %crop meta
            meta = addQC2Meta(meta); % add existing qc data
        end
        
        [~,xv,yv] =  roipoly; % draw ROI
        plot(xv,yv,'b','linewidth',2); % plot ROI
        
        % find strips overlapping ROI
        n = find(meta.xmax > min(xv) &  meta.xmin < max(xv) & ...
            meta.ymax > min(yv) &  meta.ymin < max(yv));
        
        in=zeros(size(n));
        for i=1:length(in)
            in(i) = any(inpolygon(meta.x{n(i)},meta.y{n(i)},xv,yv)) | ...
                any(inpolygon(xv,yv,meta.x{n(i)},meta.y{n(i)}));
        end
        n=n(logical(in));
        
        % select the whichever registration has the better sigma_bias (all or 1 yr)
        meta.sigma = nanmin([meta.sigma_all(:)';meta.sigma_1yr(:)'])';
        meta.sigma(meta.sigma > 1) = NaN;

        meta.avg_rmse(meta.avg_rmse == 0) = NaN;
        
        f2=figure;
        i=1;
        for i=1:length(n)
     
            % load and plot strip hillshade
            fileName=strrep(meta.f{n(i)},'meta.txt','dem_browse.tif');
            fprintf('loading %d of %d: %s\n',i,length(n),fileName);
            I=readGeotiff(fileName);
            imagesc(I.x,I.y,I.z,'alphadata',single(I.z ~= 0))
            set(gca,'color','r')
            axis xy  equal tight;
            colormap gray;
            set(gcf,'units','normalized');
            set(gcf,'position',[0.01,0.01,.35,.9])
            hold on
            % plot ROI
            plot(xv,yv,'b','linewidth',2)
            set(gca,'xlim',[min(I.x)-500 max(I.x)+500],'ylim',[min(I.y)-500 max(I.y)+500])
            
            % load qc data
            qc=load([fileparts(fileName),'/qc.mat']);
            [~,IA]=intersect(qc.fileNames, fileName);
            
            % existing qc is 3, plot the mask polys
            j=1;
            if qc.flag(IA) == 3
                hl = [];
                for j=1:length(qc.x{IA})
                   if ~isempty(qc.x{IA}{j})
                    hl =  [hl,plot(qc.x{IA}{j},qc.y{IA}{j},'g','linewidth',2)];
                   end
                end
                
            end
            
            %% specify flag
            j=1;
            while j
                
                fprintf('gcp sigma=%.2f, mean coreg RMSE=%.2f, max coreg RMSE=%.2f, current qc flag: %d\n',...
                    meta.sigma(n(i)),meta.avg_rmse(n(i)), meta.max_rmse(n(i)), qc.flag(IA));
                
                try
                    flag=input('Enter quality flag: 1=good, 2=partial, 3=manual edit, 4=poor, 5=bad or just enter to leave as is/error\n');
                catch
                    
                    fprintf('%d not recogized, try again\n');
                end
                
                if isempty(flag); break; end
                
                try
                    if isnumeric(flag)
                        if flag == 1 || flag == 2 || flag == 3 || flag == 4 || flag == 5
                            break;
                        end
                    end
                catch
                    fprintf('%d not recogized, try again\n',flag);
                end
                
                if iscell(flag); flag=flag{1}; end
                
                fprintf('%d not recogized, try again\n',flag);
                
            end
            
            if isempty(flag); clf; continue; end
            
            %% if manual selected 
            if flag ==3
                  
                
                j=1;
                if qc.flag(IA) == 3
    
                    s=input('keep existing mask polygons? (y/n)\n','s');
                    
                    if strcmpi(s,'n')
                        
                        qc.x{IA}=cell(1); qc.y{IA}=cell(1);
                        delete(hl);
                    else
                        
                        j= length(qc.x{IA})+1;
                    end
                    
                end

                fprintf('entering manual edit mode\n')
                
                 while j
                    
                     fprintf('drawing mask region %d\n',j);
                    
                     
                     figure(f2);
                     [~,qc.x{IA}{j},qc.y{IA}{j}] =  roipoly;

                    plot(qc.x{IA}{j},qc.y{IA}{j},'g','linewidth',2)
                    
                    
                    while j
                        s=input('continue editing this image? (y/n)\n','s');
                        
                        if ~strcmpi(s,'y') && ~strcmpi(s,'n')
                            fprintf('%s not recogized, try again\n',s);
                        else
                            break
                        end
                    end
                    
                    if strcmpi(s,'n'); break; end
                    
                    j=j+1;
                 end      
            end
            
            
            qc.flag(IA)=flag;
            
            if flag == 1 || flag == 2
                qc.x{IA}=cell(1); qc.y{IA}=cell(1);
                
            end
            
            if flag == 4 || flag == 5
                qc.x{IA}=cell(1); qc.y{IA}=cell(1);
                qc.x{IA}{1} = [min(I.x);min(I.x);max(I.x);max(I.x);min(I.x)];
                qc.y{IA}{1} = [min(I.y);max(I.y);max(I.y);min(I.y);min(I.y)];
            end
            
            
            
            clf
            
            save([fileparts(fileName),'/qc.mat'],'-struct','qc');
            
            
        end
        
        close(f2)
        
    elseif strcmpi(s,'c')
         clim(1) =input('enter min c value\n');
         clim(2) = input('enter max c value\n');
        
         set(gca,'clim',clim);
              
    else
        break
    end
    
end

if deleteFlag

    s=input('Delete edited tile file (y/n)\n','s');
    
    
    if strcmpi(s,'y')
        
         [filePath,fileName] = fileparts(tileFile);
         delete([filePath,'/',fileName(1:5),'*']) 
    end
end

clf
