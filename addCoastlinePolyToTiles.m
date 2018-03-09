coastlineShape ='E:\basemap\world\noaa_bound\GSHHS_shp\f\GSHHS_f_L1_3413.shp';

tilefile = 'V:/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/PGC_Imagery_Mosaic_Tiles_Arctic_coast.mat'; %PGC/NGA Tile definition file
tiles=load(tilefile);

S = shaperead(coastlineShape);

xmin = nan(size(S));
xmax = xmin;
ymin = xmin;
ymax = xmin;

i=1;
for i=1:length(S)
    
    xmin(i) = S(i).BoundingBox(1);
    xmax(i) = S(i).BoundingBox(2);
    
    ymin(i) = S(i).BoundingBox(3);
    ymax(i) = S(i).BoundingBox(4);
    
end

vxtile = [tiles.x0,tiles.x0,tiles.x1,tiles.x1,tiles.x0];
vytile = [tiles.y0,tiles.y1,tiles.y1,tiles.y0,tiles.y0];

i=1;
coastline = cell(size(tiles.I));
for i=1:length(tiles.x0)
    surroundFlag = false;
    fprintf('tile %d of %d \n',i,length(tiles.x0));
    
    n = find(tiles.x0(i) <= xmax & tiles.x1(i) >= xmin &...
        tiles.y0(i) <= ymax & tiles.y1(i) > ymin);
    
    if ~isempty(n)
        
        j=1;
        for j=1:length(n)
            
            xv = S(n(j)).X;
            yv = S(n(j)).Y;
            
            %find poly breaks and remove nested polys
            
            nn = [0,find(isnan(xv))];
            if length(nn) > 2
                [~,dnn] = max(diff(nn));
                xv = xv(nn(dnn)+1:nn(dnn+1)-1);
                yv = yv(nn(dnn)+1:nn(dnn+1)-1);
            end
            
            % check to see if all points outside vertices tile
            
            intile = xv >= tiles.x0(i) & xv <= tiles.x1(i) &...
                yv >= tiles.y0(i) & yv <= tiles.y1(i);
            
            if ~any(intile)
                
                % is the tile surrouded?
                
                if  any([~any(xv <= tiles.x0(i) & yv <= tiles.y0(i));...
                        ~any(xv <= tiles.x0(i) & yv >= tiles.y1(i));...
                        ~any(xv >= tiles.x1(i) & yv >= tiles.y1(i));...
                        ~any(xv >= tiles.x1(i) & yv <= tiles.y0(i))])
                    
                    % any any conditions are null, this means that that at
                    % least quadant does not contain points = not
                    % surrounded = ocean
                    
                    continue;
                    
                else
                    
                    surroundFlag = true;
                    
                end
                
            end
            

            
            if surroundFlag
                
                coastline{i}{j} = [vxtile(i,:);vytile(i,:)];
                clear land
                break;
                
                
            else
                
                if ~exist('land','var')
                    % build grid
                    res = 40;
                    buff = 0;
                    
                    x = tiles.x0(i)-buff*res: res:tiles.x1(i)+buff*res;
                    y = tiles.y1(i)+buff*res:-res:tiles.y0(i)-buff*res;
                    y = y(:);
                    
                    land =  false(length(y),length(x));
                end
                
                nn =isfinite(xv) & isfinite(yv);
                land(roipoly(x,y,land,xv(nn),yv(nn)))=true;
                
                
                
            end
            
            
        end
        
        
        if exist('land','var')
            B = bwboundaries(land, 8, 'noholes'); % find data coverage boundaries
            clear land
            
            j=1;
            for j=1:length(B)
                
                xp=x(B{j}(:,2));
                yp=y(B{j}(:,1));
                
                coastline{i}{j} = [xp(:)';yp(:)']; % find data coverage boundaries
            end
            
        end
        
        
%          j=1;
%         for j=1:length(coastline{i})
%             
%             plot(coastline{i}{j}(1,:),coastline{i}{j}(2,:));
%             hold on
%             drawnow
%             
%         end
        
    end
end

tiles.coastline = coastline;

save(tilefile,'-struct','tiles');








