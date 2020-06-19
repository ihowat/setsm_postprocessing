
%coastlineShape = '/mnt/pgc/data/projects/arcticdem/coastline/GSHHS_f_L1_3413.shp';
%coastlineShape = '/mnt/pgc/data/projects/arcticdem/coastline/GSHHS_f_L1_GIMPgl_3413_clipArcDEMbuff100km.shp';
coastlineShape = 'E:\scratch\data\coastlines\GSHHS_f_L1_GIMPgl_3413_clipArcDEMbuff100km.shp';
%coastlineShape = 'E:\scratch\data\coastlines\GSHHS_f_L1_3413.shp';

%tilefile = '/mnt/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file
tilefile = 'V:/pgc/data/scratch/claire/repos/setsm_postprocessing_pgc/PGC_Imagery_Mosaic_Tiles_Arctic.mat'; %PGC/NGA Tile definition file

if exist('land','var')
    clear land
end

tiles=load(tilefile);

coastline = cell(size(tiles.I));
if isfield(tiles ,'coastline')
    coastline = tiles.coastline;
else
    coastline = cell(size(tiles.I));
end

%fix_tile_list = {};
fix_tile_list = {'58_23', '58_24', '59_22', '59_23', '59_24', '60_22', '60_23', '60_24'};

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
last_print_len=0;
for i=1:length(tiles.x0)
    surroundFlag = false;

    fprintf(repmat('\b', 1, last_print_len));
    last_print_len = fprintf('tile %d of %d: %s ',i,length(tiles.x0), tiles.I{i});
    if ~isempty(fix_tile_list) && ~any(strcmp(tiles.I{i}, fix_tile_list))
        continue
    else
        ;
%        coastline{i} = {};
    end

%    % find all polys whose bounding boxes intersect this tile
%    n = find(tiles.x0(i) <= xmax & tiles.x1(i) >= xmin &...
%        tiles.y0(i) <= ymax & tiles.y1(i) >= ymin);
    
%    if ~isempty(n)
%
%        j=1;
%        % for each poly whose bounding box intersects this tile
%        for j=1:length(n)
%
%            xv = S(n(j)).X;
%            yv = S(n(j)).Y;
%
%            %find poly breaks and remove nested polys
%
%            nn = [0,find(isnan(xv))];
%            if length(nn) > 2
%                [~,dnn] = max(diff(nn));
%                xv = xv(nn(dnn)+1:nn(dnn+1)-1);
%                yv = yv(nn(dnn)+1:nn(dnn+1)-1);
%            end
%
%            % check to see if all points outside vertices tile
%
%            % find all points in poly that fall inside this tile
%            intile = xv >= tiles.x0(i) & xv <= tiles.x1(i) &...
%                yv >= tiles.y0(i) & yv <= tiles.y1(i);
%
%            % if no points in poly fall inside this tile...
%            if ~any(intile)
%
%                % is the tile surrouded?
%
%                if  any([~any(xv <= tiles.x0(i) & yv <= tiles.y0(i));...
%                        ~any(xv <= tiles.x0(i) & yv >= tiles.y1(i));...
%                        ~any(xv >= tiles.x1(i) & yv >= tiles.y1(i));...
%                        ~any(xv >= tiles.x1(i) & yv <= tiles.y0(i))])
%
%                    % any any conditions are null, this means that that at
%                    % least quadant does not contain points = not
%                    % surrounded = ocean
%
%                    continue;
%
%                else
%
%                    surroundFlag = true;
%
%                end
%
%            end
%
%
%
%            if surroundFlag
%
%                coastline{i}{j} = [vxtile(i,:);vytile(i,:)];
%                clear land
%                break;
%
%
%            else
%
%                if ~exist('land','var')
%                    % build grid
%                    res = 40;
%                    buff = 0;
%
%                    x = tiles.x0(i)-buff*res: res:tiles.x1(i)+buff*res;
%                    y = tiles.y1(i)+buff*res:-res:tiles.y0(i)-buff*res;
%                    y = y(:);
%
%                    land =  false(length(y),length(x));
%                end
%
%                nn =isfinite(xv) & isfinite(yv);
%                land(roipoly(x,y,land,xv(nn),yv(nn)))=true;
%
%
%
%            end
%
%
%        end
%
%
%        if exist('land','var')
%            B = bwboundaries(land, 8, 'noholes'); % find data coverage boundaries
%            clear land
%
%            j=1;
%            for j=1:length(B)
%
%                xp=x(B{j}(:,2));
%                yp=y(B{j}(:,1));
%
%                coastline{i}{j} = [xp(:)';yp(:)']; % find data coverage boundaries
%            end
%
%        end
        
        
    j=1;
    for j=1:length(coastline{i})

        plot(coastline{i}{j}(1,:),coastline{i}{j}(2,:));
        hold on
        drawnow

    end
        
%    end
end

%tiles.coastline = coastline;
%
%save(tilefile,'-struct','tiles');








