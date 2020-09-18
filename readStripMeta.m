function meta=readStripMeta(metaFileName,varargin)
% readStripMeta: reads strip metadata files into a a data structure
%
% meta=readStripMeta(metaFileName) returns the strip and scene metadata in
% a structure with fields:
%     {'creation_date'         } 
%     {'strip_creation_date'   }
%     {'strip_projection_proj4'}
%     {'x'                     }
%     {'y'                     }
%     {'scene_alignment'       }
%     {'scene_meta'            } 
% where 'scene_meta' is a substructure containing the scene metadata for
% each scene comprising the strip.
%
% Ian Howat (ihowat@gmail.com)
% version 1.0: '01-Aug-2019 15:59:05'

noSceneMetaFlag=false;
if ~isempty(varargin)
	noSceneMetaFlag=true;
end

% initialize output
meta.fileName=metaFileName;

% read each line of meta file into a cell array (textscan doesnt work)
C=textread(metaFileName,'%s','delimiter','\n');

% iterate thrhough C and seperate out name/value pairs into cells.
cnt=1;
name=cell(1);
val=cell(1);
for it=1:length(C)
    c=C{it};
    n=find((c == '=') | (c == ':') );
    
    if isempty(n); continue; end
    
    name{cnt}=deblank(lower(c(1:n(1)-1)));
    val{cnt}=deblank(lower(c(n(1)+1:end)));
    
    % convert numeric vals to nums. warning off due to character length
    % warnings.
    warning off
    numVal= str2num(val{cnt});
    warning on
    if ~isempty(numVal)
        val{cnt}=numVal;
    end
    
    cnt=cnt+1;
end

% remove offending chars from names
name = strrep(name,' ','_');
name = strrep(name,'(','');
name = strrep(name,')','');
name = strrep(name,',','');
name = strrep(name,'-','_');



% Convert date strings to datenums
n=find(strcmp(name,'creation_date') | strcmp(name,'strip_creation_date'));
val(n)=strrep(val(n),'mon',''); % some weird formatting was used on scenes
val(n)=strrep(val(n),'tue','');
val(n)=strrep(val(n),'wed','');
val(n)=strrep(val(n),'thu','');
val(n)=strrep(val(n),'fri','');
val(n)=strrep(val(n),'sat','');
val(n)=strrep(val(n),'sun','');
val(n)=cellfun(@datenum,val(n),'uniformoutput',0);

% locate and break out scene metadata
scene_num=1;
first_scene_line=[];
while  scene_num
    
    line=find(strcmp(name,['scene_',num2str(scene_num),'_name']));
   
    if isempty(line); break; end
    
    %meta for this scene begins on this line
    first_scene_line(scene_num)=line; 
    
    scene_num= scene_num+1;
end

% make all the scene name fields the same
name(first_scene_line)={'scene_name'};

% number of scenes
Nscenes=length(first_scene_line);

if ~noSceneMetaFlag

	% seperate out scene lines per scene
	scene_meta_name=cell(1,Nscenes);
	scene_meta_val=cell(1,Nscenes);
	first_scene_line=[first_scene_line,length(name)];
	for it=1:Nscenes
    		scene_meta_name{it}=name(first_scene_line(it):first_scene_line(it+1)-1);
    		scene_meta_val{it}=val(first_scene_line(it):first_scene_line(it+1)-1);    
	end

end

% remove scene meta from main cells
name(first_scene_line(1):end)=[];
val(first_scene_line(1):end)=[];

% reformat strip data into structure
for it=1:length(name)
    if isstr(val{it})
        val{it}=strrep(val{it},'''','');
        val{it}=strrep(val{it},'(','');
        val{it}=strrep(val{it},')','');
        eval(['meta.',name{it},'=''',val{it},''';'])
    else
        eval(['meta.',name{it},'=[',num2str(val{it}),'];'])
    end
end

% calculate and save strip area
meta.A = polyarea(meta.x, meta.y);

% get scene alignment stats
r=find(contains(C,'Mosaicking Alignment Statistics'));
r=r+2;
A=[];
while r
    if isempty(C{r}); break; end
    A=[A;sscanf(C{r},'%*s %f %f %f %f %f %f %f')'];
    r=r+1;
end
meta.scene_alignment.rmse=A(:,1);
meta.scene_alignment.dz=A(:,2);
meta.scene_alignment.dx=A(:,3);
meta.scene_alignment.dy=A(:,4);
meta.scene_alignment.dze=A(:,5);
meta.scene_alignment.dxe=A(:,6);
meta.scene_alignment.dye=A(:,7);

if ~noSceneMetaFlag
	for n=1:Nscenes
    		for it=1:length(scene_meta_name{n})
        		if  isstr(scene_meta_val{n}{it})
            			scene_meta_val{n}{it}=strrep(scene_meta_val{n}{it},'''','');
            			eval(['meta.scene_meta(',num2str(n),').',scene_meta_name{n}{it},'=''',scene_meta_val{n}{it},''';'])
        		else
            			eval(['meta.scene_meta(',num2str(n),').',scene_meta_name{n}{it},'=[',num2str(scene_meta_val{n}{it}),'];'])
        		end
    		end
	end
end
