function meta=readSceneMeta(metaFile)
%READMETA read SETSM output metafile for an individual DEM
%
% meta=readSceneMeta(metaFile) reads each line of the metaFile into the
% structure meta. The field neames are the same as in the metaFile but with
% underscores for spaces and all lowercase.
%
% Ian Howat, ihowat@gmail.com
% 24-Jul-2017 14:57:49

% read each line into a cellstr
c=textread(metaFile,'%s','delimiter','\n');

% initialize output strucure
meta=struct();

% loop through each line in c and import
for i=1:length(c)
    
    % initialize name/value cell
    s=cell(1,2);
    
    % find location of 1st =. strsplit doest let you split by only 1st
    % appearance.
    n = find(c{i}=='=');
    
    % no equal sign means no name/value pair so skip
    if isempty(n); continue; end
    
    % parse line by n
    s{1}= c{i}(1:n-1);
    s{2}= c{i}(n+1:end);
    
    % remove trailing spaces in field name
    s{1}=deblank(s{1});
    
    % make spaces underscores
    s{1}=lower(strrep(s{1},' ','_'));
    
    % second occurence of this field, add "image_1_" to existing and
    % "image_2_" to this one.
    if isfield(meta,s{1})
        eval(['meta.image_1_',s{1},'=meta.',s{1}]);
        meta=rmfield(meta,s{1});
        s{1}=['image_2_',s{1}];
    end
    
    % preserve single quotes (needed for proj string)
    s{2} = strrep(s{2},'''','''''');
    
    % remove trailing spaces in value
    s{2}=deblank(s{2});
    
    % test if character value and add single quotes to import
    if (isempty(str2num(s{2})) || any(isnan(str2num(s{2})))) && ~strcmpi(s{2},'nan')
        eval(['meta.',s{1},'=''',s{2},''';']);
    else
        % otherwise import as numeric scalar or 1 by n vector
        eval(['meta.',s{1},'=[',s{2},'];']);
    end
    
end


