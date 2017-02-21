function [meta] = stripSearch(meta,tilex0,tilex1,tiley0,tiley1)
% stripSearch filters strip metastruct to strips overlapping tile bounds
%
% [meta] = stripSearch(meta,tilex0,tilex1,tiley0,tiley1) where meta is the
%   arcticDEM meta structure that contains footprint vertices in .x and .y
%   and ranges in minx, maxx, miny , maxy. Non-overlapping strips are
%   removed from meta. NOTE:all meta fields must be record-per-row or else
%   they will be screwed up by the (n,:) removal.
%
% Ian Howat, ihowat@gmail.com, Ohio State

% make tile boundary polygon
tilevx = [tilex0;tilex0;tilex1;tilex1;tilex0];
tilevy = [tiley0;tiley1;tiley1;tiley0;tiley0];


% quick search: find strips within range of this tile. This does not
% account for background area of around strips but just pairs them down to
% speed the poly intersection loop
n = meta.xmax > tilex0 & meta.xmin < tilex1 & ...
    meta.ymax > tiley0 & meta.ymin < tiley1;

% if no overlap, return
if ~any(n); fprintf('no strip overlap\n'); meta=[]; return; end

% par down database structure to possible overlapping tiles
meta = structfun(@(x) ( x(n,:) ), meta, 'UniformOutput', false);

% search for all strip footprints overlapping this tile
n=false(size(meta.f));
for i=1:length(n)
    n(i) = any(inpolygon(meta.x{i},meta.y{i},tilevx,tilevy)) | ...
        any(inpolygon(tilevx,tilevy,meta.x{i},meta.y{i}));
end

% if no overlap, return
if ~any(n); fprintf('no strip overlap'); meta=[]; return; end

% remove files with no overlap
meta = structfun(@(x) ( x(n,:) ), meta, 'UniformOutput', false);

