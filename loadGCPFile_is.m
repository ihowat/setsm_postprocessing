function [gcp]=loadGCPFile_is(gcpfile)
% loadGCPFile_is load icesat csv file into structure

[~,name,ext]=fileparts(gcpfile);
fprintf('loading gcp file: %s\n',gcpfile)
switch ext
    case '.csv'
        fid=fopen(gcpfile);
        %-2171380.33353,868114.550919,469.9240,2003-10-14 18:58:57.794,0.0840
        D=textscan(fid,'%f32%f32%f32%u32-%u8-%u8 %u8:%u8:%f32%f32','delimiter',',');
        fclose(fid);
        gcp.x=D{1};
        gcp.y=D{2};
        gcp.z=D{3};
        gcp.t=datenum([double(D{4}) double(D{5}) double(D{6}) double(D{7}) double(D{8}) double(D{9})]);
        gcp.i_gmc=D{10};
        clear D;
        
        % remove missing data lines
        n = isnan(gcp.x) | isnan(gcp.y) | isnan(gcp.z);
        gcp = structfun(@(x) ( x(~n) ), gcp, 'UniformOutput', false);
        clear n
        
        % find and remove duplicates
        [~,n] = unique([gcp.x,gcp.y],'rows');
        gcp = structfun(@(x) ( x(n) ), gcp, 'UniformOutput', false);
        clear n
        
        gcp.dataset=name;
        
    otherwise
        error('gcp file must be .csv or .mat')
        
end
