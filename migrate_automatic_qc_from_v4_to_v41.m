
qc_dir_in  = '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tile_qc/tile_qc_v13e';
qc_dir_out = '/mnt/pgc/data/elev/dem/setsm/REMA/mosaic/v2/tile_qc/tile_qc_v13e_strips_v4.1';
remerge_info_matfile = '/mnt/pgc/data/projects/earthdem/strip_databases/REMAdatabase4_2m_v4.1_20220511_remerge_info.mat';


fprintf('Loading remerge info matfile: %s\n', remerge_info_matfile);
remerge_struct = load(remerge_info_matfile);


qc_file_listing = dir([qc_dir_in,'/*.mat']);
qc_files_cellarr = fullfile({qc_file_listing.folder}.', {qc_file_listing.name}.');

for qc_file_idx = 1:length(qc_files_cellarr)
    qc_file = qc_files_cellarr{qc_file_idx};
    fprintf('Loading qc file: %s\n', qc_file);
    qc = load(qc_file);

    [~,qc_fname,qc_ext] = fileparts(qc_file);
    qc_out_file = [qc_dir_out,'/',qc_fname,qc_ext];
    qc_out = struct();
    qc_out.stripID = {};
    qc_out.seg = double([]);
    qc_out.flag = uint8([]);
    qc_out.x = {};
    qc_out.y = {};

    stripID_cellarr = unique(qc.stripID);
    for stripID_idx = 1:length(stripID_cellarr)
        stripID = stripID_cellarr{stripID_idx};

        qc_this_stripID_idx_arr = find(strcmp(qc.stripID, stripID));
        qc_this_stripID_seg_arr = qc.seg(qc_this_stripID_idx_arr);

        stripID_in_remerge_info_idx = find(strcmp(remerge_struct.stripID, stripID));
        if length(stripID_in_remerge_info_idx) == 0
            fprintf('WARNING: stripID from qc file not found in remerge info matfile, skipping; %s, stripID: %s\n', qc_file, stripID);
            continue;
        elseif length(stripID_in_remerge_info_idx) > 1
            error('More than 1 match for stripID from qc file in remerge info matfile; %s, stripID: %s', qc_file, stripID);
        end

        remerge_info_file = remerge_struct.remerge_info_file{stripID_in_remerge_info_idx};

        for i = 1:length(qc_this_stripID_idx_arr)
            qc_idx = qc_this_stripID_idx_arr(i);

            qc_seg = qc.seg(qc_idx);
            qc_flag = qc.flag(qc_idx);
            qc_x = qc.x{qc_idx};
            qc_y = qc.y{qc_idx};

            if ~isempty(qc_x)
                if qc_flag == 3 || qc_flag == 5
                    ;
                else
                    error('qc exists with flag that is not 3 or 5; %s, stripID: %s, seg: %d, flag: %d', qc_file, stripID, qc_seg, qc_flag);
                end
            end
        end

%        if ~isfile(remerge_info_file)
%            fprintf('WARNING: remerge.info file does not exist for stripID from qc file, skipping; %s, stripID: %s, remerge.info file: %s\n', qc_file, stripID, remerge_info_file);
%            continue;
%        end
%        [v4_seg_arr, v41_seg_arr] = readRemergeInfoFile(remerge_info_file);

        if isfile(remerge_info_file)
            [v4_seg_arr, v41_seg_arr] = readRemergeInfoFile(remerge_info_file);
        else
            [strip_dir,~,~] = fileparts(remerge_info_file);
            v4_nseg = length(dir([strip_dir,'/*_meta.txt']));
            if v4_nseg == 0
                continue;
            end
            v4_seg_arr_temp = [];
            v41_seg_arr_temp = [];
            for i = 1:v4_nseg
                v4_seg_arr_temp = [v4_seg_arr_temp; qc_this_stripID_seg_arr];
                v41_seg_arr_temp = [v41_seg_arr_temp; i * ones(length(qc_this_stripID_seg_arr), 1)];
            end
            v4_seg_arr = v4_seg_arr_temp;
            v41_seg_arr = v41_seg_arr_temp;
        end

        v41_seg_unique_arr = unique(v41_seg_arr);
        for v41_seg_unique_dx = 1:length(v41_seg_unique_arr)
            v41_seg_seg = v41_seg_unique_arr(v41_seg_unique_dx);
            if v41_seg_seg == -1
                continue;
            end

            v4_seg_map_to_this_v41_seg = v4_seg_arr(find(v41_seg_arr == v41_seg_seg));
            [~,ia,~] = intersect(qc_this_stripID_seg_arr, v4_seg_map_to_this_v41_seg);
            qc_this_stripID_and_v41_seg_idx_arr = qc_this_stripID_idx_arr(ia);

            if isempty(qc_this_stripID_and_v41_seg_idx_arr)
                continue;
            end

            v41_seg_flag = 0;
            v41_seg_x = [];
            v41_seg_y = [];

            num_seg_with_qc = 0;
            for i = 1:length(qc_this_stripID_and_v41_seg_idx_arr)
                qc_idx = qc_this_stripID_and_v41_seg_idx_arr(i);

                v4_seg_stripID = qc.stripID{qc_idx};
                v4_seg_seg = qc.seg(qc_idx);
                v4_seg_flag = qc.flag(qc_idx);
                v4_seg_x = qc.x{qc_idx};
                v4_seg_y = qc.y{qc_idx};

%                if ~isempty(v4_seg_x)
%                    num_seg_with_qc = num_seg_with_qc + 1;
%                    if num_seg_with_qc > 1
%                        fprintf('Found v4 qc for multiple segments of strip; %s, stripID: %s\n', qc_file, stripID);
%                    end
%                end

                if ~strcmp(v4_seg_stripID, stripID)
                    error('v4 segment stripID (%s) is not equal to expected working stripID (%s)', v4_seg_stripID, stripID);
                end
                if ~ismember(v4_seg_seg, v4_seg_map_to_this_v41_seg)
                    error('v4 segment number (%d) is not parted of expected working segment set', v4_seg_seg);
                end

                if v4_seg_flag == 3
                    v41_seg_flag = 3;
                else
                    if ~isempty(v4_seg_x) && v4_seg_flag ~= 5
                        error('qc exists with flag=%d; %s, stripID: %s, seg: %d', v4_seg_flag, qc_file, stripID, v4_seg_seg);
                    end
                    if v41_seg_flag ~= 3
                        % FIXME: The following comparison should be inverted.
                        % -t    But in addition, a QC mask should be added
                        % -t    with the footprint of any '5' flag v4 segments.
                        if v4_seg_flag > v41_seg_flag
                            v41_seg_flag = v4_seg_flag;
                        end
                    end
                end

                v41_seg_x = [v41_seg_x, v4_seg_x];
                v41_seg_y = [v41_seg_y, v4_seg_y];

            end

            qc_out.stripID{end+1} = stripID;
            qc_out.seg(end+1) = v41_seg_seg;
            qc_out.flag(end+1) = v41_seg_flag;
            qc_out.x{end+1} = v41_seg_x;
            qc_out.y{end+1} = v41_seg_y;

        end
    end

    fprintf('Writing %d qc records to v41 qc file: %s\n',length(qc_out.stripID),qc_out_file);
    save(qc_out_file,'-struct','qc_out','-v7.3');
end


function [v4_seg, v41_seg] = readRemergeInfoFile(fileName)

    fid = fopen(fileName, 'r');
    nlines = 0;
    tline = fgetl(fid);
    while ischar(tline)
        nlines = nlines + 1;
        tline = fgetl(fid);
    end
    fclose(fid);

    fid = fopen(fileName, 'r');
    arr = textscan(fid, 'seg%f>%s', 'CommentStyle','(');
    fclose(fid);

    v4_seg = arr{1};
    v41_seg = arr{2};
    v41_seg = strrep(v41_seg, 'seg', '');
    v41_seg = strrep(v41_seg, 'none', '-1');
    v41_seg = cellfun(@(x) str2num(x), v41_seg);

    if ~isequal(size(v4_seg), size(v41_seg))
        error('Unequal size of left and right segment mapping parsed from remerge info file: %s', fileName);
    end
    if isempty(v4_seg)
        error('Could not parse segment numbers from remerge info file: %s', fileName);
    end
    if length(v4_seg) ~= nlines
        error('Number of segments parsed from remerge info file (%d) does not match number of lines (%d): %s', length(v4_seg), nlines, fileName);
    end

end
