clear all;

currentPath = pwd;
pala_local = false;
% data_folder = [currentPath,'\localization\'];

if pala_local
    data_folder = [currentPath,'\pala\'];
else
    data_folder = [currentPath,'\localization\'];
end
localfiles = dir([data_folder '*.mat']);

max_blood_speed = 500; % mm/s;
pixel_size = 0.1; %mm
time_interval = 0.001; %s
max_link_dis = max_blood_speed*time_interval/pixel_size;
%local = load([localfiles(1).folder filesep localfiles(1).name],'mattracking').mattracking;

Nbuffers = numel(localfiles);
Track_tot = {};
parfor hhh = 1:min(999,Nbuffers)
    local = load([localfiles(hhh).folder filesep localfiles(hhh).name],'mattracking').mattracking;
    
    %if pala local
    if pala_local
        local(:,1:2) = local(:,1:2)*10;
    end
    
    
    sizelocal = size(local);

    numberOfFrames = 800;

    index_frames=arrayfun(@(i) find(local(:,3)==i),[1:numberOfFrames],'UniformOutput',false);

    points=arrayfun(@(i) [local(index_frames{i},1),local(index_frames{i},2)],[1:numberOfFrames],'UniformOutput',false);

    n_cells = cellfun(@(x) size(x, 1), points);
    n_slices = numel(points);
    %max_link_dis = 5;
    adjacency_track = bg_tracking(local,points,max_link_dis);

    min_length = 15;
    n_tracks=numel(adjacency_track);
    count=1;Tracks_raw = {};
    for i_track = 1:n_tracks
        track_id = adjacency_track{i_track};
        track_points = local(track_id,:);
        if length(track_points(:,1))>min_length
            Tracks_raw{count}=track_points;
            count=count+1;
        end
    end

    interp_factor = 1/max_link_dis*.8;
    smooth_factor = 20;
    Tracks_out = {};
    for i_track = 1:size(Tracks_raw,2)
        track_points=double(Tracks_raw{1,i_track});
        xi=track_points(:,2);
        zi=track_points(:,1);
        zu=interp1(1:length(zi),smooth(zi,smooth_factor),1:interp_factor:length(zi));
        xu=interp1(1:length(xi),smooth(xi,smooth_factor),1:interp_factor:length(xi));

        if length(zi)>min_length
            Tracks_out{i_track,1}=cat(2,zu',xu');
        end
    end
    
    Track_tot{hhh} = Tracks_out;
    disp(hhh);
end

Track_tot = cat(1,Track_tot{:});

sizeout = [780,1180]+[1,1];
MatOut = zeros(sizeout);

for itrack = 1:numel(Track_tot)
    % round pixel position into [z,x] index
    pos_z_round = round(Track_tot{itrack}(:,1));
    pos_x_round = round(Track_tot{itrack}(:,2));

    % remove out of grid bubbles (ie. the grid is too small)
    outP = zeros(numel(pos_z_round),1);
    outP = or(outP,pos_z_round<1);outP = or(outP,pos_z_round>sizeout(1));
    outP = or(outP,pos_x_round<1);outP = or(outP,pos_x_round>sizeout(2));

    ind=(sub2ind(sizeout,pos_z_round(~outP),pos_x_round(~outP)));
    % Tracks are counted only once per pixel, the unique keeps only 1 point per pixel
    ind = unique(ind);

    for j=1:size(ind,1)
        % increment the intensity count of the pixel crossed by a track by +1
        MatOut(ind(j))=MatOut(ind(j))+1;
    end
end

imshow(MatOut);
if pala_local
    save('matout_pala_local.mat','MatOut');
else
    save('matout.mat','MatOut');
end


function adjacency_track = bg_tracking(local,points, max_link_dis)
    sizelocal = size(local);
    A = sparse(sizelocal(1),sizelocal(1));
    
    n_cells = cellfun(@(x) size(x, 1), points);
    n_slices = numel(points);
    
    ecount = zeros(n_slices-1,1);
    current_slice_index = 0;
    for i = 1 : n_slices-1
        edge = {};
        next_pset = points{i+1};
        for j = 1 : n_cells(i)
            microb = points{i}(j,:);
            for k = 1:n_cells(i+1)
                dis = norm(microb-next_pset(k,:));
                if dis <= max_link_dis
                    edge(end+1,:) = {current_slice_index+j, current_slice_index+n_cells(i)+k, dis};
                end
            end
        end
        
        edge_sort = sortrows(edge,3);
        esize = size(edge_sort);
        
        flag = zeros(n_cells(i)+n_cells(i+1),1);
        for j = 1 : esize(1)
            e = edge_sort(j,:);
            cindex = e{1} - current_slice_index;
            nindex = e{2} - current_slice_index;
            if flag(cindex) ==0 && flag(nindex) == 0
                ecount(i) = ecount(i)+1;
                flag(cindex) = 1;
                flag(nindex) = 1;
                A(e{1}, e{2}) = 1;
            end
        end
        
        current_slice_index = current_slice_index+n_cells(i);
    end
    
    microb_without_source = [];
    for j = 1 : size(A, 2)
        if isempty(find(A(:,j)==1, 1))
            microb_without_source = [ microb_without_source ; j ];
        end
    end
    n_tracks = numel(microb_without_source);
    adjacency_tracks = cell(n_tracks, 1);

    AT = A';
    for i = 1 : n_tracks
        
        tmp_holder = NaN(sizelocal(1), 1);
        
        target = microb_without_source(i);
        index = 1;
        while ~isempty(target)
            tmp_holder(index) = target;
            target = find( AT(:, target), 1, 'first' );
            index = index + 1;
        end
        
        adjacency_tracks{i} = tmp_holder ( ~isnan(tmp_holder) );
    end
    adjacency_track = adjacency_tracks;
end

