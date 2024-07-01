function plot_simulated_cell_custom(pv)
    % Plots rate map of simulated cell from simulate_cell_FR.m
    % Specify desired parameters within simulate_cell_FR.m before calling
    % this function.
    
    % read in behavioural data from vmpv object
    stc = pv.data.sessionTimeC;
    tbin_size = 0.05;  % size of time bin
    
    %%% Specify environment bin geometry here %%%
    % can be changed
    floor_width = 40;
    wall_height = 8;
    pillar_height = 5;
    num_hd_bins = 60;
    
    % do not change, dependent on other variables/environmental geometry/default setting
    wall_width = floor_width; pillar_width = floor_width/5;
    viewbin_offset = 2;
    num_place_bins = floor_width^2;
    num_view_bins = viewbin_offset + 2*floor_width^2 + wall_height*wall_width*4 + 4*pillar_height*pillar_width*4;
    
    % Code adapted from plotgridmap.m
    exploded = false;  % exploded view, i.e. floor and ceiling are separated from walls/pillars
    
    floor_x = repmat(0:floor_width, floor_width+1, 1);
    floor_y = flipud(repmat([0:floor_width]', 1, floor_width+1));
    floor_z = zeros(floor_width+1,floor_width+1);

    ceiling_x = floor_x;
    ceiling_y = floor_y;

    walls_x = repmat([0.*ones(1,wall_width) 0:wall_width-1 wall_width.*ones(1,wall_width) wall_width:-1:0], wall_height+1, 1);
    walls_y = repmat([0:wall_width-1 wall_width.*ones(1,wall_width) wall_width:-1:1 0.*ones(1,wall_width+1)], wall_height+1, 1);

    P1_x = repmat([3*pillar_width.*ones(1,pillar_width) 3*pillar_width:4*pillar_width-1 4*pillar_width.*ones(1,pillar_width) 4*pillar_width:-1:3*pillar_width], pillar_height+1, 1);
    P1_y = repmat([pillar_width:2*pillar_width-1 2*pillar_width.*ones(1,pillar_width) 2*pillar_width:-1:pillar_width+1 pillar_width.*ones(1,pillar_width+1)], pillar_height+1, 1);

    P2_x = repmat([pillar_width.*ones(1,pillar_width) pillar_width:2*pillar_width-1 2*pillar_width.*ones(1,pillar_width) 2*pillar_width:-1:pillar_width], pillar_height+1, 1);
    P2_y = P1_y;

    P3_x = P1_x;
    P3_y = repmat([3*pillar_width:4*pillar_width-1 4*pillar_width.*ones(1,pillar_width) 4*pillar_width:-1:3*pillar_width+1 3*pillar_width.*ones(1,pillar_width+1)], pillar_height+1, 1);

    P4_x = P2_x;
    P4_y = P3_y;
    
    if exploded
        ceiling_z = 3*wall_height.*ones(floor_width+1,floor_width+1);
        walls_z = repmat([2*wall_height:-1:wall_height]', 1, wall_width*4 + 1);
        PX_z = repmat([wall_height+pillar_height:-1:wall_height]', 1, pillar_width*4 + 1);
    else
        ceiling_z = wall_height.*ones(floor_width+1,floor_width+1);
        walls_z = repmat([wall_height:-1:0]', 1, wall_width*4 + 1);
        PX_z = repmat([pillar_height:-1:0]', 1, pillar_width*4 + 1);
    end

    floor_last_bin = floor_width^2 + viewbin_offset;
    floor = flipud(reshape(viewbin_offset+1:floor_last_bin, floor_width, floor_width)');

    % ceiling follows floor mapping, top down view
    ceiling_last_bin = floor_last_bin + floor_width^2;
    ceiling = flipud(reshape(floor_last_bin+1:ceiling_last_bin, floor_width, floor_width)');

    % from top down, slit walls at bottom left corner, open outwards.
    % start from row closest to ground, rightwards, then climb rows
    walls_last_bin = ceiling_last_bin + 4*wall_width*wall_height;
    walls = flipud(reshape(ceiling_last_bin+1:walls_last_bin, wall_width*4, wall_height)');

    % BL - bottom left, and so on, from top view, same slicing as walls
    % pillar width 8, height 5
    P1_last_bin = walls_last_bin + 4*pillar_width*pillar_height;
    P2_last_bin = P1_last_bin + 4*pillar_width*pillar_height;
    P3_last_bin = P2_last_bin + 4*pillar_width*pillar_height;
    P1_BR = flipud(reshape(walls_last_bin+1:P1_last_bin, pillar_width*4, pillar_height)');
    P2_BL = flipud(reshape(P1_last_bin+1:P2_last_bin, pillar_width*4, pillar_height)');
    P3_TR = flipud(reshape(P2_last_bin+1:P3_last_bin, pillar_width*4, pillar_height)');
    P4_TL = flipud(reshape(P3_last_bin+1:num_view_bins, pillar_width*4, pillar_height)');

    % mark out place bins under pillars
    under_pillars = nan(4*pillar_width^2, 1); j = 1;
    for i = [pillar_width+1:2*pillar_width, 3*pillar_width+1:4*pillar_width]
        under_pillars(j:j+2*pillar_width-1) = floor_width*(i-1) + [pillar_width+1:2*pillar_width, 3*pillar_width+1:4*pillar_width];
        j = j+2*pillar_width;
    end
    
    % Preprocessing on vmpv data (binning into tbin_size time bins + applying vel/minObs filters)
    % using vel filter only for now
    ThresVel = 1;
    conditions = ones(size(stc,1),1);
    conditions = conditions & get(pv,'SpeedLimit',ThresVel);

    % don't touch for now
    UseMinObs = false;
    if UseMinObs
        place_bins_sieved = pv.data.place_good_bins;
        view_bins_sieved = pv.data.view_good_bins;
        conditions = conditions & (pv.data.pv_good_rows);
    else
        place_bins_sieved = 1:(40 * 40);
        view_bins_sieved = 1:5122;
    end

    % Construct new stc, with each row representing a time bin
    bin_stc = nan(2*size(stc,1),4);
    bin_stc(1,1:4) = stc(2,1:4);

    current_tbin = 1; % refers to bin_stc row being filled
    stc_last = 2;

    dstc = diff(stc(:,1));
    cstc = find(dstc ~= 0) + 1;
    cstc(1:end-1,2) = cstc(2:end,1) - 1;
    cstc(end,2) = size(stc,1);
    cstc_track = 1;

    while bin_stc(current_tbin, 1) + tbin_size <= stc(end,1)
        if ~conditions(stc_last) % do not allow any bin to include any ~condition stc rows
            bin_stc(current_tbin, 1) = stc(stc_last+1,1);
            stc_last = stc_last + 1;
            continue
        else
            while bin_stc(current_tbin, 1) < stc(stc_last+1,1)
                while ~any(cstc(cstc_track,1):cstc(cstc_track,2) == stc_last)
                    cstc_track = cstc_track + 1;
                end
                match_idx = cstc(cstc_track,1):cstc(cstc_track,2); % to account for multiple simultaneously occupied bins
                match_b_idx = current_tbin:current_tbin-1+length(match_idx);
                if length(match_idx) > 1
                    bin_stc(match_b_idx, 1) = bin_stc(current_tbin, 1);
                    for i = 1:length(match_idx)
                        bin_stc(current_tbin+i-1, 2:4) = stc(match_idx(i), 2:4);
                    end
                    current_tbin = current_tbin + length(match_idx) - 1; % update index of latest filled tbin
                else
                    bin_stc(current_tbin, 2:4) = stc(stc_last, 2:4);
                end

                bin_stc(current_tbin+1, 1) = bin_stc(current_tbin, 1) + tbin_size;
                current_tbin = current_tbin + 1;
            end
            while stc_last+1 <= size(stc,1) && stc(stc_last+1,1) <= bin_stc(current_tbin, 1)
                stc_last = stc_last + 1;
            end
        end
    end

    bin_stc = bin_stc(~isnan(bin_stc(:,1)),:);
    % place/view bin sieving
    rows_remove = [];
    for k = 1:size(bin_stc,1)
        if ~any(place_bins_sieved == bin_stc(k,2)) || ~any(view_bins_sieved == bin_stc(k,4))
            rows_remove = [rows_remove k];
        end
    end
    bin_stc(rows_remove,:) = [];

    bin_dstc = diff(bin_stc(:,1));
    bin_cstc = find(bin_dstc ~= 0) + 1;
    bin_cstc(1:end-1,2) = bin_cstc(2:end,1) - 1;
    bin_cstc(end,2) = size(bin_stc,1);
    
    % simulate firing rates on behavioural data
    place_ratemap = zeros(num_place_bins,1); place_durmap = zeros(num_place_bins,1);
    hd_ratemap = zeros(num_hd_bins,1); hd_durmap = zeros(num_hd_bins,1);
    view_ratemap = zeros(num_view_bins,1); view_durmap = zeros(num_view_bins,1);
    
    firing_rates = simulate_cell_FR_custom(bin_stc(:,2:4));
    for k = 1:length(firing_rates)
        place_ratemap(bin_stc(k,2)) = place_ratemap(bin_stc(k,2)) + firing_rates(k)*tbin_size;
        place_durmap(bin_stc(k,2)) = place_durmap(bin_stc(k,2)) + tbin_size;
        hd_ratemap(bin_stc(k,3)) = hd_ratemap(bin_stc(k,3)) + firing_rates(k)*tbin_size;
        hd_durmap(bin_stc(k,3)) = hd_durmap(bin_stc(k,3)) + tbin_size;
        view_ratemap(bin_stc(k,4)) = view_ratemap(bin_stc(k,4)) + firing_rates(k)*tbin_size;
        view_durmap(bin_stc(k,4)) = view_durmap(bin_stc(k,4)) + tbin_size;
    end
        
    place_nonzero = find(place_durmap > 0);
    hd_nonzero = find(hd_durmap > 0);
    view_nonzero = find(view_durmap > 0);
    place_ratemap(place_nonzero) = place_ratemap(place_nonzero) ./ place_durmap(place_nonzero);
    hd_ratemap(hd_nonzero) = hd_ratemap(hd_nonzero) ./ hd_durmap(hd_nonzero);
    view_ratemap(view_nonzero) = view_ratemap(view_nonzero) ./ view_durmap(view_nonzero);
    
    
    % Generate plots for each spatial variable
    % Place plot
    fp = figure('Name','Place plot');
    surf(floor_x, floor_y, floor_z, flipud(reshape(place_ratemap(1:1600), 40, 40)'));
    alpha 1; shading flat;
    zlim([0,1]);
    view(-35,20);
    colormap jet;
    colorbar;

    rectangle('Position', [8, 8, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position', [8, 24, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position', [24, 8, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position', [24, 24, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);

    % Head direction plot
    fh = figure('Name','Head direction plot');
    pax = polaraxes();
    polarplot(deg2rad((0:60)*360/60), [hd_ratemap; hd_ratemap(1)]);
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';

    % View plot
    fv = figure('Name','View plot');
    
    % Plot floor
    surf(floor_x, floor_y, floor_z, flipud(reshape(view_ratemap(3:1600+3-1), 40, 40)'));
    alpha 0.35; shading flat;
    hold on;

    % Plot ceiling and walls
    surf(ceiling_x, ceiling_y, ceiling_z, flipud(reshape(view_ratemap(1603:1603+1600-1), 40, 40)'));
    alpha 0.35; shading flat;
    surf(walls_x, walls_y, walls_z, flipud(reshape(view_ratemap(3203:3203+1280-1), 40*4, 8)'));      
    alpha 0.35; shading flat;

    % Plot pillars
    surf(P1_x, P1_y, PX_z, flipud(reshape(view_ratemap(4483:4483+160-1), 8*4, 5)'));
    alpha 0.35; shading flat;
    surf(P2_x, P2_y, PX_z, flipud(reshape(view_ratemap(4643:4643+160-1), 8*4, 5)'));
    alpha 0.35; shading flat;
    surf(P3_x, P3_y, PX_z, flipud(reshape(view_ratemap(4803:4803+160-1), 8*4, 5)'));
    alpha 0.35; shading flat;
    surf(P4_x, P4_y, PX_z, flipud(reshape(view_ratemap(4963:4963+160-1), 8*4, 5)'));
    alpha 0.35; shading flat; 
    view(-35,20);
    colormap jet;
    colorbar;

    rectangle('Position', [8, 8, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position', [8, 24, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position', [24, 8, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);
    rectangle('Position', [24, 24, 8, 8], 'EdgeColor', 'k', 'LineWidth', 1);

end
