cases = 's'; % Assuming 'cases' is a character array, not a cell array
fsavefigs = true;

for casename = cases
    
    datafolder = fullfile('data', casename);
    figfolder = fullfile('zzzz', casename);
    if fsavefigs && ~isdir(figfolder)
        mkdir(figfolder);
    end
    
    dataext = load(fullfile('diag11.mat'));
    
    qhc = coarsenqh(datafolder, 0.02); % load data, second argument is junk
    
    Ntruth = sqrt(size(qhc, 1));
    k = [0:Ntruth/2, -Ntruth/2+1:-1]';
    [KX, KY] = meshgrid(k, k);
    Knorm = hypot(KX, KY);
    Knorminv = 1 ./ Knorm;
    Knorminv(1, 1) = 0;
    
    kp = Knorm(1, 1:Ntruth/2+1);
    
    psi_hat_qhc = -Knorminv(:) .* qhc / Ntruth^2;

    t_columns = size(psi_hat_qhc,2)
    t = 10^5;
    dt = 0.1;
    t/dt;
    tslice= 100/dt ;

    psi_hat_qhc(:,1:tslice)
    qhcm = abs(mean(psi_hat_qhc, 2)).^2;
    qhcv = var(psi_hat_qhc, 0, 2);
    
    qhcm = reshape(qhcm, Ntruth, Ntruth);
    qhcv = reshape(qhcv, Ntruth, Ntruth);
    
    ev = zeros(Ntruth/2+1, 1);
    em = zeros(Ntruth/2+1, 1);
    for jj = 1:Ntruth
        for ii = 1:Ntruth
            k = sqrt(KX(ii, jj).^2 + KY(ii, jj).^2);
            if ceil(k) <= Ntruth/2
                r = k-floor(k);
                % plus 1 because of MATLAB indexing
                ev(floor(k)+1) = ev(floor(k)+1) + (1-r) * abs(qhcv(ii, jj));
                ev(ceil(k)+1) = ev(ceil(k)+1) + r * abs(qhcv(ii, jj));
                
                em(floor(k)+1) = em(floor(k)+1) + (1-r) * abs(qhcm(ii, jj));
                em(ceil(k)+1) = em(ceil(k)+1) + r * abs(qhcm(ii, jj));
            end
        end
    end
    ev = 0.5 * ev;
    em = 0.5 * em;
    
    % Divide data into three segments
    seg_length = floor(length(kp) / 7);
    kp_segments = {kp(1:seg_length), kp(1:2*seg_length), kp(1:3*seg_length),kp(1:4*seg_length), kp(1:5*seg_length), kp(1:6*seg_length),kp(1:end)};
    em_segments = {em(1:seg_length), em(1:2*seg_length), em(1:3*seg_length),em(1:4*seg_length), em(1:5*seg_length), em(1:6*seg_length),em(1:end)};
    ev_segments = {ev(1:seg_length), ev(1:2*seg_length), ev(1:3*seg_length),ev(1:4*seg_length), ev(1:5*seg_length), ev(1:6*seg_length),ev(1:end)};
    
    % for seg_idx = 1:3
        handle = figure;
        loglog(kp_segments{1}, em_segments{1} + ev_segments{1}, 'LineWidth', 1.25);
        hold on;

        loglog(kp_segments{2}, em_segments{2} + ev_segments{2}, 'LineWidth', 1.25);
        loglog(kp_segments{3}, em_segments{3} + ev_segments{3}, 'LineWidth', 1.25);
        loglog(kp_segments{4}, em_segments{4} + ev_segments{4}, 'LineWidth', 0.75);
        loglog(kp_segments{5}, em_segments{5} + ev_segments{5}, 'LineWidth', 1.25);
        loglog(kp_segments{6}, em_segments{6} + ev_segments{6}, 'LineWidth', 1);
        loglog(kp_segments{7}, em_segments{7} + ev_segments{7}, 'LineWidth', 0.5);

        loglog(kp_segments{1}, ev_segments{1}); 
        loglog(kp_segments{2}, ev_segments{2}); 
        loglog(kp_segments{3}, ev_segments{3}); 
        loglog(kp_segments{4}, ev_segments{4}); 
        loglog(kp_segments{5}, ev_segments{5}); 
        loglog(kp_segments{6}, ev_segments{6}); 
        loglog(kp_segments{7}, ev_segments{7}); 


        loglog(kp_segments{1}, em_segments{1});
        loglog(kp_segments{2}, em_segments{2});
        loglog(kp_segments{3}, em_segments{3});
        loglog(kp_segments{4}, em_segments{4});
        loglog(kp_segments{5}, em_segments{5});
        loglog(kp_segments{6}, em_segments{6});
        loglog(kp_segments{7}, em_segments{7});

        legend('total7', 'variance7', 'mean7');
        
        if strcmp(dataext.params.forcingtype, 'constant')
            kkp = (1:35)';
            alpha = -5/3;
        else
            kkp = (1:45)';
            alpha = -3;
        end
        seg_idx =1;
        [~, ~, ikp] = intersect(kkp, 30);
        kkp_interp = kkp(ikp); % Interpolate kkp to match dimensions with ikp
        ev_interp = interp1(kp_segments{seg_idx}, ev_segments{seg_idx}, kkp_interp); % Interpolate ev_segments to match dimensions with kkp_interp
        amin = fminsearch(@(a) mean(log(abs(a * kkp_interp.^alpha - ev_interp))), 1);
                
        loglog(kkp, amin * kkp.^alpha, 'k--');
        text(2, amin * kkp(1).^alpha, strcat('k^{', strtrim(rats(alpha)), '}'));
        
        xlabel('wavenumber');
        ylabel('energy');
      
        

    % end
      hold off; 
        drawnow; 
        if fsavefigs
            resizefonts(handle);
            export_fig(fullfile(figfolder, sprintf('spectrum_segment_seven')), '-q101', '-p.01', '-m1', '-pdf');
        end
                % Pause to allow observation of the figure before proceeding
        pause;
end
