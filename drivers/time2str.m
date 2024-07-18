function timestr = time2str(t)
    if t == Inf
        timestr = sprintf(' --:--:--');
    else
        [hh, mm, tt] = sec2hhmmss(t);
        timestr = sprintf('%02d:%02d:%02d', hh, mm, tt);
    end
end