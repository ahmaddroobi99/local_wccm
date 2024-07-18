function [hh, mm, tt] = sec2hhmmss(seconds)
    h = fix(mod(seconds, [0, 3600, 60]) ./ [3600, 60, 1]);
    hh = h(1);
    mm = h(2);
    tt = h(3);
end