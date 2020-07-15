function t = LM_getTime(opt,Fs,type)

switch type
    case 'forward'
        t =  (-opt.maxLag:-opt.minLag) / Fs;
    case 'backward'
        t =  (opt.minLag:opt.maxLag) / Fs;
end

end
%
%