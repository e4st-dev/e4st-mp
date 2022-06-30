function root_dir = e4st_root()
    cur_fn = fileparts(mfilename('fullpath'));
    root_dir = fullfile(cur_fn, '..', '..', '..');
end


