function v = lctoolboxversion()

% Run git describe and capture output
thispath = fileparts(mfilename('fullpath'));
[status, out] = system(sprintf('cd "%s" && git describe --long --always', thispath)); % make sure to navigate to the folder of version first!

if status == 0
    v = strtrim(out);
    fprintf("Detected version from git: %s\n", v);
else
    error("Git command failed.\nWorking directory: %s\nError: %s", pwd, out);
end

end