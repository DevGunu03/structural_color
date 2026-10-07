function p = repo_root()
%REPO_ROOT Folder that holds data/, refs/ and matlab/.
p = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
