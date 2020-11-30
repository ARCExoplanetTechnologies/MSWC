function [runDir dmDir imagDir] = make_dirs(runNum)
    runDir = ['run' num2str(runNum) '\'];
    mkdir(runDir)

    dmDir = [runDir 'DM\'];
    mkdir(dmDir);

    imagDir = [runDir 'IMAGES\'];
	mkdir(imagDir);
    
end

