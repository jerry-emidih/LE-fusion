function textprogressbar(c)
    % textprogressbar() displays a text-based progress bar on the console 
    % without the need of a graphical interface. This is a good option for
    % no-GUI instances of MATLAB. It also displays the estimated finishing 
    % time for the process and recalculates it on each update.
    %
    % EXAMPLE OF USE:
    % textprogressbar('Process'); % initializes the progress bar with a
    % custom string
    % textprogressbar(P); % updates the progress bar at P and recalculates
    % the estimated time
    % textprogressbar(100); % finalizes the current progress bar and resets
    %
    % Author: Miguel A. Lago (e-mail: miguelangel.lago (at) gmail (dot) com)
    % Version: 1.0
    %
    % Revisions:
    % - 2019-Feb-08: Initial function
    %
    % Based on by: textprogressbar by Paul Proteus (https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar)
    % and progressbar by Steve Hoelzer (https://www.mathworks.com/matlabcentral/fileexchange/6922-progressbar)
    %% Initialization
    persistent deleteString stringProcess initTime prevC;           %   Carriage return pesistent variable
    % Visualization parameters
    strPercentageLength = 10;   %   Length of percentage string (must be >5)
    strDotsMaximum      = 25;   %   The total number of dots in a progress bar
    %% Main 
    if isempty(stringProcess) && ~ischar(c)
        % If there is no initialization string, use default
        stringProcess='Progress ';
        deleteString=[];
        initTime=now;
        prevC=0;
    elseif isempty(stringProcess) && ischar(c)
        % Initialize with a string
        stringProcess=c;
        deleteString=[];
        initTime=now;
        prevC=0;
    elseif isnumeric(c)
        % Update progress
        elapsedTime=now-initTime;
        etaS=86400*100*elapsedTime/(c-prevC); %in seconds
        etaS=etaS-c*etaS/100;
        if c==100 % Finalizing string
            etaStr=' DONE\n';
        else % Estimated time for humans
            etaStr=[' ' sec2timestr(etaS)];
        end
        percentageOut = [num2str(round(c,2)) '%%'];
        percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
        nDots = floor(c/100*strDotsMaximum);
        dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
        strOut = [stringProcess percentageOut dotOut etaStr];
        if prevC>c % If the script is being run again, do not delete previous characters
            deleteString=[];
        end
        fprintf([deleteString strOut]);
        
        % Update carriage return
        deleteString = repmat('\b',1,length(strOut)-1);
        
        prevC=c;
        initTime=now;
        if c==100
            % Finalize progress bar
            deleteString = [];  
            prevC=0;
        end
    else
        % Any other unexpected input
        error('Unsupported argument type');
    end
end
function timestr = sec2timestr(sec)
    % Convert a time measurement from seconds into a human readable string.
    % Convert seconds to other units
    w = floor(sec/604800); % Weeks
    sec = sec - w*604800;
    d = floor(sec/86400); % Days
    sec = sec - d*86400;
    h = floor(sec/3600); % Hours
    sec = sec - h*3600;
    m = floor(sec/60); % Minutes
    sec = sec - m*60;
    s = floor(sec); % Seconds
    % Create time string
    if w > 0
        if w > 9
            timestr = sprintf('%d week', w);
        else
            timestr = sprintf('%d week, %d day', w, d);
        end
    elseif d > 0
        if d > 9
            timestr = sprintf('%d day', d);
        else
            timestr = sprintf('%d day, %d hr', d, h);
        end
    elseif h > 0
        if h > 9
            timestr = sprintf('%d hr', h);
        else
            timestr = sprintf('%d hr, %d min', h, m);
        end
    elseif m > 0
        if m > 9
            timestr = sprintf('%d min', m);
        else
            timestr = sprintf('%d min, %d sec', m, s);
        end
    else
        timestr = sprintf('%d sec', s);
    end
end

    