function txtprogressbar(fraction_done)
    % show simulation progress and remaining time as text
    % This program is useful when using Matlab in console mode (under Linux)
    % for long simulations (usually several hours). Launch Matlab in 
    % the background with the following command:
    % 
    % nohup matlab -nojvm -nosplash -nodesktop < my_program.m &
    % 
    % All text output is directed to nohup.out file, which can be tested 
    % periodically with: "tail -f nohup.out". Using "nohup" has also the 
    % advantage that Matlab continues to run after the user is logged out.
    % Demo:
    %   n = 1000;
    %   txtprogressbar % Set starting time
    %   for i = 1:n
    %       pause(0.01) % Do something important
    %       txtprogressbar(i/n) % Update text
    %   end
    % author: Bogdan Cristea
    % revision date: 13.05.2007
    % revision date: 23.05.2007 - help added, now simulation progress has one 
    % digit after the decimal point
    % thanks: R. S. Schestowitz, University of Manchester (for text version) and
    %         Steve Hoelzer (for graphical version)
    persistent start_time
    % Set defaults if fraction_done not passed in
    if nargin < 1
        fraction_done = 0;
        percent_done = 0;
        start_time = clock;
    end
    % Update progress message
    if (fraction_done == 0)
        progress_msg = ' 0%';
    else
        percent_done = 100*fraction_done;
        run_time = etime(clock, start_time);
        time_left = run_time/fraction_done - run_time;
        progress_msg = sprintf('%2.1f%%    %s remaining', percent_done, sec2timestr(time_left));
    end
    % Display text
    clc;
    disp(progress_msg);
    % disp('|==================================================|');%uncomment these lines if you want a progress bar
    % disp(['|' char('#'*ones(1, floor(percent_done/2))) '|'])
    % disp('|==================================================|');
    % If task completed, clear vars, then exit
    if percent_done==100 % Task completed
        clear start_time% Clear persistent vars
        return
    end
    % ------------------------------------------------------------------------------
    function timestr = sec2timestr(sec)
    % Convert a time measurement from seconds into a human readable string.
    % Convert seconds to other units
    d = floor(sec/86400); % Days
    sec = sec - d*86400;
    h = floor(sec/3600); % Hours
    sec = sec - h*3600;
    m = floor(sec/60); % Minutes
    sec = sec - m*60;
    s = floor(sec); % Seconds
    % Create time string
    if d > 0
        if d > 9
            timestr = sprintf('%d day',d);
        else
            timestr = sprintf('%d day, %d hr',d,h);
        end
    elseif h > 0
        if h > 9
            timestr = sprintf('%d hr',h);
        else
            timestr = sprintf('%d hr, %d min',h,m);
        end
    elseif m > 0
        if m > 9
            timestr = sprintf('%d min',m);
        else
            timestr = sprintf('%d min, %d sec',m,s);
        end
    else
        timestr = sprintf('%d sec',s);
    end
    