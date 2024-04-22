function elapsedDuration = getElapsedDuration(startTime)
% Get the stopwatch time elapsed as a duration in the format HH:MM:SS.
    elapsedDuration = duration(seconds(toc(startTime)), 'Format', 'hh:mm:ss');
end
