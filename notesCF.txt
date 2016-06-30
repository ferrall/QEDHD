- I've split it into a simple main file and then a humdim file to
be imported.

-Eventually humdim should be split into shorter files.

-The stuff that was commented out in main() has been split into
    separate functions. Some errors will occur because they probably
    refer to variables still local to main() ... to get the code to run
    they will have to be declared globally.

-Back in Ox v4 you could not reference a row vector like "x[3]".
    You had to treat it like a matrix, so "x[3][0]". This meant a lot
    of extra [0] and []. Now in Ox v7, if x is currently a vector, x[3]
    is fine. And if one of the dimension is not 1 Ox will produce a
    warning. So I'm getting rid of all the extra brackets, but if I get
    rid of one that is needed you'll get a warning. (You might have to
    find the line in human

-Once it looks a little cleaner I will start re-organizing all the
    separate parameters and outcomes into arrays and/or objects. I
    won't finish the job, just suggest a path and show the Ox syntax
    for stuff.

-Cosmetic changes (indentation, moved comments)
