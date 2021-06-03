function out = col2trip(color)
% function that maps names of colors to RGB triplets

if color == "orange"
    out = [255,106,0]/255;
elseif color == "teal"
    out = [0,255,255];
elseif color == "purple"
    out = [127,0,110]/255;
elseif color == "darkpurple"
    out = [127,0,55]/255;
elseif color == "pink"
    out = [255, 0, 220]/255;
elseif color == "darkgreen"
    out = [38 127 0]/255;
elseif color == "lightgreen"
    out = [0 255 72]/255;
elseif color == "lightblue"
    out = [155 127 255]/255;
elseif color == "darkblue"
    out = [33 0 107]/255;
elseif color == "darkred"
    out = [127 0 0]/255;
elseif color == "grey"
    out = [96 96 96]/255;
elseif color == "darkgrey"
    out = [48 48 48]/255;
elseif color == "lightgrey"
    out = [160 160 160]/255;
    
    
elseif color == "black"
    out = "black";
elseif color == "red"
    out = "red";
elseif color == "blue"
    out = "blue";
elseif color == "green"
    out = "green";
elseif color == "black"
    out = "black";
elseif color == "magenta"
    out = "magenta";
end

end
   