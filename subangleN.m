
function [angleN]=subangleN(x,y)

    if x >0 
        if y>=0
            angleN=atan(y/x);
        else
            angleN=atan(y/x)+2*pi;
        end
    end
    if (x==0 & y>=0)
        angleN=pi/2;
    end
    if (x==0 & y<=0)
        angleN=3*pi/2;
    end
    if x<0
        angleN=atan(y/x)+pi;
    end
    