function [x,y] = andrewAH(x,y,N1,N2,w,tol,n) 
u = zeros(N2,N1);

xp = x;
yp = y;
%SOR loop
for k = 1:n

    for i = 2:N2-1
        for j = 2:N1-1
            x(i,j) = (1/4)*(xp(i,j+1)+x(i,j-1)+xp(i+1,j)+x(i-1,j));
            y(i,j) = (1/4)*(yp(i,j+1)+y(i,j-1)+yp(i+1,j)+y(i-1,j));


            x(i,j) = w*x(i,j) + (1-w)*xp(i,j);
            y(i,j) = w*y(i,j) + (1-w)*yp(i,j);



        end
    end

    maxerrx = 0;
    maxerry = 0;
    for i = 2:N2-1
        for j = 2:N1-1
            xtol = x(i,j)-xp(i,j);
            ytol = y(i,j)-yp(i,j);
            
            if(maxerrx > xtol)
                maxerrx = xtol;
            end
            if(maxerry > ytol)
                maxerry = ytol;
            end
        end
    end
    
    if(sqrt(xtol^2 + ytol^2)) < tol
        disp(k)
        break
    end
    
    
    xp = x; 
    yp = y;
    orient landscape;
    mesh(x,y,u,'EdgeColor','black');
    grid on;
    xlabel('x')
    ylabel('y')
    axis square;
    view(0,90);
    pause(0.00001)

end