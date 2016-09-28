function [x,y] = andrewWINSLOW(x,y,N1,N2,w,tol,n)  
u = zeros(N1, N2);

xp = x;
yp = y;
%SOR loop
for k = 1:n
    
    
    for i = 1:N1-1
        for j = 2:N2-1
            
            if i == 1
                x_eta = (xp(i,j+1) - x(i,j-1))/2;
                y_eta = (yp(i,j+1) - y(i,j-1))/2;

                y_psi = (yp(i+1,j)-y(N1-1,j))/2;
                x_psi = (xp(i+1,j)-x(N1-1,j))/2;

                alpha = x_eta^2 + y_eta^2;
                beta = x_psi*x_eta + y_psi*y_eta;
                gamma = x_psi^2 + y_psi^2;


                x(i,j) = (1/(2*(alpha+gamma)))*(alpha*(xp(i+1,j)+xp(N1-1,j))+ ...
                gamma*(xp(i,j+1)+x(i,j-1)) - (1/2)*beta*(xp(i+1,j+1)-x(i+1,j-1)-xp(N1-1,j+1)+xp(N1-1,j-1)));
                y(i,j) = (1/(2*(alpha+gamma)))*(alpha*(yp(i+1,j)+yp(N1-1,j))+ ...
                gamma*(yp(i,j+1)+y(i,j-1)) - (1/2)*beta*(yp(i+1,j+1)-y(i+1,j-1)-yp(N1-1,j+1)+yp(N1-1,j-1)));


            else
                x_eta = (xp(i,j+1) - x(i,j-1))/2;
                y_eta = (yp(i,j+1) - y(i,j-1))/2;

                y_psi = (yp(i+1,j)-y(i-1,j))/2;
                x_psi = (xp(i+1,j)-x(i-1,j))/2;

                alpha = x_eta^2 + y_eta^2;
                beta = x_psi*x_eta + y_psi*y_eta;
                gamma = x_psi^2 + y_psi^2;


                x(i,j) = (1/(2*(alpha+gamma)))*(alpha*(xp(i+1,j)+x(i-1,j))+ ...
                gamma*(xp(i,j+1)+x(i,j-1)) - (1/2)*beta*(xp(i+1,j+1)-x(i+1,j-1)-x(i-1,j+1)+x(i-1,j-1)));
                y(i,j) = (1/(2*(alpha+gamma)))*(alpha*(yp(i+1,j)+y(i-1,j))+ ...
                gamma*(yp(i,j+1)+y(i,j-1)) - (1/2)*beta*(yp(i+1,j+1)-y(i+1,j-1)-y(i-1,j+1)+y(i-1,j-1)));


                
            end
            
            
            x(i,j) = w*x(i,j) + (1-w)*xp(i,j);
            y(i,j) = w*y(i,j) + (1-w)*yp(i,j);

        end
    end

    for j = 2:N2-1
        x(N1,j) = x(1,j);
        y(N1,j) = y(1,j);
    end

    
    maxerrx = 0;
    maxerry = 0;
    for i = 2:N1-1
        for j = 2:N2-1
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
%     orient landscape;
%     mesh(x,y,u,'EdgeColor','black');
%     grid on;
%     xlabel('x')
%     ylabel('y')
%     axis square;
%     view(0,90);
%     pause(0.0001)

end