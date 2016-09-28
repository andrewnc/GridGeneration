function [ptotal, p2] = final(x,y,N1,N2,n,k,omega,deltat)  
image = sqrt(-1);

r = 6;

firsttol = 10^-7;
secondtol = 10^-2;

%Matrix Initialization
p0 = zeros(N1,N2);
p1 = zeros(N1,N2);
p2 = zeros(N1,N2);

%initial condition for t = 0
time = 0;
for i = 1:N1-1
    p1(i,1) = -(exp(image*k*x(i,1)));

end
    %wrap up
    p1(N1,1) = p1(1,1);

%Time step
for time = 1:n
    
    %Computation at the obstacle boundary
    for i = 1:N1-1
        p2(i,1) = -(exp(image*k*x(i,1))*exp(-image*omega*time*deltat));
    end
    %wrap up
    p2(N1,1) = p2(1,1);
    
    %Computation of interior points
    for i = 1:N1-1
        for j = 2:N2-1
            xeta = (x(i,j+1) - x(i,j-1))/2;
            yeta = (y(i,j+1) - y(i,j-1))/2;
            

            c = 1;
            
            if i == 1
                ypsi = (y(i+1,j)-y(N1-1,j))/2;
                xpsi = (x(i+1,j)-x(N1-1,j))/2;
                
                jacobian = xpsi*yeta-xeta*ypsi;

                alpha = xeta^2 + yeta^2;
                beta = xpsi*xeta + ypsi*yeta;
                gamma = xpsi^2 + ypsi^2;
            
                g = ((c*deltat)/jacobian)^2;
                
                p2(i,j) = p1(N1-1,j)*(alpha*g) + p1(i+1,j)*(alpha*g)+p1(i,j-1)*(gamma*g)+p1(i,j+1)*(gamma*g) ...
                - 2*g*beta*(p1(i+1,j+1)-p1(i+1,j-1)-p1(N1-1,j+1)+p1(N1-1,j-1))/4 - p0(i,j) + p1(i,j)*((-2*alpha*g)-(2*gamma*g)+2);
            else
                ypsi = (y(i+1,j)-y(i-1,j))/2;
                xpsi = (x(i+1,j)-x(i-1,j))/2;  
                
                jacobian = xpsi*yeta-xeta*ypsi;
                
                alpha = xeta^2 + yeta^2;
                beta = xpsi*xeta + ypsi*yeta;
                gamma = xpsi^2 + ypsi^2;
                
                g = ((c*deltat)/jacobian)^2;

                p2(i,j) = p1(i-1,j)*(alpha*g) + p1(i+1,j)*(alpha*g)+p1(i,j-1)*(gamma*g)+p1(i,j+1)*(gamma*g) ...
                - 2*g*beta*(p1(i+1,j+1)-p1(i+1,j-1)-p1(i-1,j+1)+p1(i-1,j-1))/4 - p0(i,j) + p1(i,j)*((-2*alpha*g)-(2*gamma*g)+2);
            end
        end
%         p2(N1,j) = p2(1,j);
    end
    %wrap up
    for j = 2:N2-1
        p2(N1,j) = p2(1,j);
    end
    
    %Computation at artificial boundary
    for i = 1:N1-1
        
        %Backwards scheme to eliminate ghost points
        xeta = (3/2)*x(i,N2) - 2*x(i,N2-1) + (1/2)*x(i,N2-2);
        yeta = (3/2)*y(i,N2) - 2*y(i,N2-1) + (1/2)*y(i,N2-2);
        
        c = 1;
        
        
        if i == 1
            ypsi = (y(i+1,N2)-y(N1-1,N2))/2;
            xpsi = (x(i+1,N2)-x(N1-1,N2))/2;
            
            jacobian = xpsi*yeta-xeta*ypsi;
            g = ((c*deltat)/jacobian)^2;

            alpha = xeta^2 + yeta^2;
            beta = xpsi*xeta + ypsi*yeta;
            gamma = xpsi^2 + ypsi^2;
            
            %Variables used to ease calculation
            delt = c/(r*jacobian);
            lambda = y(i,N2)*xpsi-x(i,N2)*ypsi;
            kay = x(i,N2)*yeta - y(i,N2)*xeta;
            D = (1+(gamma*g)/(deltat*delt*lambda));
            capC = (-2*alpha*g-2*gamma*g+2);
            B1 = ((1/4)*(3*p1(i+1,N2) - 4*p1(i+1,N2-1) + p1(i+1,N2-2) - 3*p1(N1-1,N2) + 4*p1(N1-1,N2-1) - p1(N1-1,N2-2)));
            R1 = p0(i,N2)/(deltat*delt*lambda);
            R2 = 2*lambda*kay*(p1(i+1,N2) - p1(N1-1,N2))/2;
            R3 = c*p1(i,N2)/(delt*lambda*r);
            R4 = p1(i,N2-1);
        
            %Actual boundary equation
            p2(i,N2) = (g*alpha/D)*p1(i+1,N2) + (g*alpha/D)*p1(N1-1,N2)-(2*g*beta/D)*B1+(gamma*g/D)*R1 - (gamma*g/D)*R2 - (gamma*g/D)*R3 ...
                + (gamma*g/D)*R4+(gamma*g/D)*p1(i,N2-1)-p0(i,N2)/D+(capC/D)*p1(i,N2);
        else
        	ypsi = (y(i+1,N2)-y(i-1,N2))/2;
            xpsi = (x(i+1,N2)-x(i-1,N2))/2;
            
            jacobian = xpsi*yeta-xeta*ypsi;
            g = ((c*deltat)/jacobian)^2;

            alpha = xeta^2 + yeta^2;
            beta = xpsi*xeta + ypsi*yeta;
            gamma = xpsi^2 + ypsi^2;
            
            %Variables used to ease calculation
            delt = c/(r*jacobian);
            lambda = y(i,N2)*xpsi-x(i,N2)*ypsi;
            kay = x(i,N2)*yeta - y(i,N2)*xeta;
            D = (1+(gamma*g)/(deltat*delt*lambda));
            capC = (-2*alpha*g-2*gamma*g+2);
            B1 = ((1/4)*(3*p1(i+1,N2) - 4*p1(i+1,N2-1) + p1(i+1,N2-2) - 3*p1(i-1,N2) + 4*p1(i-1,N2-1) - p1(i-1,N2-2)));
            R1 = p0(i,N2)/(deltat*delt*lambda);
            R2 = 2*lambda*kay*(p1(i+1,N2) - p1(i-1,N2))/2;
            R3 = c*p1(i,N2)/(delt*lambda*r);
            R4 = p1(i,N2-1);
          
            %Actual boundary equation
            p2(i,N2) = (g*alpha/D)*p1(i+1,N2) + (g*alpha/D)*p1(i-1,N2)-((2*g*beta)/D)*B1+(gamma*g/D)*R1 - (gamma*g/D)*R2 - (gamma*g/D)*R3 ...
                + (gamma*g/D)*R4+(gamma*g/D)*p1(i,N2-1)-p0(i,N2)/D+(capC/D)*p1(i,N2);
         end
    end
        %wrapping up
        p2(N1,N2) = p2(1,N2);
    

    %Stop criteria
    
    fullvalue = 0;

    for i = 1:N1
        for j = 1:N2
            diff = abs(abs(p2(i,j))-abs(p1(i,j)));
            fullvalue = fullvalue + diff;
        end
    end
    
    average = fullvalue / (N1*N2);
    
    if time <=2000
        if(average < firsttol)
            disp(time)
            break
        end
    elseif time > 2000
        if(average < secondtol)
            disp(time)
            break
        end
    end
    
    
    
    %Update pressure matrix
    p0 = p1;
    p1 = p2;
    
    
    if mod(time,100) == 0
        disp(time)
    end
end

%creation of incident wave for entire domain, to graph
pinc = zeros(N1,N2);
for new = 1:time
    for i = 1:N1
        for j = 1:N2
            pinc(i,j) = exp(sqrt(-1)*k*x(i,j))*exp(-sqrt(-1)*omega*new*deltat);
        end
    end
end



ptotal = p2 + pinc;