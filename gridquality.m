function gridquality(x,y,N1,N2)

MDO = 0;
NDO = 0;
minJac = 1000;
jacobian = 1;
disp('Check of Grid Quality')
for i = 2:N1-1
    for j = 2:N2-1
        xeta = (x(i,j+1) - x(i,j-1))/2;
        yeta = (y(i,j+1) - y(i,j-1))/2;

        ypsi = (y(i+1,j)-y(i-1,j))/2;
        xpsi = (x(i+1,j)-x(i-1,j))/2;
        
        alpha = xeta^2 + yeta^2;
        beta = xpsi*xeta + ypsi*yeta;
        gamma = xpsi^2 + ypsi^2;
        angle = acos(beta/((gamma*alpha)^(1/2)));
        
        
        if(jacobian < 0)
            val = 1;
        elseif(jacobian > 0)
            val = 0;               
        end
        
        jacobian = xpsi*yeta-xeta*ypsi;
        if i ~= 2 && j~= 2
            if(jacobian < 0 && val == 0)
                displaypoint = ['jacobian zero at', '(',num2str(i),',',num2str(j),')'];
                disp(displaypoint);
            elseif(jacobian > 0 && val == 1)
                displaypoint = ['jacobian zero at', '(',num2str(i),',',num2str(j),')'];
                disp(displaypoint);
            end
        end
        
        if(abs(jacobian) < minJac)
            minJac = abs(jacobian);
        end
        
        
        if(abs(90-angle) > MDO)
            MDO = abs(90-angle);
        end
        
        NDO = NDO + abs(90-angle);
        
    end
end

NDO = NDO/(N1-2);
NDO = NDO/(N2-2);


gridsizetodisp = [num2str(N1), 'X', num2str(N2)];
disp(gridsizetodisp)
ndotodisp = ['NDO: ', num2str(NDO)];
disp(ndotodisp)
mdotodisp = ['MDO: ', num2str(MDO)];
disp(mdotodisp)
absjacobitodisp = ['Absolute minimum value of the Jacobian: ', num2str(minJac)];
disp(absjacobitodisp)

