
% Amplitude Scattered field Solution and Amplitude Total field Solution 
% evaluated at grid points x(\xi,\eta) y(\xi,\eta)


function ExactvsNumerical(x,y,p2,k,r0,theta,n1,n2)

for j=1:n2
    for i=1:n1
        r(i,j)=sqrt(x(i,j)^2+y(i,j)^2);
        th(i,j)=subangleN(x(i,j),y(i,j));
    end
end
% th(1,n2) = 2*pi;
% th(n1,n2) = 0;
sum=0;
for j=1:n2    
    for i=1:n1
        m=1;
        sum=-besselj(0,k*r0)*besselh(0,k*r(i,j))/besselh(0,k*r0); %first term of actual solution
        while m<21
            nextterm=-2*(sqrt(-1))^m*besselj(m,k*r0)*besselh(m,k*r(i,j))...
                       *cos(m*th(i,j))/besselh(m,k*r0);
            sum=sum+nextterm;
            m=m+1;
        end
        pinc=exp(sqrt(-1)*(k*r(i,j)*cos(th(i,j)-theta))); %incident pressure
        pscactual(i,j)=abs(sum);
        actual(i,j)=abs(sum+pinc); %steady state
    end
end

figure;
surf(x,y,actual);
shading interp;
title('Exact Solution');
view(-4,86);

% Computing errors

%---------This part of the code finds the max and average norm error 
%                         for the steady state
er=0;
totalerr=0;
for i=1:n1
    for j=1:n2
        err(i,j)=abs(pscactual(i,j) - abs(p2(i,j)));
        if err(i,j)>er
            er=err(i,j);
            indexi=i;
            indexj=j;
        end
        totalerr=totalerr+err(i,j);
    end
end
avgerr=totalerr/(n1*n2);   
fprintf('Maximum error:%d\n ',er);
fprintf('Max error occurs on ray:%i and ring:%i\n',indexi,indexj);
fprintf('Which is the point:%d %d\n',x(indexi,indexj),y(indexi,indexj));
disp('Average Maximum Norm error for p_sc');
avgerr  

figure;
surf(x,y,err);
title('Error graph at each point');
view(-4,86);


%********************L2 Error at Artificial Boundary**********************
fprintf('\n*****************Errors at Artificial Boundary r=r(N)***************************\n');
L2error = 0;
L2pscactual = 0;
dtheta = 2*pi/(n1-1);
for i=1:n1
%     L2error = L2error + diff(N,j)^2*r(N);
    L2error = L2error + err(i,n2)^2;
    L2pscactual = L2pscactual + (abs(pscactual(i,n2)))^2;
end

L2error = sqrt(L2error*dtheta);
L2pscactual = sqrt(L2pscactual*dtheta);
RelL2error = L2error/L2pscactual;

fprintf('\nL2 Error for u_sc on r = R: %d\n', L2error);
fprintf('\nRelative L2 Error for u_sc on r = R: %d\n', RelL2error);


figure;
plot(th(1:n1,n2),pscactual(1:n1,n2),'-k',th(1:n1,n2),abs(p2(1:n1,n2)),'-r');
str={['Comparison at Artificial Boundary: Analytical vs Approximated  ']...
    ,[num2str(n1) 'x'  num2str(n2)   '    k= ' num2str(k) ]};
xlabel('theta');
ylabel('Field at Artificial Boundary');
legend('Analytical', 'Numerical');
title(str);

figure;
polar(th(1:n1,n2),pscactual(1:n1,n2),'-k'); hold on;
polar(th(1:n1,n2),abs(p2(1:n1,n2)),'-r'); hold off;
str={['Comparison at Artificial Boundary: Analytical vs Approximated  ']...
    ,[num2str(n1) 'x'  num2str(n2)   '    k= ' num2str(k) ]};
xlabel('theta');
ylabel('Field Artificial Boundary');
legend('Analytical', 'Numerical');
title(str);


