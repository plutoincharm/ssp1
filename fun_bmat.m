function B = fun_bmat(x,u,dt)
%%% For finding B matrix 
global A B Nx Nu pert MI L m nx ny tx ty g  ra    
for i = 1:Nx

        for j = 1:(Nu)
            xnew = x + fun_xdot(x,u,dt)*dt/2; % is this constant ??
            %xplus = x;
            %xplus(j) =  x(j) + pert; %%% pert if iside diffferentiations
            uplus = u;
            uplus(j) =  u(j) + pert;
       
            if i<= Nx/2

          
                    %%% Jacobian of  Velocity with  u
                    f = xnew((Nx/2)+i);
                    F = x(i) + f*dt;
                    qddot = fun_qddot(x,uplus,dt);
                    fplus = x((Nx/2)+i) + qddot(i)*dt/2;  
                    Fplus = x(i) + fplus*dt; 
                    B(i,j) = (Fplus - F)/pert;


           
            elseif i>Nx/2
                %%%  Jacobian of acceleration  with  q  and qdot
                       %{
                         if j<=Nx/2
                                %%% pert in q for jac with acceleration
                                xplusr =    xnew;
                                xplusr(j) = xnew(j) + pert; 
                                
                        elseif j > Nx/2        
                
                            %%% pert in qdot for jac with acceleration
                            xplusr = xnew;
                            xplusr(j) =  xnew(j) + pert; 
                            xplusr(j-(NX/2)) =  xnew(j-(NX/2)) + pert*dt/2;
                        end
                       %}
                
                fr = fun_qddot(xnew,u,dt);
                Fr = x(i) + fr(i-(Nx/2))*dt;
                fplusr = fun_qddot(xnew,uplus,dt); 
                Fplusr = x(i) + fplusr(i-(Nx/2))*dt; 
                B(i,j) = (Fplusr - Fr)/pert;
            end    
                
        end  

   

end