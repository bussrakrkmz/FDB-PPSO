function [] = fdb_ppso() 
 
    [npop, nvar, maxit, xmin, xmax] = problem_terminate();    
   
    dx=xmax-xmin;               
    vmax=0.5*dx;
   
    it=1; %Function of Evaluation Counter
     
    %Main Looop
    while( it<=maxit) 
        
        % Initialization
        if it==1
            
            gbestcost(1)=inf;  
            
            for i=1:npop
                
                velocity(i,:)=zeros(1,nvar);  
                
                delta(i)=unifrnd(0,2*pi); 
                               
                position(i,:)=xmin +(xmax-xmin)*rand(nvar,1);
    
                cost(i)=problem( position(i,:)');

                pbest(i,:)=position(i,:); 
                pbestcost(i)=cost(i);
           
                if pbestcost(i)<gbestcost(it)
                    gbest=pbest(i,1:nvar); 
                    gbestcost(it)=pbestcost(i);
                end
            end
        else
            gbestcost(it)=gbestcost(it-1);

            for i=1:npop

                aa=2*(sin(delta(i)));
                bb=2*(cos(delta(i)));
                ee=abs(cos(delta(i)))^aa;
                tt=abs(sin(delta(i)))^bb;

                velocity(i,:)=(ee)*(pbest(i,:)-position(i,:)) +(tt)*(gbest-position(i,:));

                velocity(i,:)=min(max(velocity(i,:),-vmax),vmax);

                position(i,:)=position(i,:)+velocity(i,:);
                
                position(i,:)=min(max(position(i,:),xmin),xmax);
                
                cost(i)=problem( position(i,:)' );
                
                delta(i)=delta(i)+(abs(aa+bb)*(2*pi));

                vmax=(abs(cos(delta(i)))^2)*dx;
               
                if cost(i)<pbestcost(i)
                    pbest(i,:)=position(i,:);
                    pbestcost(i)=cost(i);
                    
                    if pbestcost(i)<gbestcost(it)
                        gbest=pbest(i,:);
                        gbestcost(it)=pbestcost(i);
                    end
                end
            end
        end
        it=it+1;
    end       
    Cost_Result=gbestcost(end);
    
    fprintf('Best Fitness: %d\n', Cost_Result);
    disp('Best Solution:'); 
    disp(gbest);
end
