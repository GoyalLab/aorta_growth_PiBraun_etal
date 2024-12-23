clear all;
%addpath ..a
load total_and_dividing_cells.mat;
trial_max=1000;
guess_param=zeros(trial_max,9);
optim_param_list=zeros(trial_max,9);
error_list=zeros(trial_max,9);
parfor (trial=1:trial_max)
    guess_param(trial,:)=[0.5*rand() 10*rand() 0.5+3.5*rand() 10+15*rand() 0.5*rand() 10*rand() 0.5+3.5*rand() rand() rand()];
   [optim_param,optim_error]=run_optimization(guess_param(trial,:),T_list,N,F);
    optim_param_list(trial,:)=optim_param;
    error_list(trial)=optim_error;
    trial
end
 save ('results_pruning_expo.mat','guess_param','optim_param_list','error_list');
 %%
function [optim_param,optim_error]=run_optimization(old_param,T_list,N,F)
    old_error=find_error (old_param,T_list,N,F);
    step_s=0.1;
    new_param=zeros(size(old_param));
    while (step_s>0.0009)
    for (step=1:10000)
        for (i=1:length(old_param))
            new_param(i)=old_param(i)*(1+step_s*randn());
	        if (new_param(4)>30) new_param(4)=30; end
            if (new_param(1)<0) new_param(1)=0.0001; end
            if (new_param(5)<0) new_param(5)=0.0001; end
            if (new_param(3)<0.5) new_param(3)=0.5001; end
	        if (new_param(7)<0.5) new_param(7)=0.5001; end 
            if (new_param(8)<0) new_param(8)=0.0001; end  
            if (new_param(8)>1) new_param(8)=0.9999; end 
            if (new_param(9)<0) new_param(9)=0.0001; end  
            if (new_param(9)>1) new_param(9)=0.9999; end 
      
         end
            new_error=find_error (new_param,T_list,N,F);
        if (new_error<old_error)
                old_param=new_param;
                old_error=new_error;
        end
    end
         step_s=step_s*0.1;
    end
    optim_param=old_param;
    optim_error=old_error;

end
%%
% for (i=1:length(T_list))
%         [cell_no,dividing_frac]=count_number(old_param,T_list(i));
%         cell_no_list(i)=cell_no;
%         dividing_list(i)=dividing_frac;
% end
% cell_no_normalized=cell_no_list(2:end)/cell_no_list(1);
% figure (1);
% plot (T_list(2:end),N(:,1),'color','k','Linewidth',2);
% hold on
% plot (T_list(2:end),cell_no_list(2:end)/cell_no_list(1),'color','b','Linewidth',2);
% hold off
% figure (2);
% plot (T_list,F(:,1),'color','k','Linewidth',2);
% hold on
% plot (T_list,dividing_list,'color','b','Linewidth',2);
% hold off;
% %count_number(param,60)
%%



function error= find_error (param,T_list,N,F)
    cell_no_list=zeros(1,length(T_list));
    dividing_list=zeros(1,length(T_list));
    for (i=1:length(T_list))
        [cell_no,dividing_frac]=count_number(param,T_list(i));
        if (size(cell_no))cell_no_list(i)=cell_no(1);end
	if (size(cell_no))dividing_list(i)=dividing_frac(1);end
    end
       cell_no_normalized=cell_no_list(2:end)/cell_no_list(1);
       error_list=([(cell_no_normalized' -N(:,1))./N(:,2); (dividing_list'-F(:,1))./F(:,2)]);
       error=norm(error_list);
end
