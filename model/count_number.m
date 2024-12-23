function [cell_no,dividing_frac]=count_number(param,T)
    if(T==0) T=0.0001; end
    f_1=param(1);
    T_1=param(2);
    N_cyc_1=round(param(3));
    T_transit=param(4);
    f_2=param(5);
    T_2=param(6);
    N_cyc_2=round(param(7));
    p_1=param(8);
    p_2=param(9);
    cell_no=0;
    dividing_cells_total=0;
    T_wv1=min([T, T_transit]);
    if (N_cyc_1)
    for (N=1:N_cyc_1)
        T_c=T_wv1/N;
        T_c_2=T_wv1/(N-1);
        prob=expcdf(T_c_2,T_1)-expcdf(T_c,T_1);
        [total_cell,dividing_cells]=count_number_cycles(N-1,f_1,p_1);
        cell_no=cell_no+prob*total_cell;
        dividing_cells_total=dividing_cells_total+prob*dividing_cells;
    end 
        prob=expcdf(T_c,T_1);
        [total_cell,dividing_cells]=count_number_cycles(N,f_1,p_1);
        cell_no=cell_no+prob*total_cell;
    else
     cell_no=1;
    end 

    if (T>T_transit && N_cyc_2)
        cell_no_initial=cell_no;
        cell_no=0;    
        for (N=1:N_cyc_2)
            T_c=(T-T_transit)/N;
            T_c_2=(T-T_transit)/(N-1);
            prob=expcdf(T_c_2,T_2)-expcdf(T_c,T_2);
            [total_cell,dividing_cells]=count_number_cycles(N-1,f_2,p_2);
            cell_no=cell_no+prob*total_cell*cell_no_initial;
            dividing_cells_total=dividing_cells_total+prob*dividing_cells;
        end 
        prob=expcdf(T_c,T_2);
        [total_cell,dividing_cells]=count_number_cycles(N,f_2,p_2);
        cell_no=cell_no+prob*total_cell*cell_no_initial;
    end

    dividing_frac=dividing_cells_total./cell_no;
    
end

function [total_cell,dividing_cells]=count_number_cycles(comp_cyc,f,p)
   total_cell=1;
   for (N=1:comp_cyc)
       total_cell=total_cell*(1+f*(1-p));
   end
   dividing_cells=total_cell*f;
end

