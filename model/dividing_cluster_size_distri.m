function csd=dividing_cluster_size_distri(param,T)
        f_1=param(1);
        T_1=param(2);
        N_cyc_1=round(param(3));
        p=param(8);
        max_rounds=4;
        prob_rounds=zeros(max ([max_rounds, N_cyc_1])+1,1);
        for(N=1:N_cyc_1)
            T_c=T/N;
            T_c_2=T/(N-1);
            prob_rounds(N)=expcdf(T_c_2,T_1)-expcdf(T_c,T_1);
        end
            prob_rounds(N+1)=expcdf(T_c,T_1);
            cluster_size_cycles=zeros(2^max_rounds,max_rounds+1);
            cluster_size_cycles(1,1)=1;
            for (N=1:max_rounds)
                if(N==1)
                    f=1;
                else
                    f=f_1;
                end
            for (k=1:2^(N-1))
                cluster_contri=possiblle_progeny(k,f,p);
                cluster_size_cycles(1:2*k,N+1)=cluster_size_cycles(1:2*k,N+1)+cluster_size_cycles(k,N)*cluster_contri';             
            end        
            end
            csd=zeros(2^max_rounds,1);
            for (N=1:max_rounds+1)
                csd=csd+prob_rounds(N)*cluster_size_cycles(:,N);
            end
end