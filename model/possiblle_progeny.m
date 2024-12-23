function cluster_contri=possiblle_progeny (initial_size,f,p)
   p_div=f*(1-p);
   for (i=0:initial_size)
       cluster_contri(initial_size+i)=nchoosek(initial_size,i)*p_div^i*(1-p_div)^(initial_size-i);
   end
end