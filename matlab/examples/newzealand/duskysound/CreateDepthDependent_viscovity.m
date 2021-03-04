function out=CreateDepthDepedent_viscovity(shearzonemodel,etaM)
    fid  = fopen(shearzonemodel,'r');
    data = textscan(fid,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f','commentstyle','#');
    fclose(fid);
    nshz = length(data{1});
    leng = data{5};width = data{6};thickness = data{7};
    lo   = data{10};wo   = data{11};to       = data{12};
    ndepth = sum(thickness./to);
    nd0 = 0;
    out = [];
    snc = 0;
    for i = 1:length(data{1})
         m   = data{5}(i)/data{10}(i);  
         n   = data{6}(i)/data{11}(i);
         nl  = data{7}(i)/data{12}(i);
         for k = 1:n
             nd = nd0;
             for j =1:nl
                 nd = nd +1;
                 sn0 = m*nl*(k-1) + (j-1)*m+snc;
                 t0 = 1 +sn0;
                 t1 = t0 + m -1;
                 if k ==1
%                     disp(sprintf('lvl = %d    n0 = %4d    n1 = %4d',i,t0,t1))
                 end
                 out(t0:t1) = etaM(nd);
             end
         end
         snc = t1;
         nd0 = nd;
    end
    out = out';
    
end 