function f = mse(G,Q,p,v,t)
    m = length(Q);
    f = 0;
    for i=1:m
        p_i = G.exp(p,v,t(i));
        q_i = Q{i};
        f = f + grarc(q_i,p_i);
        %fprintf('d(p_%d,q_%d)=%.10f\n', i,i,cdc_grDistSq(q_i,p_i));
    end
    f=f/m;
end