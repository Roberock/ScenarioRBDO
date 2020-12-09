function norm_data=Normalizein01(bla)
% normalize in [0,1]
 norm_data = (bla - min(bla)) ./ ( max(bla) - min(bla) ); 
end