function V3 = L21norm(V3_tmp, tao) 
for i = 1:size(V3_tmp,2)
   V3(:,i) = softVector(V3_tmp(:,i),tao); 
end
end


function x = softVector(x,tao)
tmp = max(norm(x)-tao,0);
x = tmp/(tmp+tao)*x;
end