function Lx=get_LaplacianMatrix(r,Qx)
Ax_bar=zeros(r,r);
Dx_bar=zeros(r,r);
for i =1:1:length(Qx)-1
    if(Qx(i)~=Qx(i+1))
        Ax_bar(Qx(i),Qx(i+1))=1; %°ëÕý¶¨¾ØÕó
        Ax_bar(Qx(i+1),Qx(i))=1; 
    end
end
for j=1:1:r
   Dx_bar(j,j)=sum(Ax_bar(j,:)); 
end
Lx=Dx_bar-Ax_bar;